#%% [packages]
# all packages are imported initially
import pandas as pd
import numpy as np
from Bio.PDB import *
from bs4 import BeautifulSoup
import urllib
import requests
import urllib.request, os, requests as req, pickle, time
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import operator
from Bio.Align.Applications import ClustalOmegaCommandline
import os.path
from os import listdir
from os.path import isfile, join
import io
from scipy import stats
import subprocess
from numbers import Number
from scipy.stats import hypergeom
import re

#%% [pre-processing]
# all phosphorylatino sites are read from dbPAF database
os.chdir("file location")

all_data=pd.read_excel("41587_2019_344_MOESM4_ESM.xlsx", sheet_name="annotated_phosphoproteome", engine="openpyxl")

functional_scores=pd.read_excel("41587_2019_344_MOESM5_ESM.xlsx", engine="openpyxl")

all_data=all_data.merge(functional_scores, on=["uniprot","position"],how="left")

#add sequences
#read fasta file
fasta_file="sequences.fasta"

with open(fasta_file,"r") as fastafile:

    seqs={}

    for line in fastafile:
        line=line.strip()
        if line.startswith(">"):
            name=line[1:]
            seqs[name]=""
        else:
             seqs[name] = seqs[name] + line

fasta=pd.DataFrame.from_dict(seqs,orient="index")
fasta.columns=["sequence"]
fasta.insert(0,"ids",fasta.index)
fasta.reset_index(inplace=True,drop=True)

#seperate ids part to get uniprot ids
fasta_divided=pd.DataFrame(fasta["ids"].str.split("|").tolist(),columns=("indicator","Protein Group","description"))

# match them with sequences
fasta=fasta.join(fasta_divided,how="left")
fasta=fasta.drop(["ids","indicator"],axis=1)
#create a small df to merge with sequence and add it to the main data
temp=all_data.copy()

temp=temp.merge(fasta,right_on="Protein Group",left_on="uniprot",how="left")

all_data["Sequence"]=list(temp["sequence"])

all_data.dropna(subset=["Sequence"], inplace=True)
# phospho windows (13 size) are created from total sequences based on their position
phospho_windows=[]
#loop throughout all phosphorylation
for i,row in all_data.iterrows():
    #location of phospho-site
    place=int(row["position"])
    # if it is in beginning or end of the protein, fill remaining resiudes with "_"
    if place<7:
        phospho_windows.append((7-place)*"_"+row["Sequence"][0:place+6])
    elif (len(row["Sequence"])-place)<7:
        phospho_windows.append(row["Sequence"][place-7::]+(6-len(row["Sequence"])+place)*"_")
    else:
        phospho_windows.append(row["Sequence"][place-7:place+6])

all_data["phospho window"]=phospho_windows

# eliminate phosphorylation sites with wrong sequence
all_data=all_data.loc[all_data["phospho window"].str[6].isin(["S","T","Y"])]


#%% [obtaining structural information]
# PDB results are found by webscrapping of uniprot due to unavailability of coverage and resolution information from PDB RCSB and uniprot bulk
pdb_and_coverage=[]
count=1
test=0
# loop through all phospho_sites
for uniprot in all_data["uniprot"]:
    #connecting to the xml page of each protein and connect to the page
    url = "http://www.uniprot.org/uniprot/"+uniprot+".xml"
    print(url)
    session = requests.Session()
    retry = Retry(connect=4, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)

    response = ""

    while response == "":
        # some pages are deleted or there are some connection errors
        try:
            #to read xml files beautifulsoup is used with searching motifs in PDB information part
            response = session.get(url, timeout=5)
            thedata = response.content
            soup=BeautifulSoup(thedata, "xml")
            #information part of uniprot page
            references=soup.find_all("dbReference")
            # PDB informations
            references_2=[str(x) for x in references if "type=\"PDB\"" in str(x) and "value=" in str(x)]
            temp=[]
            #obtain all PDB informations, id,chainsi resolution and coverage
            for ref in references_2:
                if ("id=\"" in str(ref)) & ("\"chains\" value=\"" in str(ref)) & ("resolution\" value=" in str(ref)):
                    temp.append([str(ref)[str(ref).find("id=\"")+4:str(ref).find("\" type")],str(ref)[str(ref).find("\"chains\" value=\"")+16:str(ref).find("\"/>\n</dbReference>")],float(str(ref)[str(ref).find("resolution\" value=\"")+19:str(ref).find("\"/>\n<property type=\"chains\"")])])
                else:
                    pass
            pdb_and_coverage.append(temp)
        except urllib.error.HTTPError:
            # if no page is found
            pdb_and_coverage.append([])
            response="*"
        except:
            # uniprot cut connection once in a while, so reconnect after a while
            print("ZZzzzz...")
            print(url)
            time.sleep(5)
            test=test+1
            continue
    print(str(count) + " is done")
    count=count+1


all_data["PDB_information"]=pdb_and_coverage

#export phosphorylation sites with PDB information due to time-consuming webscrapping
all_data_2=all_data.copy()
pdb_and_coverage_2=[]
for i in pdb_and_coverage:
    temp=[]
    for x in i:
        temp.append(":".join(list(map(str,x))))
    pdb_and_coverage_2.append(";".join(temp))

all_data_2["PDB_information"]=pdb_and_coverage_2

all_data_2.to_csv("all_phospho_data_with_pdb_id.csv")
all_data_2.to_excel("all_phospho_data_with_pdb_id.xlsx")

"""
#%% [optional read again from saved dataset]
PDB_information_2=pd.read_excel("all_phospho_data_with_pdb_id.xlsx")
PDB_information_2=list(PDB_information_2["PDB_information"])

new=[x.split(";") if type(x)==str else x for x in PDB_information_2 ]
new = [[s.split(":") if type(s)==str else s for s in l] if type(l)==list else l for l in new]

all_data["PDB_information"]=new
"""

#%% [filtering PDB information with best resolution and highest coverage]
phospho_check=[]
#iterate thorugh each row
for i,row in all_data.iterrows():
    #check empty ones and store them as nan
    if type(row["PDB_information"])==float:
        phospho_check.append([])
    # if there are pdb files store best ones
    else:
        # create temp file to save shape
        temp=[]
        #save real places of phospho sites if it contains more than one save it as int list
        phospho_site=row["position"]
        # loop through each pdb file
        for pdbs in row["PDB_information"]:
            #check if pdb have more than one place
            if "," in pdbs[1]:
                #then split them
                for each in pdbs[1].split(", "):
                    #use it if no match occurs
                    temp2=3
                    #store numbers
                    numbs=list(map(int,each.split("=")[1].split("-")))
                    if (numbs[0]<np.array(phospho_site)) & (np.array(phospho_site)<numbs[1]):
                        temp2=pdbs
                    else:
                        pass
                if temp2==3:
                    pass
                else:
                    temp.append(temp2)
            #if there is one interval, do the same
            else:
                numbs=list(map(int,pdbs[1].split("=")[1].split("-")))
                if (numbs[0]<np.array(phospho_site)) & (np.array(phospho_site)<numbs[1]):
                    temp.append(pdbs)
                else:
                    pass
        #find pdb with best resolution
        if temp==[]:
            pass
        else:
            x=[float(i[2]) for i in temp]
            temp=temp[x.index(min(x))]

        phospho_check.append(temp)

all_data["PDB_filtered"]=phospho_check

list_of_pdb_filtered=list(set([a[0] for a in all_data["PDB_filtered"] if len(a)>0]))

with open("pdb_files.txt", 'w') as output:
    for row in list_of_pdb_filtered:
        output.write(str(row) + '\n')

#%% [model PDBs from modbase]

os.mkdir(r"location\Modbase")
os.chdir(r"location\Modbase")

#model list are taken
empty_proteins=list(set(all_data.loc[all_data["PDB_filtered"].map(lambda d: len(d)) == 0,"uniprot"]))
#downloaded from modbase
for protein in empty_proteins:
    os.system("curl -L \"https://salilab.org/modbase/retrieve?databaseID={}\" -o {}.xml".format(protein,protein))


#%%  [filtration of models based on their modpipe quality score]

mypath=r"location\Modbase"
model_pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#read all models and select with eligible ones and store their names and chain
pdb_model_names={}
model_dict={}
for xml in model_pdbfiles:
    with open(xml, "r") as file:
        # Read each line in the file, readlines() returns a list of lines
        content = file.readlines()
        # Combine the lines in the list into a string
        content = "".join(content)
    #save them in dictionary by their uniprot name and read as xml
    model_dict[xml.split(".")[0]]=content
    bs_content = BeautifulSoup(content, "lxml")
    split_model=bs_content.find_all("pdbfile")
    #split each model
    split_model_2=[str(x) for x in split_model if "REMARK 220" in str(x)]
    quality=[]
    names=[]
    location=[]
    #filter models based on score and obtain all information
    for models in split_model_2:
        for line in models.splitlines():
            if "ModPipe Quality Score" in line:
                if float(line.replace(" ","").split(":")[1])>=1.1:
                    quality.append(models)
                    names.append(re.findall("TEMPLATE PDB.*$",models,re.MULTILINE)[0].replace(" ","").split(":")[1]\
                                 +":"+re.findall("TEMPLATE CHAIN.*$",models,re.MULTILINE)[0].replace(" ","").split(":")[1]
                                 +"="+re.findall("TEMPLATE BEGIN.*$",models,re.MULTILINE)[0].replace(" ","").split(":")[1]\
                                 +":"+re.findall("TEMPLATE END.*$",models,re.MULTILINE)[0].replace(" ","").split(":")[1])
                else:
                    pass
    model_dict[xml.split(".")[0]]=quality
    pdb_model_names[xml.split(".")[0]]=names
    del model_dict[xml.split(".")[0]]
    print(xml)

#make names upper and if it is duplicates eliminate them
for protein in pdb_model_names.keys():
    pdb_model_names[protein]=list(set([x.upper() for x in pdb_model_names[protein]]))
# transform to the dataframe
pdb_model_names=pd.DataFrame(list(pdb_model_names.items()),columns=["uniprot","Model PDB"])

# add to the main dataframe
all_data=all_data.merge(pdb_model_names, on="uniprot", how="left")


#%% [filtering model PDB information with best resolution and highest coverage]

phospho_check_2=[]
#iterate thorugh each row
for i,row in all_data.iterrows():
    #check empty ones and store them as nan
    if type(row["Model PDB"])==float:
        phospho_check_2.append([])
    # if there are pdb files store best ones
    else:
        # create temp file to save shape
        temp=[]
        #save real places of phospho sites if it contains more than one save it as int list
        phospho_site=row["position"]
        # loop through each pdb file
        for pdbs in row["Model PDB"]:
            #check if pdb have more than one place
            if "," in pdbs[1]:
                #then split them
                for each in pdbs[1].split(", "):
                    #use it if no match occurs
                    temp2=3
                    #store numbers
                    numbs=list(map(int,each.split("=")[1].split(":")))
                    if (numbs[0]<np.array(phospho_site)) & (np.array(phospho_site)<numbs[1]):
                        temp2=pdbs
                    else:
                        pass
                if temp2==3:
                    pass
                else:
                    temp.append(temp2)
            #if there is one interval, do the same
            else:
                try:
                    numbs=list(map(int,pdbs.split("=")[1].split(":")))
                    if (numbs[0]<np.array(phospho_site)) & (np.array(phospho_site)<numbs[1]):
                        temp.append(pdbs)
                    else:
                        pass
                except:
                    pass
        if temp==[]:
            phospho_check_2.append(temp)
        else:
            phospho_check_2.append(temp[0])

#if it is empty replace it with nan
all_data["model_PDB_filtered"]=phospho_check_2

list_of_model_pdb_filtered=list(set([a.split("=")[0].split(":")[0] for a in all_data["model_PDB_filtered"] if len(a)>0]))

with open("model_pdb_files.txt", 'w') as output:
    for row in list_of_model_pdb_filtered:
        output.write(str(row) + '\n')

#%%

#download all pdb files PDBdownload.py
os.chdir("location")
import PDBdownload


#%% [transform PDB files to monomer state which phospho_site locates]

#get all the pdb files
mypath="location\\PDB"
lesspdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
lesspdbfiles=[x for x in lesspdbfiles if x.endswith(".ent")]
#filter filtered list for empty lists
filtered_2=[x for x in phospho_check if x!=[]]

os.chdir("location\\PDB")

#create files only with desired chain
for less in lesspdbfiles:
    #choose the chain and create pdb file contains only chain atoms
    chain=[t for t in filtered_2 if (t[0]==less.split(".")[0].split("-")[0].upper())]
    for alp in chain:
        with open(less,"r") as f:
            with open(less.split(".")[0]+"_"+alp[1].split("=")[0].split("/")[0]+"_less.ent","w") as f1:
                test=f.readline()
                while test:
                    if test.split()[0]=="ATOM":
                        if test.split()[4][0]==alp[1].split("=")[0].split("/")[0]:
                            f1.write(test)
                        else:
                            pass
                    else:
                        f1.write(test)
                    test=f.readline()
    print(less)


all_data.reset_index(drop=True,inplace=True)

#%%

#download all model pdb files PDBdownload.py
os.chdir("location")
import PDBdownload

#%% [transform model PDB files to monomer state which phospho_site locates]

#get all the pdb files
mypath="location\\Modbase models"
lesspdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
lesspdbfiles=[x for x in lesspdbfiles if x.endswith(".ent")]
#filter filtered list for empty lists
filtered_2=[x for x in phospho_check_2 if x!=[]]

os.chdir("location\\Modbase models")

#create files only with desired chain
for less in lesspdbfiles:
    #choose the chain and create pdb file contains only chain atoms
    chain=[t for t in filtered_2 if (t.split(":")[0]==less.split(".")[0].split("-")[0].upper())]
    for alp in chain:
        with open(less,"r") as f:
            with open(less.split(".")[0]+"_"+alp.split(":")[1].split("=")[0].split("/")[0]+"_less.ent","w") as f1:
                test=f.readline()
                while test:
                    if test.split()[0]=="ATOM":
                        if (test.split()[4][0]==alp.split(":")[1].split("=")[0].split("/")[0]) | (alp.split(":")[1].split("=")[0].split("/")[0]==""):
                            f1.write(test)
                        else:
                            pass
                    else:
                        f1.write(test)
                    test=f.readline()
    print(less)

#%%
#run freesasa program

#%% [monomer RSA score acquisition]

#get just chain freesasa result
os.chdir(r"location\PDB")

#read all freesasa results
mypath=r"location\PDB"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]

#change name of columns to make it appliable for all datasets
data_for_PDB=all_data[["PDB_filtered","PDB_places","position","phospho window"]]
data_for_PDB.rename(columns={"position":"Real Places"},inplace=True)
data_for_PDB.rename(columns={"phospho window":"Sequence Window"},inplace=True)

#3 letter to 1 letter dictionary
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def RSA_score_acquisiton (pdb_type, pdbfiles):

    freesasa_results=[]
    # loop through all phospho_sites
    for o, row in data_for_PDB.iterrows():
        print(o)
        #if there is no PDB
        if row[pdb_type]==[]:
            freesasa_results.append("")
        else:
            x=row[pdb_type]:
            #save name of pdb
            pdb_name=x[0].lower()+"_"+x[1].split("=")[0].split("/")[0]+"_less"
            # no file
            if not "pdb{}.ent.out".format(pdb_name) in pdbfiles:
                freesasa_results.append("")
                continue
            else:
                #if there is a match with phospho site
                lines=[]
                #read the file
                for line in pd.read_csv("pdb{}.pdb.out".format(pdb_name),skiprows=[0,1,2,3,4,5,6], engine="python",chunksize=1, encoding="utf-8"):
                    lines.append(line.iloc[0,0])
                lines_2=[]
                #format location and residue of freesasa result without space seperation to fit on a dataframe
                for bla in lines:
                    if "REM" in bla:
                        lines_2.append(bla)
                    elif "RES" in bla:
                        blob=bla.split()
                        if len(blob[2])==5:
                            blob[2]=blob[2][0]+" "+blob[2][1::]
                            blob=" ".join(blob)
                            lines_2.append(blob)
                        else:
                            lines_2.append(bla)
                    else:
                        lines_2.append(bla)
                # fit to a dataframe
                pdb_dict=pd.read_csv(io.StringIO('\n'.join(lines_2)), delim_whitespace=True)
                # remove last rows to prevent columns disturbance
                count=0
                for a in pdb_dict.index:
                    if a==(('END', 'Absolute', 'sums')):
                        ind=count
                        break
                    else:
                        count=count+1
                        pass
                pdb_dict=pdb_dict[0:ind]
                # reset index to create 3 columns with position and residue information
                resetted=pdb_dict.reset_index()
                seq=resetted["level_1"].tolist()
                # due to the difference in sequences of freesasa result and protein,
                # created a sequence from freesasa result to find exact location
                seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
                found=0
                #search for phospho_site with 7 residue window when found, obtain RSA value
                for i in range(0,7):
                    place=seqs.find(row["Sequence Window"][i:i+7])
                    if place==-1:
                        continue
                    else:
                        found=place+6-i
                        break
                if found==0:
                    freesasa_results.append("")
                else:
                    freesasa_results.append(resetted.iloc[found]["REL"])
    return freesasa_results

all_data["freesasa result"]=RSA_score_acquisiton ("PDB_filtered", pdbfiles)


#%% [model monomer RSA score acquisition]

#get just chain freesasa result
os.chdir(r"location\Modbase models")

#read all freesasa results
mypath=r"location\Modbase models"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]

all_data["model freesasa result"]=RSA_score_acquisiton ("model_PDB_filtered", pdbfiles)


#%% [whole protein RSA score acquisition]

#get just chain freesasa result
os.chdir(r"location\PDB")

mypath=r"location\PDB"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if not "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]

all_data["freesasa intact result"]=RSA_score_acquisiton ("PDB_filtered", pdbfiles)

#%% [model whole protein RSA score acquisition]

#get just chain freesasa result
os.chdir(r"location\Modbase models")

#read all freesasa results
mypath=r"location\Modbase models"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if not "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]

all_data["model freesasa intact result"]=RSA_score_acquisiton ("model_PDB_filtered", pdbfiles)

#%% [categorization based on RSA]

def position_categorization(RSA):
    y=[]
    posit=[]
    for i,row in all_data.iterrows():
        if row[RSA]=="":
            posit.append("no_data")
        elif float(row[RSA])<=5:
            posit.append("core")
            y.append(i)
        elif (50>=float(row[RSA])>5):
            posit.append("intermediate")
            y.append(i)
        elif float(row[RSA])>50:
            posit.append("surface")
            y.append(i)
    return posit
all_data["position Just Chain"]=position_categorization("freesasa result")
all_data["position Complex"]=position_categorization("freesasa intact result")
all_data["Model position Just Chain"]=position_categorization("model freesasa result")
all_data["Model position Complex"]=position_categorization("model freesasa intact result")

#%% [final location detection]

#label Position according to just chain and complex

all_data.loc[((all_data["Position Just Chain"]=="surface") | (all_data["Position Just Chain"]=="intermediate"))&\
               (all_data["Position Complex"]=="core"), "Phos Position"]="Interface"

all_data.loc[(all_data["Position Just Chain"]=="core") & (all_data["Position Complex"]=="core"), "Phos Position"]="core"

all_data.loc[(all_data["Position Just Chain"]=="intermediate") & (all_data["Position Complex"]=="intermediate"), "Phos Position"]="intermediate"
all_data.loc[(all_data["Position Just Chain"]=="surface") & (all_data["Position Complex"]=="surface"), "Phos Position"]="surface"
all_data.loc[(all_data["Position Just Chain"]=="intermediate") & (all_data["Position Complex"]=="surface"), "Phos Position"]="surface"
all_data.loc[(all_data["Position Just Chain"]=="surface") & (all_data["Position Complex"]=="intermediate"), "Phos Position"]="surface"

#label Model Position according to just chain and complex
all_data.loc[(all_data["Model Position Just Chain"]=="core") & (all_data["Model Position Complex"]=="core"), "Model Position"]="core"
all_data.loc[(all_data["Model Position Just Chain"]=="intermediate") & (all_data["Model Position Complex"]=="intermediate"), "Model Position"]="intermediate"
all_data.loc[(all_data["Model Position Just Chain"]=="surface") & (all_data["Model Position Complex"]=="surface"), "Model Position"]="surface"
all_data.loc[((all_data["Model Position Just Chain"]=="surface") | (all_data["Model Position Just Chain"]=="intermediate"))&\
               (all_data["Model Position Complex"]=="core"), "Model Position"]="Interface"
all_data.loc[(all_data["Model Position Just Chain"]=="intermediate") & (all_data["Model Position Complex"]=="surface"), "Model Position"]="surface"
all_data.loc[(all_data["Model Position Just Chain"]=="surface") & (all_data["Model Position Complex"]=="intermediate"), "Model Position"]="surface"

all_data["Final position"]=all_data["Phos position"]
all_data["Final position"].fillna(all_data["Model position"],inplace=True)

#%% [quantification of locations]

for_count=all_data[["uniprot","position Just Chain","Model position Just Chain"\
                      ,"Model position Complex","position Complex","Phos position","Model position","Final position","PDB_ID"]]
for_count.rename(columns={"uniprot":"Protein Group"}, inplace=True)

for_count.groupby("position Just Chain")["Protein Group"].count()
for_count.groupby("position Complex")["Protein Group"].count()
for_count.groupby("Phos position")["Protein Group"].count()

for_count[["position Just Chain","Protein Group"]].groupby("position Just Chain").nunique()
for_count[["position Complex","Protein Group"]].groupby("position Complex").nunique()
for_count[["Phos position","Protein Group"]].groupby("position").nunique()

for_count.groupby("Model position Just Chain")["Protein Group"].count()
for_count.groupby("Model position Complex")["Protein Group"].count()
for_count.groupby("Model position")["Protein Group"].count()

for_count.groupby("Final position")["Protein Group"].count()

#%% [extraction of data]

all_data.to_excel("all_phospho_paper_result_data.xlsx")
