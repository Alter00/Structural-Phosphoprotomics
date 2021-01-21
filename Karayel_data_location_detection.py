#%% [packages]

import pandas as pd
import numpy as np
from Bio.PDB import *
from bs4 import BeautifulSoup
import urllib
import requests
import xmltodict
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

#%% [structural determination of Karayel data]
os.chdir(r"location")

#read file as dataframe for 3 cases
all_data=pd.read_excel("Supplemantary 4 for processing.xlsx",header=[0, 1],sheetname="Supplementary Table 4",index_col=None)

#add state column
all_data.loc[(all_data["Interphase/Cytokinesis"]==True) &\
             (all_data["Mitosis/Cytokinesis"]==False),"State"]="dynamic cytokinesis"
all_data.loc[all_data["Mitosis/Interphase"]==True &\
             (all_data["Mitosis/Cytokinesis"]==False),"State"]="dynamic mitosis"
all_data.loc[(all_data["Mitosis/Cytokinesis"]==True),"State"]="dynamic both"

all_data.loc[(all_data["Mitosis/Cytokinesis"]==False) & (all_data["Interphase/Cytokinesis"]==False)\
 & (all_data["Mitosis/Interphase"]==False),"State"]="stable"


 #%%
 #read fasta file
 fasta_file="sprot_2014_08.fasta"

 with open(fasta_file,"r") as fastafile:

     seqs={}
     s=""
     name=None

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


 fasta_divided=pd.DataFrame(fasta["ids"].str.split("|").tolist(),columns=("indicator","Protein Group","description"))

 fasta=fasta.join(fasta_divided,how="left")
 fasta=fasta.drop(["ids","indicator"],axis=1)
 #create a small df to merge with sequence and add it to the main data
 temp=all_data[["Protein Group","Modification"]]

 temp=temp.merge(fasta,on="Protein Group",how="left")

 all_data["Protein Sequence"]=list(temp["sequence"])


#create phospho windows with phosphorylation higher than 50% probability
phospho_windows=[]
for z in range(len(all_data.index)):
    #store each row's phospho site
    places_scores=np.array(all_data.iloc[z]["Phospho_Modification with PhosphoRS probability"].replace("\"","").replace("; ",":").replace("(",":").replace(")",":").replace(": ","").split(":"))
    if places_scores[1::3][places_scores[2::3].astype(float)>50].size:
        temp=[]
        for place in places_scores[1::3][places_scores[2::3].astype(float)>50].astype(int):
            length=len(all_data["Protein Sequence"][z])
            fin=all_data["Protein Sequence"][z].find(all_data["Sequence"][z])
            if (length-place-fin)<7:
                temp.append(all_data["Protein Sequence"][z][(fin+place-7):(length)]+(7-length+fin+place)*"_")
            elif (fin+place<8):
                temp.append((((7-fin-place))*"_")+all_data["Protein Sequence"][z][0:(fin+place+6)])
            else:
                temp.append(all_data["Protein Sequence"][z][(fin+place-1-6):(fin+place+6)])
        phospho_windows.append(temp)
    else:
        phospho_windows.append([])

#make single list string
phospho_windows=["".join(x) if len(x)==1 else x for x in phospho_windows]
#put new phospho window
all_data["Seq Window"]=phospho_windows
#filter proteins without any seq window
all_data=all_data[all_data["Seq Window"].map(lambda d: len(d))>0]
all_data.reset_index(inplace=True,drop=True)

all_data=all_data.explode('Seq Window').reset_index(drop=True)
#create small dataframe to store phospho places
data_for_place=all_data.copy()
data_for_place.rename(index=str,columns={"Protein Sequence":"Total Sequence"},inplace=True)

"""
#filter phosphodata which have lower than 70% probability
phospho_windows["Phospho_Modification with PhosphoRS probability"].apply(lambda x: x.replace("; ",":").replace("(",":").replace(")",":").replace(": ","").split(":")[2::3])
"""

#get real places
real_places=[]
for i in range(len(data_for_place)):
    if data_for_place.iloc[i]["Seq Window"].startswith("_"):
            number_3=data_for_place.iloc[i]["Seq Window"].count("_")
            real_places.append(data_for_place.iloc[i]["Total Sequence"].find(data_for_place.iloc[i]["Seq Window"].replace("_","")[0])+1+6-number_3)
    elif data_for_place.iloc[i]["Seq Window"].endswith("_"):
        number_4=data_for_place.iloc[i]["Seq Window"].count("_")
        real_places.append(data_for_place.iloc[i]["Total Sequence"].find(data_for_place.iloc[i]["Seq Window"].replace("_",""))+1+6)
    else:
        real_places.append(data_for_place.iloc[i]["Total Sequence"].find(data_for_place.iloc[i]["Seq Window"])+1+6)

#arrange real places in proper format
formatted_real_place=[]
for i in range(len(real_places)):
    if type(real_places[i])==list:
        temp=[]
        for each in real_places[i]:
            temp.append(data_for_place.iloc[i]["Total Sequence"][each-1]+":"+str(each))
        formatted_real_place.append(temp)
    else:
        formatted_real_place.append(data_for_place.iloc[i]["Total Sequence"][real_places[i]-1]+":"+str(real_places[i]))
all_data["Real Places"]=formatted_real_place
#create an excel file to obtain PDB ID's from uniprot
all_data["Position"]=all_data["Real Places"].apply(lambda x: int(x.split(":")[1]))
all_data["Residue"]=all_data["Real Places"].apply(lambda x: x.split(":")[0])

#%%

# get all the information from uniprot page one by one if there is no site, put empty and if no answer wait 5 second for it
pdb_and_coverage=[]
count=1
for uniprot in data_for_place["Protein Group"]:
    url = "http://www.uniprot.org/uniprot/"+uniprot+".xml"
    print(url)
    session = requests.Session()
    retry = Retry(connect=4, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)

    response = ""

    while response == "":

        try:
            response = session.get(url, timeout=5)
            thedata = response.content
            soup=BeautifulSoup(thedata, "xml")
            references=soup.find_all("dbReference")
            references_2=[str(x) for x in references if "type=\"PDB\"" in str(x) and "value=" in str(x)]
            temp=[]
            for ref in references_2:
                if ("id=\"" in str(ref)) & ("\"chains\" value=\"" in str(ref)) & ("resolution\" value=" in str(ref)):
                    temp.append([str(ref)[str(ref).find("id=\"")+4:str(ref).find("\" type")],str(ref)[str(ref).find("\"chains\" value=\"")+16:str(ref).find("\"/>\n</dbReference>")],float(str(ref)[str(ref).find("resolution\" value=\"")+19:str(ref).find("\"/>\n<property type=\"chains\"")])])
                else:
                    pass
            pdb_and_coverage.append(temp)
        except urllib.error.HTTPError:
            pdb_and_coverage.append([])
            response="*"
        except:
            print("ZZzzzz...")
            print(url)
            time.sleep(5)
            continue
    print(str(count) + " is done")
    count=count+1

all_data["PDB_information"]=pdb_and_coverage
#create new dataframe for storing all the information (it is essential to save as string not list)
all_data_2=all_data.copy()
pdb_and_coverage_2=[]
for i in pdb_and_coverage:
    temp=[]
    for x in i:
        temp.append(":".join(list(map(str,x))))
    pdb_and_coverage_2.append(";".join(temp))

all_data_2["PDB_information"]=pdb_and_coverage_2

all_data_2["PDB_information"].to_csv("karayel_data_with_pdb_id.csv")
all_data_2["PDB_information"].to_excel("karayel_data_with_pdb_id.xlsx")


#%%
PDB_information_2=pd.read_excel("karayel_data_with_pdb_id.xlsx",header=[0])
PDB_information_2=list(PDB_information_2["PDB_information"])

new=[x.split(";") if type(x)==str else x for x in PDB_information_2 ]
new = [[s.split(":") if type(s)==str else s for s in l] if type(l)==list else l for l in new]

all_data["PDB_information"]=new


#%%
all_data["Position"]=all_data.Position.apply(int)

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
        phospho_site=row["Position"]
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

#if it is empty replace it with nan
all_data["PDB_filtered"]=phospho_check

list_of_pdb_filtered=list(set([a[0] for a in all_data["PDB_filtered"] if len(a)>0]))



with open('karayel.txt', 'w') as f:
    for item in list_of_pdb_filtered:
        f.write("%s\n" % item)

#%%

#download all pdb files PDBdownload.py
os.chdir("location")
import PDBdownload

#%%
os.chdir(r"D:\Anaconda\Karayel Thesis\Thesis\Modbase")

#model list are taken
empty_proteins=list(set(all_data.loc[all_data["PDB_filtered"].str.len()==0]["Protein Group"]))
#downloaded from modbase
for protein in empty_proteins:
    os.system("curl -L \"https://salilab.org/modbase/retrieve?databaseID={}\" -o {}.xml".format(protein,protein))

#%%

os.chdir(r"D:\Anaconda\Karayel Thesis\Thesis\Modbase")


mypath=r"D:\Anaconda\Karayel Thesis\Thesis\Modbase"
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
    model_dict[xml.split(".")[0]]=content
    bs_content = BeautifulSoup(content, "lxml")
    split_model=bs_content.find_all("pdbfile")
    split_model_2=[str(x) for x in split_model if "REMARK 220" in str(x)]
    quality=[]
    names=[]
    location=[]
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

pdb_model_names=pd.DataFrame(list(pdb_model_names.items()),columns=["Protein Group","Model PDB"])

all_data=all_data.merge(pdb_model_names, on="Protein Group", how="left")

#%%

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
        phospho_site=row["Position"]
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

with open('karayel_model.txt', 'w') as f:
    for item in list_of_pdb_filtered:
        f.write("%s\n" % item)
#%%

#download all pdb files PDBdownload.py
os.chdir("location")
import PDBdownload

%%
#change directory

os.chdir("location")


#get all the pdb files
mypath="location"
lesspdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
lesspdbfiles=[x for x in lesspdbfiles if x.endswith(".ent")]
#filter filtered list for empty lists
filtered_2=[x for x in phospho_check if x!=[]]

os.chdir("location")

#create files only with desired chain
for less in lesspdbfiles:
    chain=[t for t in filtered_2 if (t[0]==less.split(".")[0][3::].split("-")[0].upper())]
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

#%%


#change directory

os.chdir("location")


#get all the pdb files
mypath="location"
lesspdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
lesspdbfiles=[x for x in lesspdbfiles if x.endswith(".pdb")]
#filter filtered list for empty lists
filtered_2=[x for x in phospho_check_2 if x!=[]]

os.chdir("location")

#create files only with desired chain
for less in lesspdbfiles:
    chain=[t for t in filtered_2 if (t.split(":")[0]==less.split(".")[0].split("-")[0].upper())]
    for alp in chain:
        with open(less,"r") as f:
            with open(less.split(".")[0]+"_"+alp.split(":")[1].split("=")[0].split("/")[0]+"_less.pdb","w") as f1:
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

all_data.reset_index(drop=True,inplace=True)

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
data_for_PDB=all_data[["PDB_filtered","Position","Residue","Seq Window"]]
data_for_PDB.rename(columns={"Position":"Real Places"},inplace=True)
data_for_PDB.rename(columns={"Seq Window":"Sequence Window"},inplace=True)

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
all_data["Position Just Chain"]=position_categorization("freesasa result")
all_data["Position Complex"]=position_categorization("freesasa intact result")
all_data["Model Position Just Chain"]=position_categorization("model freesasa result")
all_data["Model Position Complex"]=position_categorization("model freesasa intact result")

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

all_data["Final Position"]=all_data["Phos Position"]
all_data["Final Position"].fillna(all_data["Model Position"],inplace=True)

#%% [quantification of locations]

for_count=all_data[["Uniprot","Position Just Chain","Model Position Just Chain"\
                      ,"Model Position Complex","Position Complex","Phos Position","Model Position","Final Position","PDB_ID"]]
for_count.rename(columns={"Uniprot":"Protein Group"}, inplace=True)

for_count.groupby("Position Just Chain")["Protein Group"].count()
for_count.groupby("Position Complex")["Protein Group"].count()
for_count.groupby("Phos Position")["Protein Group"].count()

for_count[["Position Just Chain","Protein Group"]].groupby("Position Just Chain").nunique()
for_count[["Position Complex","Protein Group"]].groupby("Position Complex").nunique()
for_count[["Phos Position","Protein Group"]].groupby("Position").nunique()

for_count.groupby("Model Position Just Chain")["Protein Group"].count()
for_count.groupby("Model Position Complex")["Protein Group"].count()
for_count.groupby("Model Position")["Protein Group"].count()

for_count.groupby("Final Position")["Protein Group"].count()

#%% [extraction of data]

all_data.to_excel("karayel_phospho_result_data.xlsx")


#%% [karayel data analysis]

os.chdir(r"location")

karayel=pd.read_excel("karayel_phospho_result_data.xlsx")
karayel=karayel.loc[karayel["Final Position"].notnull()]
karayel=karayel.drop_duplicates(subset=["Protein Group","Position"])



karayel.loc[(karayel[["Interphase/Cytokinesis","Mitosis/Cytokinesis","Mitosis/Interphase"]]==1)\
            .any(axis="columns"),"State"]="Dependent"
karayel.loc[~(karayel[["Interphase/Cytokinesis","Mitosis/Cytokinesis","Mitosis/Interphase"]]==1)\
            .any(axis="columns"),"State"]="Independent"
    

karayel["3log10cytokinesis"]=np.log10(karayel["relative fold changeCytokinesis"])
karayel["2log10mitosis"]=np.log10(karayel["relative fold changeMitosis "])
karayel["1log10interphase"]=np.log10(karayel["relative fold changeInterphase"])

karayel["identifier"]=karayel["Protein Group"]+"_"+karayel["Real Places"]
# filter empty states
karayel_3=karayel.dropna(subset=["State"])
# sorted based on position
karayel_3["Final Position"]=karayel_3["Final Position"].replace({"Interface":"interface"})
box_data=karayel_3.sort_values(by=["Final Position"], ascending=False)
df_plot = box_data.groupby(['State', 'Final Position']).size().reset_index().pivot(columns='Final Position', index='State', values=0)

# percentage of phosphorylation location
df_plot=df_plot.iloc[:, ::-1]
df_perc=(100. * df_plot / df_plot.sum()).round(1)
# statistical test (fisher_exact)
stats.chi2_contingency([df_plot.loc["Independent"].fillna(0),df_plot.loc["Dependent"].fillna(0)])

stats.fisher_exact([[5,2],[1,2]])

x = ["Surface","Intermediate","Interface","Core"]
y1=np.array(df_perc.loc["Dependent"])
y2=np.array(df_perc.loc["Independent"])

#phosphorylation site bar plot based on location and state


fig, ax = plt.subplots(figsize=(5,8))
# Create a bar chart in position bar_1
ax.tick_params(axis=u'both', which=u'both',length=0)
plt.rcParams["axes.labelsize"] = "medium"
plt.setp(ax.spines.values(), linewidth=0.4)
plt.rcParams["font.family"] = "Arial"

barwidth=0.60
numbers_1=["50%\n(56)","53%\n(73)","100%\n(5)","58%\n(14)"]
numbers_2=["50%\n(57)","47%\n(65)","0%\n(0)","42%\n(10)"]


colors = plt.get_cmap('Paired')(np.linspace(0.15, 0.9, 3))
# Create a bar chart in position bar_1
bar1=ax.bar(x, y1, label='Dependent', color="#003f5c",width=barwidth, alpha=0.8)
bar2=ax.bar(x, y2 ,bottom=y1,label='Independent', color="#ffa600",width=barwidth, alpha=0.8)

for idx,rect in enumerate(bar1):
    height = rect.get_height()/2+1.5
    ax.text(rect.get_x() + rect.get_width()/2., height,
            numbers_1[idx],
            ha='center', va='top', rotation=0,fontsize=18,
            fontstyle="oblique",fontweight="medium",family="Arial", color="white")

for idx,rect in enumerate(bar2):
    height = rect.get_height()/2+bar1[idx].get_height()+1.5
    ax.text(rect.get_x() + rect.get_width()/2., height,
            numbers_2[idx],
            ha='center', va='top', rotation=0,fontsize=18,
            fontstyle="oblique",fontweight="medium",family="Arial", color="black")




plt.yticks(fontsize=22)
plt.ylim(0,100)
plt.ylabel("Phospho-sites (%)",fontsize=26,weight="roman",family="Arial")
#ax.tick_params(axis = 'x', which = 'major', labelsize =24, rotation=45,ha='center')
plt.xticks(fontsize =24, rotation=60,rotation_mode="anchor",ha="right")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
legend=ax.legend(title='     Cell Cycle', prop={'size': 24},bbox_to_anchor=(0,1.3), mode="expand",frameon=False)
legend.get_title().set_fontsize('26')
legend._legend_box.align = "left"
plt.savefig("karayel_distribution_bar_plot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [karayel dynamicity analysis]


os.chdir(r"location")
# read core dynamicity data and filtered it, also seperate RSA score as list
change_core=pd.read_excel("core_phosphosites_location_change.xlsx")
change_core.PDB_filtered=list(map(ast.literal_eval,change_core.PDB_filtered))
change_core=change_core.loc[~change_core["freesasa result"].isnull()]
change_core["freesasa result"]=change_core["freesasa result"].fillna("[]")
change_core["freesasa result"]=list(map(ast.literal_eval,change_core["freesasa result"]))
change_core=change_core.loc[change_core["phospho window"].str[6].isin(["S","T","Y"])]
change_core["just change"]=change_core["just change"].fillna("unchanged")
change_core.rename(columns={"just change": "dynamicity"}, inplace=True)

change_core["freesasa result"] = change_core["freesasa result"].apply(lambda x: [num for num in x if isinstance(num, (int,float))])
change_core=change_core.loc[change_core["freesasa result"].apply(len)>2]

karayel_4=karayel_3.loc[karayel_3['Final Position']=="core"]


karayel_4=karayel_4.merge(change_core[["uniprot","position","dynamicity","freesasa result"]], \
                          how="left", left_on=["Protein Group","Position"], right_on=["uniprot","position"])

karayel_4.dropna(subset=["dynamicity"], inplace=True)
karayel_4.groupby(["dynamicity","State"]).size()
#hyper geom test to see statistical difference in distribution
hypergeom.pmf(M=14, n=10, k=6, N=4, loc=0)

# genename acquisiton
a=karayel_4["Protein Group"]
#uniprot id to genename
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'GENENAME',
'format': 'tab',
'query': "\t".join(a)
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
rawData = pd.read_csv(io.StringIO(response.decode('utf-8')),sep="\t")

karayel_4=karayel_4.merge(rawData, left_on="Protein Group", right_on="From", how="left")
# created identifiers for karayel phosphorylation sites
karayel_4["identifier"]=karayel_4.To+"_"+karayel_4.Residue+"-"+karayel_4.Position.apply(str)

karayel_4['dynamicity'] = karayel_4['dynamicity'].map({'changed': 'Dynamic Core', 'unchanged': 'Static Core'})

karayel_4['All_states'] = karayel_4['State']+" "+karayel_4['dynamicity']
# plot the figure
sns.set(style="whitegrid", palette="tab10",font_scale = 2,font="Arial")
fig, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(14, 8))


ax1=sns.stripplot(x="freesasa result_y", y="identifier", hue="dynamicity",\
                 data=Independent.explode("freesasa result_y"),jitter=True, alpha=.5,size=7, zorder=1,ax=ax1)

ax1.set_xlabel("",size = 20,alpha=0.7)
ax1.set_ylabel("Cell Cycle \nIndependent",size = 28,labelpad=12)
ax2=sns.stripplot(x="freesasa result_y", y="identifier", hue="dynamicity",\
                 data=Dependent.explode("freesasa result_y"),jitter=True, alpha=.5,size=7, zorder=1,ax=ax2)
ax2.set_xlabel("RSA score",size = 26)
ax2.set_ylabel("Cell Cycle \nDependent",size = 28,labelpad=27)
ax2.get_legend().remove()


leg=ax1.legend(loc="upper right", prop={'size': 26},markerscale=2)
leg.set_title("", prop = {'size':'large'})

plt.subplots_adjust(wspace=0, hspace=0)
ax2.autoscale()
plt.savefig("Karayel_core_change__plot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% abundance change plot
a=karayel.groupby("Final Position",as_index=False).mean()[["Final Position","3log10cytokinesis","2log10mitosis","1log10interphase"]]
b=karayel_4.groupby("dynamicity",as_index=False).mean()[["dynamicity","3log10cytokinesis","2log10mitosis","1log10interphase"]]
b.rename(columns={"dynamicity":"Final Position"},inplace=True)
result=pd.concat([a, b], ignore_index=True)

result=pd.melt(result,value_vars=["3log10cytokinesis","2log10mitosis","1log10interphase"], id_vars=["Final Position"])


result.sort_values(by="variable",inplace=True)
di={"1log10interphase":"Interphase","2log10mitosis":"Mitosis","3log10cytokinesis":"Cytokinesis"}
result.replace({"variable": di}, inplace=True)
die={"changed":"Dynamic Core","unchanged":"Static Core","surface":"Surface","intermediate":"Intermediate"}
result.replace({"Final Position": die}, inplace=True)
result=result.loc[result["Final Position"]!="core"]

fig, ax = plt.subplots(figsize=(18, 9))
plt.rcParams['font.family'] = 'Arial'

sns.set_context("paper", font_scale=2.3)
g=sns.lineplot(data=result, x="variable", y="value", hue="Final Position",ax=ax, linewidth=5, palette="tab10",legend=False)
plt.legend(frameon=False, loc='lower center', ncol=5, markerscale=1000)

ax.axhline(0, ls='--', color="black")
g.set(xlim=(-0.2,2.2))
g.set(ylim=(-0.25,1))
g.set(xlabel=None)
g.tick_params(bottom=False)
g.set(xticklabels=[])
g.set(ylabel="Median Fold Change")

ax.text(-0.1,-0.1,"Interphase",fontsize="medium")
ax.text(0.9,-0.1,"Mitosis",fontsize="medium")
ax.text(1.85,-0.1,"Cytokinesis",fontsize="medium")

plt.savefig("all_cluster.pdf", dpi=1500)
plt.savefig("all_cluster.svg", dpi=1500)


#%% file containing all abundance info
c=karayel.dropna(subset=["Final Position","3log10cytokinesis","2log10mitosis","1log10interphase"])[["identifier","Final Position","3log10cytokinesis","2log10mitosis","1log10interphase"]]
d=karayel_4.dropna(subset=["dynamicity","3log10cytokinesis","2log10mitosis","1log10interphase"])[["identifier","dynamicity","3log10cytokinesis","2log10mitosis","1log10interphase"]]

d.rename(columns={"dynamicity":"Final Position"},inplace=True)
result_2=pd.concat([c, d], ignore_index=True)
result_2["Final Position"]=result_2["Final Position"].replace({"changed":"dynamic","unchanged":"static"})

result_2.to_excel("Karayel_change.xlsx", index=False)