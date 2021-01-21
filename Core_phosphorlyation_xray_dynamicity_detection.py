#%% [packages]
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
import ast
from itertools import compress
from itertools import chain
import seaborn as sns
import matplotlib.pyplot as plt

#%% [data acquisiton]
os.chdir(r"location")

all_data=pd.read_excel("all_phospho_paper_result_data.xlsx")
all_data=all_data.loc[all_data["Final Position"].notnull()]
all_data=all_data.loc[all_data["Final Position"]=="core"].reset_index()

all_data.PDB_information=list(map(ast.literal_eval,all_data.PDB_information))
all_data["Model PDB"]=all_data["Model PDB"].fillna("[]")
all_data["Model PDB"]=list(map(ast.literal_eval,all_data["Model PDB"]))

#%% [pdb filtration]

#check all pdb files which match with phospho_site
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
            temp=list(compress(temp, [y for y in x]))

        phospho_check.append(temp)

# select the ones coverages close to each other by 90%
closes=[]
for i in phospho_check:
    if i==[]:
        closes.append(i)
    else:
        cov=np.array([int(y[1].split(",")[0].split("=")[1].split("-")[1])\
                      -int(y[1].split(",")[0].split("=")[1].split("-")[0]) for y in i])
        closes.append(list(compress(i, (cov)>=max(cov)*0.9)))

all_data["PDB_filtered"]=closes

filtered=all_data["PDB_filtered"][all_data["PDB_filtered"].apply(len) > 0]

list_of_pdb_filtered=list(set([a[0] for a in list(chain.from_iterable(filtered))]))

#%% [writed all pdbs]
os.chdir(r"location")
with open('all_pdb_file.txt', 'w') as f:
    for item in list_of_pdb_filtered:
        f.write("%s\n" % item)

#%%

#download all pdb files PDBdownload.py
os.chdir("location")
import PDBdownload


#%% [create files for monomer format]
#get all the pdb files
mypath=r"D:\Anaconda\Thesis\All_phospho\PDBs"
lesspdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
lesspdbfiles=[x for x in lesspdbfiles if x.endswith(".ent")]
#filter filtered list for empty lists
filtered_2=[x for x in phospho_check if x!=[]]
filtered_2=[item for sublist in filtered_2 for item in sublist]

os.chdir(r"location\PDBs")

#create files only with desired chain
for less in lesspdbfiles:
    chain=[t for t in filtered_2 if (t[0]==less.split(".")[0].split("pdb")[1].upper())]
    for alp in chain:
        with open(less,"r") as f:
            with open(less.split(".")[0]+"_"+alp[1].split("=")[0].split("/")[0]+"_less.pdb","w") as f1:
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

#%% [monomer rsa scores]

#get just chain freesasa result
os.chdir(r"location\PDBs")

mypath=r"location\PDBs"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]


data_for_PDB=all_data[["PDB_filtered","PDB_places","position","phospho window"]]
data_for_PDB.rename(columns={"position":"Real Places"},inplace=True)
data_for_PDB.rename(columns={"phospho window":"Sequence Window"},inplace=True)

#3 letter to 1 letter dictionary
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



c=[]
freesasa_results=[]
for o, row in data_for_PDB.iterrows():
    print(o)
    #if there is no PDB
    if row["PDB_filtered"]==[]:
        freesasa_results.append("")
        c.append([o,0])
    else:
        temp=[]
        for x in row["PDB_filtered"]:
            #save name of pdb
            pdb_name=x[0].lower()+"_"+x[1].split("=")[0].split("/")[0]+"_less"

            if not "pdb{}.pdb.out".format(pdb_name) in pdbfiles:
                temp.append("")
                c.append([o,2])
                continue
            else:
                #if there is no match with phospho site
                lines=[]
                for line in pd.read_csv("pdb{}.pdb.out".format(pdb_name),skiprows=[0,1,2,3,4,5,6], engine="python",chunksize=1, encoding="utf-8"):
                    lines.append(line.iloc[0,0])
                lines_2=[]
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
                pdb_dict=pd.read_csv(io.StringIO('\n'.join(lines_2)), delim_whitespace=True)
                count=0
                for a in pdb_dict.index:
                    if a==(('END', 'Absolute', 'sums')):
                        ind=count
                        break
                    else:
                        count=count+1
                        pass
                pdb_dict=pdb_dict[0:ind]

                resetted=pdb_dict.reset_index()
                seq=resetted["level_1"].tolist()
                seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
                found=0
                for i in range(0,7):
                    place=seqs.find(row["Sequence Window"][i:i+7])
                    if place==-1:
                        continue
                    else:
                        found=place+6-i
                        break
                if found==0:
                    temp.append("")
                else:
                    temp.append(resetted.iloc[found]["REL"])
        freesasa_results.append(temp)


def to_float(x):
    if x=="":
        return x
    else:
        return float(x)

freesasa_results=["" if x==["",""] else x for x in freesasa_results]
freesasa_results=[list(map(to_float,i)) if (i!="") else i for i in freesasa_results]

all_data["freesasa result"]=freesasa_results



#%% [whole protein RSA scores]

#get just chain freesasa result
os.chdir(r"location\PDBs")

mypath=r"location\PDBs"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]
pdbfiles=[x for x in pdbfiles if not "less" in x]
pdbfiles=[x for x in pdbfiles if not "bundle" in x]


data_for_PDB=all_data[["PDB_filtered","PDB_places","position","phospho window"]]
data_for_PDB.rename(columns={"position":"Real Places"},inplace=True)
data_for_PDB.rename(columns={"phospho window":"Sequence Window"},inplace=True)

#3 letter to 1 letter dictionary
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



c=[]
freesasa_results=[]
for o, row in data_for_PDB.iterrows():
    print(o)
    #if there is no PDB
    if row["PDB_filtered"]==[]:
        freesasa_results.append("")
        c.append([o,0])
    else:
        temp=[]
        for x in row["PDB_filtered"]:
            #save name of pdb
            pdb_name=x[0].lower()

            if not "pdb{}.ent.out".format(pdb_name) in pdbfiles:
                temp.append("")
                c.append([o,2])
                continue
            else:
                #if there is no match with phospho site
                lines=[]
                for line in pd.read_csv("pdb{}.ent.out".format(pdb_name),skiprows=[0,1,2,3,4,5,6], engine="python",chunksize=1, encoding="utf-8"):
                    lines.append(line.iloc[0,0])
                lines_2=[]
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
                pdb_dict=pd.read_csv(io.StringIO('\n'.join(lines_2)), delim_whitespace=True)
                count=0
                for a in pdb_dict.index:
                    if a==(('END', 'Absolute', 'sums')):
                        ind=count
                        break
                    else:
                        count=count+1
                        pass
                pdb_dict=pdb_dict[0:ind]

                resetted=pdb_dict.reset_index()
                seq=resetted["level_1"].tolist()
                seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
                found=0
                for i in range(0,7):
                    place=seqs.find(row["Sequence Window"][i:i+7])
                    if place==-1:
                        continue
                    else:
                        found=place+6-i
                        break
                if found==0:
                    temp.append("")
                else:
                    temp.append(resetted.iloc[found]["REL"])
        freesasa_results.append(temp)


def to_float(x):
    if x=="":
        return x
    else:
        return float(x)

freesasa_results=["" if x==["",""] else x for x in freesasa_results]
freesasa_results=[list(map(to_float,i)) if (i!="") else i for i in freesasa_results]

all_data["freesasa intact result"]=freesasa_results


#%% [categorize according to dynamicity]
cores=all_data.loc[all_data["Final Position"]=="core"].reset_index()

total_change=0

for i,row in cores.iterrows():
    if row["freesasa result"]=="":
        pass
    elif "" in row["freesasa result"]:
        c=[x for x in row["freesasa result"] if x!=""]
        if any(np.array(c)>5):
            cores.loc[i,"just change"]="changed"
    else:
        if any(np.array(row["freesasa result"])>5):
            cores.loc[i,"just change"]="changed"

for i,row in cores.iterrows():
    if row["freesasa intact result"]=="":
        pass
    elif "" in row["freesasa intact result"]:
        c=[x for x in row["freesasa intact result"] if x!=""]
        if any(np.array(c)>5):
            cores.loc[i,"intact change"]="changed"
    else:
        if any(np.array(row["freesasa intact result"])>5):
            cores.loc[i,"intact change"]="changed"

cores.to_excel("core_phosphosites_location_change.xlsx", index=False)


#%% [dynamicity analysis]

os.chdir(r"location")

# read the file
change_core=pd.read_excel("core_phosphosites_location_change.xlsx")
# seperate PDB names as list
change_core.PDB_filtered=list(map(ast.literal_eval,change_core.PDB_filtered))
# eliminate empty RSA scores
change_core=change_core.loc[~change_core["freesasa result"].isnull()]
change_core["freesasa result"]=change_core["freesasa result"].fillna("[]")
change_core["freesasa result"]=list(map(ast.literal_eval,change_core["freesasa result"]))
# fill-out remaining ones with unchanged tag
change_core["just change"]=change_core["just change"].fillna("unchanged")
# if there are less than 3 structure, we filtered them
change_core["freesasa result"] = change_core["freesasa result"].apply(lambda x: [num for num in x if isinstance(num, (int,float))])
change_core=change_core.loc[change_core["freesasa result"].apply(len)>2]
# finded maximum RSA change between structures
for o, row in change_core.iterrows():
    change_core.loc[o,"change amount"]=float(max(row["freesasa result"]))-float(min(row["freesasa result"]))
    change_core.loc[o,"freesasa_mean"]=np.mean(row["freesasa result"])
# change tags with proper naming
change_core['just change'] = change_core['just change'].map({'changed': 'Dynamic\nCore', 'unchanged': 'Static\nCore'})

#%% [swarm plot from maximum RSA score change between structures]
sns.set(style="whitegrid", palette="Accent",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(6, 9))
# Draw a categorical scatterplot to show each observation
sns.swarmplot(x="Final Position",y="change amount", hue="just change",
              palette=["r", "c"], data=change_core)
sns.despine(left=True)
ax.set_xlabel("Core Phosphorylations",size = 20,alpha=0.7)
ax.set_xticklabels("")
ax.set_ylabel("Maximum RSA Score Change",size = 20,alpha=0.7)
ax.legend(title='', prop={'size': 12})
plt.savefig("swarmplot_Xray_core_phosphorylations.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [optimal RSA score comparison]

all_paper_data=pd.read_excel("all_phospho_paper_result_data.xlsx")


merged_data=change_core.merge(all_paper_data, how="left", on=["uniprot","position"])

sns.set(style="whitegrid", palette="Paired_r",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(5, 5))
ax = sns.boxplot(x=merged_data['just change'], y=merged_data["freesasa result_y"], fliersize=4, notch=True, width=0.4)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7)
ax.set_ylabel("Optimal Structure RSA Score",size = 20,alpha=0.7)

plt.savefig("core_phosphorylations_relative_sasa_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [normalized b-factor comparison]


os.chdir(r"locations\changed_bfactor")

filtered=change_core["PDB_filtered"][change_core["PDB_filtered"].apply(len) > 0]

list_of_pdb_filtered=list(set([a[0] for a in list(chain.from_iterable(filtered))]))

with open("pdbs.txt", "w+") as f:
    for item in list_of_pdb_filtered:
        f.write("%s\n" % item)

pdbs=pd.read_csv("pdbs.txt",sep=" ",header=None)[0].tolist()

# download b-factor files
mypath=r"locations\changed_bfactor"
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for pdb in pdbs:
    print(pdb)
    if "{}_bfactor.txt".format(pdb) in pdbfiles:
        pass
    else:
        try:
            url = 'http://ignm.ccbb.pitt.edu/work/{}/{}.bfactor'.format(pdb,pdb)

            response = urllib.request.urlopen(url)
            webContent = response.read()

            f = open('{}_bfactor.txt'.format(pdb), 'wb')
            f.write(webContent)
            f.close
        except:
            continue

#%%
os.chdir(r"locations\changed_bfactor")

mypath=r"locations\changed_bfactor"
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "bfactor" in x]
pdbfiles=[x for x in pdbfiles if "txt" in x]


# if is there any number in a list
def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# collect all b-factor results
c=[]
bfactor=[]
for o, row in change_core.iterrows():
    print(o)
    #if there is no PDB
    if row["PDB_filtered"]==[]:
        bfactor.append(np.nan)
        c.append([o,0])
    else:
        temp=[]
        for x in row["PDB_filtered"]:
            #save name of pdb
            pdb_name=x[0].upper()
            if not "{}_bfactor.txt".format(pdb_name) in pdbfiles:
                temp.append(np.nan)
                c.append([o,2])
                continue
            else:
                out_file="{}_bfactor.txt".format(pdb_name)
                #if there is no match with phospho site
                with open(out_file, "r") as f:
                    lines=[]
                    for line in f:
                        if hasNumbers(line.split()[1]):
                            blob=line.split()
                            if len(blob[1])==7:
                                blob[1]=blob[1][0:3]+" "+blob[1][3::]
                                blob=" ".join(blob)
                                lines.append(blob)
                            else:
                                lines.append(line)
                        else:
                            lines.append(line)
                bfile=pd.read_csv(io.StringIO('\n'.join(lines)), delim_whitespace=True, names=["chain identifier", "residue type", "node index", "GNMbpredicted B- factor",  "X-ray crystallographic B factor"])
                chain=x[1].split("=")[0].split("/")[0]
                bfile=bfile.loc[bfile["chain identifier"].apply(str)==chain]
                seq=bfile["residue type"].tolist()
                seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
                found=0
                for i in range(0,7):
                    place=seqs.find(row["phospho window"][i:i+7])
                    if place==-1:
                        continue
                    else:
                        found=place+6-i
                        break
                if bfile.iloc[found]["X-ray crystallographic B factor"]==0:
                    print(o)
                    temp.append(np.nan)
                elif len(set(bfile["X-ray crystallographic B factor"].tolist()))==1:
                    print(o,2)
                    temp.append(np.nan)
                else:
                    temp.append(bfile.iloc[found]["X-ray crystallographic B factor"]/np.mean(bfile["X-ray crystallographic B factor"]))

        bfactor.append(temp)

change_core["x_ray_bfactor"]=bfactor


c=[]
bfactor=[]
for o, row in change_core.iterrows():
    print(o)
    #if there is no PDB
    if row["PDB_filtered"]==[]:
        bfactor.append(np.nan)
        c.append([o,0])
    else:
        temp=[]
        for x in row["PDB_filtered"]:
            #save name of pdb
            pdb_name=x[0].upper()
            if not "{}_bfactor.txt".format(pdb_name) in pdbfiles:
                temp.append(np.nan)
                c.append([o,2])
                continue
            else:
                out_file="{}_bfactor.txt".format(pdb_name)
                #if there is no match with phospho site
                with open(out_file, "r") as f:
                    lines=[]
                    for line in f:
                        if hasNumbers(line.split()[1]):
                            blob=line.split()
                            if len(blob[1])==7:
                                blob[1]=blob[1][0:3]+" "+blob[1][3::]
                                blob=" ".join(blob)
                                lines.append(blob)
                            else:
                                lines.append(line)
                        else:
                            lines.append(line)
                bfile=pd.read_csv(io.StringIO('\n'.join(lines)), delim_whitespace=True, names=["chain identifier", "residue type", "node index", "GNMbpredicted B- factor",  "X-ray crystallographic B factor"])
                chain=x[1].split("=")[0].split("/")[0]
                bfile=bfile.loc[bfile["chain identifier"].apply(str)==chain]
                seq=bfile["residue type"].tolist()
                seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
                found=0
                for i in range(0,7):
                    place=seqs.find(row["phospho window"][i:i+7])
                    if place==-1:
                        continue
                    else:
                        found=place+6-i
                        break
                if bfile.iloc[found]["GNMbpredicted B- factor"]==0:
                    print(o)
                    temp.append(np.nan)
                elif len(set(bfile["GNMbpredicted B- factor"].tolist()))==1:
                    print(o,2)
                    temp.append(np.nan)
                else:
                    temp.append(bfile.iloc[found]["GNMbpredicted B- factor"]/np.mean(bfile["GNMbpredicted B- factor"]))

        bfactor.append(temp)

change_core["gnm_bfactor"]=bfactor
# get mean b-factor
change_core["mean_gnm"]=change_core.gnm_bfactor.apply(np.nanmean)
change_core["mean_xray"]=change_core.x_ray_bfactor.apply(np.nanmean)

#statistical analysis
stats.kstest(change_core.functional_score.dropna(),"norm")
stats.kstest(change_core.functional_score.dropna().apply(np.log2),"norm")

stats.ttest_ind(changed["mean_gnm"],unchanged["mean_gnm"])
stats.ttest_ind(changed["mean_xray"],unchanged["mean_xray"], equal_var=False)

data=change_core.dropna(subset=["mean_xray"])
# boxplot for b-factor
sns.set(style="whitegrid", palette="Paired_r",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(5, 8))
ax = sns.boxplot(x=data['just change'], y=data.mean_xray, fliersize=4, notch=True, width=0.4)
ax = sns.swarmplot(x=data['just change'], y=data.mean_xray, color=".35", size=3.7)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7)
ax.set_ylabel("Mean(Normalized B-factor)",size = 20,alpha=0.7)

plt.savefig("core_phosphorylations_bfactor_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [disorder score analysis]

os.chdir(r"location")

phospho=pd.read_csv("All_cores_secondary.csv")
phospho["uniprot"]=phospho.id.apply(lambda x: x.split("_")[1])

known_phos=change_core.merge(phospho, how="left", left_on=["uniprot","position"], right_on=["uniprot","n"])

changed=known_phos.loc[known_phos["just change"]=="Dynamic\nCore"]
unchanged=known_phos.loc[known_phos["just change"]=="Static\nCore"]


stats.kstest(known_phos["disorder"].dropna(),"norm")
stats.kstest(known_phos["disorder"].dropna().apply(np.log2).replace([np.inf, -np.inf], np.nan).dropna(),"norm")

stats.ttest_ind(changed["logdisorder"].dropna(),unchanged["logdisorder"].dropna(), equal_var=False)

sns.set(style="whitegrid", palette="Paired_r",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(5, 8))
ax = sns.boxplot(x=known_phos['just change'], y=known_phos.logdisorder, fliersize=4, notch=True, width=0.4)
ax = sns.swarmplot(x=known_phos['just change'], y=known_phos.logdisorder, color=".35", size=3.7)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7)
ax.set_ylabel("log2(Disorder Score)",size = 20,alpha=0.7)

plt.savefig("core_phosphorylations_disorder_scores_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [secondary structure analysis]

# statistical analysis
for i in range(6):
    dy=known_phos.groupby(["just change","q8"]).size()["Dynamic\nCore"]
    st=known_phos.groupby(["just change","q8"]).size()["Static\nCore"]
    a,b=stats.fisher_exact(np.array([[dy[i],len(changed)-dy[i]],[st[i],len(unchanged)-st[i]]]))
    if b<0.05:
        print(b)
        print(known_phos.groupby(["just change","q8"]).size().index[i])


category_names = ['Coil', '\u03B2-Sheet',
                 r'$ 3_{10}-helix $', r'$ \alpha-helix $', r'$ \pi-helix $', 'Bend',"Turn"]

# secondary structure information is obtained in two lists
known_phos.groupby(["just change","q8"]).size()

results = {
    'Dynamic Core (%)': [57/149*100, 37/149*100, 5/149*100, 40/149*100, 0/149*100, 6/149*100, 4/149*100],
    'Static Core (%)': [38/173*100, 73/173*100, 5/173*100, 45/173*100, 2/173*100, 2/173*100, 8/173*100],
}


def survey(results, category_names):

    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(18, 2))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.6,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.1 else 'dimgray'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            if int(c)==0:
                pass
            else:
                ax.text(x, y, str(int(c)), ha='center', va='center',
                        color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small', handlelength=0.8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    return fig, ax


survey(results, category_names)
plt.savefig("core_phosphorylations_secondary_structure.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [knwon phosphorylation sites detection]

os.chdir(r"location")

phospho=pd.read_csv("Phosphorylation_site_dataset",sep="\t" ,skiprows=[0,1])
phospho=phospho.loc[phospho.ORGANISM=="human"]
phospho["position"]=phospho.MOD_RSD.apply(lambda x: int(x.split("-")[0][1::]))

known_phoshorylations=change_core.merge(phospho, how="left", left_on=["uniprot","position"], right_on=["ACC_ID","position"])

known_phoshorylations["HL"]=np.sum(known_phoshorylations[["MS_LIT","MS_CST"]],axis=1)

known_phoshorylations.groupby("just change")[["HL","LT_LIT"]].agg(sum)

# bar plot for known phosphorylation site distribution


category_names = ['Entry in Phosphositeplus', 'Absent in Phosphositeplus']
results = {
    'Dynamic Core': [115/149*100, 34/149*100],
    'Static Core': [110/173*100, 63/173*100]
}


def survey(results, category_names, size):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.2, 0.75, data.shape[1]))

    fig, ax = plt.subplots(figsize=size)
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax.text(x, y, str(int(c))+" %", ha='center', va='center',
                    color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small', handlelength=0.9)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return fig, ax


survey(results, category_names,(9, 1))
plt.savefig("known_phospho_bar_plot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [funtional score analysis]

os.chdir(r"location")
result=pd.read_excel("all_phospho_paper_result_data.xlsx")
result=result.loc[result["Final Position"].notnull()]

changed=change_core.loc[change_core["just change"]=="Dynamic\nCore"]
unchanged=change_core.loc[change_core["just change"]=="Static\nCore"]

all_comp=pd.concat([changed[["just change","functional_score"]],unchanged[["just change","functional_score"]],\
                   result.loc[result["Final Position"]!="core"][["Final Position","functional_score"]]], ignore_index=True)

all_comp["Final Position"]=all_comp["Final Position"].fillna(all_comp["just change"])

# box plot for functional score comparison
sns.set(style="whitegrid", palette="Set3",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(8, 8))
ax = sns.boxplot(x=all_comp['Final Position'], y=all_comp.functional_score, fliersize=4, notch=True, width=0.6)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7)
ax.set_ylabel("Functional Score",size = 20,alpha=0.7)
ax.set(xticklabels=["Dynamic\nCore","Static\nCore","Intermediate","Surface","Interface"])
plt.savefig("core_phosphorylations_functional_scores_comparison_withall_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)


# statistical analysis


df = [np.array(x["functional_score"]) for _, x in all_comp.groupby('Final Position')]
stats.f_oneway(df[0],df[1],df[2],df[3],df[4])

mc1=multi.MultiComparison(all_comp["functional_score"],all_comp["Final Position"]\
                          .replace({"Dynamic\nCore":"Dynamic","Static\nCore":"Static"}))
res1 = mc1.tukeyhsd()
print(res1.summary())

#%% [correlation between b-factor and rsa]

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

# to find each structure RSA and corresponding b-factor
def gnm_free_corr(case):
    case_df=pd.DataFrame(columns=["freesasa","gnm_bfactor","xray_bfactor"])

    for x in case["PDB_filtered"]:
        temp=[]
        if (not "{}_bfactor.txt".format(x[0].upper()) in pdbfiles) or (not "pdb{}_less.pdb.out".format(x[0].lower()+"_"+x[1].split("=")[0].split("/")[0]) in pdbfiles):
            pass
        else:
             lines=[]
             for line in pd.read_csv("pdb{}.pdb.out".format(x[0].lower()+"_"+x[1].split("=")[0].split("/")[0]+"_less"),skiprows=[0,1,2,3,4,5,6], engine="python",chunksize=1, encoding="utf-8"):
                 lines.append(line.iloc[0,0])
             lines_2=[]
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
             pdb_dict=pd.read_csv(io.StringIO('\n'.join(lines_2)), delim_whitespace=True)
             count=0
             for a in pdb_dict.index:
                 if a==(('END', 'Absolute', 'sums')):
                     ind=count
                     break
                 else:
                     count=count+1
                     pass
             pdb_dict=pdb_dict[0:ind]

             resetted=pdb_dict.reset_index()
             seq=resetted["level_1"].tolist()
             seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
             found=0
             for i in range(0,7):
                place=seqs.find(case["phospho window"][i:i+7])
                if place==-1:
                    continue
                else:
                    found=place+6-i
                    break
             if found==0:
                 free=np.nan
                 print(x[0])
             else:
                 free=float(resetted.iloc[found]["REL"])

             out_file="{}_bfactor.txt".format(x[0].upper())
             #if there is no match with phospho site
             with open(out_file, "r") as f:
                 lines=[]
                 for line in f:
                    if hasNumbers(line.split()[1]):
                         blob=line.split()
                         if len(blob[1])==7:
                             blob[1]=blob[1][0:3]+" "+blob[1][3::]
                             blob=" ".join(blob)
                             lines.append(blob)
                         else:
                             lines.append(line)
                    else:
                         lines.append(line)
             bfile=pd.read_csv(io.StringIO('\n'.join(lines)), delim_whitespace=True, names=["chain identifier", "residue type", "node index", "GNMbpredicted B- factor",  "X-ray crystallographic B factor"])
             chain=x[1].split("=")[0].split("/")[0]
             bfile=bfile.loc[bfile["chain identifier"].apply(str)==chain]
             seq=bfile["residue type"].tolist()
             seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq])
             found=0
             for i in range(0,7):
                 place=seqs.find(case["phospho window"][i:i+7])
                 if place==-1:
                     continue
                 else:
                     found=place+6-i
                     break
             if bfile.iloc[found]["GNMbpredicted B- factor"]==0:
                 gnm=np.nan
             elif len(set(bfile["GNMbpredicted B- factor"].tolist()))==1:
                 gnm=np.nan
             else:
                 gnm=bfile.iloc[found]["GNMbpredicted B- factor"]/np.mean(bfile["GNMbpredicted B- factor"])



             if bfile.iloc[found]["X-ray crystallographic B factor"]==0:
                 xray=np.nan
             elif len(set(bfile["X-ray crystallographic B factor"].tolist()))==1:
                 xray=np.nan
             else:
                 xray=bfile.iloc[found]["X-ray crystallographic B factor"]/np.mean(bfile["X-ray crystallographic B factor"])

             case_df = case_df.append({"freesasa":free,"gnm_bfactor":gnm,"xray_bfactor":xray}, ignore_index=True)

    case_df.dropna(inplace=True)

    return case_df
# save them in a dataframe
case_study=pd.DataFrame(columns=["freesasa","gnm_bfactor","xray_bfactor"])
for i,row in change_core.iterrows():
    case_df_2=gnm_free_corr(row)
    if case_df_2.empty:
        pass
    else:
        case_study=pd.concat([case_study,case_df_2[["freesasa","gnm_bfactor","xray_bfactor"]]],ignore_index=True)

case_study.dropna(inplace=True)
# look at the correlation
sns.set(font_scale = 1.6, style="white")
g = sns.jointplot(x=case_study.freesasa, y=case_study.xray_bfactor,
                  kind="reg", truncate=False,
                  color="#32a887", height=10)
g.set_axis_labels("RSA score", "Normalized B-factor", fontsize=20)
g.annotate(stats.pearsonr)

plt.savefig("bfactor_freesasa_corr.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)
