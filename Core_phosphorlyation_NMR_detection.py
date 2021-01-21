#%% [load packages]

import numpy as np
import pandas as pd
from scipy import stats
import os
import matplotlib.pyplot as plt
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
import re
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import statsmodels.stats.multicomp as multi
import seaborn as sns

#%% [uploaded location dataset]
os.chdir(r"location")

data_for_place=pd.read_excel("all_phospho_paper_result_data.xlsx")
data_for_place=data_for_place.loc[data_for_place["Final Position"].notnull()]
data_for_place=data_for_place.loc[data_for_place["Final Position"]=="core"].reset_index()

#%% [got all pdb information from uniprot]
pdb_and_coverage=[]
count=1
for uniprot in data_for_place["uniprot"]:
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
                if ("id=\"" in str(ref)) :
                    temp.append([str(ref)[str(ref).find("id=\"")+4:str(ref).find("\" type")],(str(ref)[str(ref).find("method\" value=\"")+15:str(ref).find("\"/>\n<property type=\"resolution\"")])])
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
# selected NMR structures
check=[[y  for y in x if ("NMR" in y[1])] for x in pdb_and_coverage ]

# got all NMR files list
list_of_pdb_filtered=list(set([a[0][0] for a in check if len(a)>0]))

nmrs=[a[0][0] if len(a)>0 else [] for a in check]

data_for_place["nmr"]=nmrs
# saved for further evaluation
data_for_place.to_csv("nmr.csv")
"""
#reload again
data_for_place=pd.read_csv("nmr.csv")
data_for_place=data_for_place.loc[data_for_place["phospho window"].str[6].isin(["S","T","Y"])]

nmrs=list(set(data_for_place.nmr.tolist()))
"""

#download all nmr files PDBdownload.py
os.chdir("location")
import PDBdownload


#%% [seperate NMRs based on models]

#read all files and store them in a dataframe
for less in pdbfiles:
    with open(less, 'r') as f:
        data = f.readlines()[1::]
        data="".join(data)

    found = re.findall(r'\n*(MODEL.*?\nENDMDL)\n*', data, re.M | re.S)

    [open("location\\models\\"+less+str(i)+'.pdb', 'w').write(found[i-1]) for i in range(1, len(found)+1)]

#%% [read all models]

os.chdir(r"location\models")

mypath=r"location\models"
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "out" in x]



#read all nmr files and store them in a dictionary
pdb_dict={}
for out_file in pdbfiles:
    lines=[]
    for line in pd.read_csv(out_file,skiprows=[0,1,2,3,4,5,6], engine="python",chunksize=1, encoding="utf-8"):
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
    look=pd.read_csv(io.StringIO('\n'.join(lines_2)), delim_whitespace=True)
    pdb_dict.update({out_file.split("pdb")[1].split(".")[0]+"_"+out_file.split("pdb")[1].split(".")[1].split("ent")[1] : look})
    print(out_file)


d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

data_for_PDB=data_for_place[["nmr","position","residue","phospho window"]]
data_for_PDB.rename(columns={"position":"Real Places"},inplace=True)
data_for_PDB.rename(columns={"phospho window":"Sequence Window"},inplace=True)

# read all nmr values for phospho_sites
c=[]
freesasa_results=[]
for o, row in data_for_PDB.iterrows():
    #if there is no PDB
    if row["nmr"]==[]:
        freesasa_results.append("")
        c.append([o,0])
    else:
        #save name of pdb
        if not row["nmr"].lower()+"_1" in pdb_dict.keys():
            freesasa_results.append("")
            c.append([o,2])
            continue
        else:
            filtered_dict = {k:v for k,v in pdb_dict.items() if row["nmr"].lower() in k}
            #if there is no match with phospho site
            temp=[]
            for pdb_name in filtered_dict.keys():
                print(pdb_name)
                resetted=pdb_dict[pdb_name].reset_index()
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


freesasa_results=["" if all(y=="" for y in x) else x for x in freesasa_results]
freesasa_results=[float(i) if ((i!="") & (type(i)!=list)) else i for i in freesasa_results]

data_for_place["nmr result"]=freesasa_results

#%% [visualization and analysis of NMR data]

nmrs_data=data_for_place.loc[data_for_place["nmr result"]!=""]
nmrs_data["nmr result"]=[list(map(float,i))  for i in b["nmr result"]]
nmrs_data=nmrs_data.loc[b["Final Position"]=="core"].reset_index()

# gene_names
a=nmrs_data.uniprot
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

nmrs_data=nmrs_data.merge(rawData, left_on="uniprot", right_on="From", how="left")



# label phosphosite based on nmr result and exclude if no model is in core
for o, row in nmrs_data.iterrows():
    if any([float(y)>5 for y in row["nmr result"]]) and any([float(y)<5 for y in row["nmr result"]]):
        nmrs_data.loc[o,"change"]="Dynamic Core"
    elif not any([float(y)<5 for y in row["nmr result"]]):
        nmrs_data.loc[o,"change"]="Exclude"
    else:
        nmrs_data.loc[o,"change"]="Static Core"

nmrs_data=nmrs_data.loc[nmrs_data.change!="Exclude"]
#look at change amount
for o, row in nmrs_data.iterrows():
    nmrs_data.loc[o,"change amount"]=float(max(row["nmr result"]))-float(min(row["nmr result"]))
    nmrs_data.loc[o,"nmr_mean"]=np.mean(row["nmr result"])
#create identifier for stripplot
nmrs_data["identifier"]=nmrs_data.To+"_"+nmrs_data.residue+"-"+nmrs_data.position.apply(str)
nmrs_data[["nmr result","change","identifier"]].explode("nmr result")


# plot stripplot
sns.set(style="whitegrid", palette="Accent",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(9, 7))
ax=sns.stripplot(x="nmr result", y="identifier", hue="change", data=nmrs_data.explode("nmr result"))
sns.despine(left=True)
ax.set_xlabel("NMR RSA score",size = 20,alpha=0.7)
ax.set_ylabel("Phosphorylation with NMR Models",size = 20,alpha=0.7)
ax.legend(title='', loc="lower right")
plt.savefig("NMR_change_plot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)
# save nmr results
nmrs_data.to_csv("nmr_result.csv", index=False)
