
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

#%%
#obtain all features for machine learning system

os.chdir(r"location")
result=pd.read_excel("all_phospho_paper_result_data_with_functional_score.xlsx")
result=result.loc[result["phospho window"].str[6].isin(["S","T","Y"])]
result=result.loc[result["Final Position"].notnull()]



os.chdir(r"location")

change_core=pd.read_excel("core_phosphosites_location_change.xlsx") 
change_core.PDB_filtered=list(map(ast.literal_eval,change_core.PDB_filtered))
change_core=change_core.loc[~change_core["freesasa result"].isnull()]
change_core["freesasa result"]=change_core["freesasa result"].fillna("[]")
change_core["freesasa result"]=list(map(ast.literal_eval,change_core["freesasa result"]))
change_core=change_core.loc[change_core["phospho window"].str[6].isin(["S","T","Y"])]
change_core["just change"]=change_core["just change"].fillna("unchanged")

change_core["freesasa result"] = change_core["freesasa result"].apply(lambda x: [num for num in x if isinstance(num, (int,float))])
change_core=change_core.loc[change_core["freesasa result"].apply(len)>2]

for o, row in change_core.iterrows():
    change_core.loc[o,"change amount"]=float(max(row["freesasa result"]))-float(min(row["freesasa result"]))
    change_core.loc[o,"freesasa_mean"]=np.mean(row["freesasa result"])

change_core['just change'] = change_core['just change'].map({'changed': 'Dynamic\nCore', 'unchanged': 'Static\nCore'})


mutation=change_core.merge(result, how="left", on=["uniprot","position"])
change_core["RSA"]=mutation["freesasa result_y"].values


os.chdir(r"D:\Anaconda\Thesis\Bfactor")

mypath=r"D:\Anaconda\Thesis\Bfactor"     
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "bfactor" in x] 
pdbfiles=[x for x in pdbfiles if "txt" in x] 

                
#read all files and store them in a dataframe
def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

bfactor_dict={}
for out_file in pdbfiles:
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
    look=pd.read_csv(io.StringIO('\n'.join(lines)), delim_whitespace=True, names=["chain identifier", "residue type", "node index", "GNMbpredicted B- factor",  "X-ray crystallographic B factor"])
    bfactor_dict.update({out_file.split("_")[0] : look})

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

bfactor=[]
for o,pdb in change_core.iterrows():
    if (pdb["PDB_ID"].split(":")[0] not in bfactor_dict.keys()):
        bfactor.append(np.nan)
    else:
        bfile=bfactor_dict[pdb["PDB_ID"].split(":")[0]]
        chain=pdb["PDB_ID"].split(":")[1]
        bfile=bfile.loc[bfile["chain identifier"]==chain]
        seq=bfile["residue type"].tolist()
        seqs="".join([d[x] if (x!="UNK" and (x in d.keys())) else "." for x in seq]) 
        found=0
        for i in range(0,7):
            place=seqs.find(pdb["phospho window"][i:i+7])
            if place==-1:
                continue
            else:
                found=place+6-i
                break
        if bfile.iloc[found]["X-ray crystallographic B factor"]==0:
            print(o)
            bfactor.append(np.nan)
        elif len(set(bfile["X-ray crystallographic B factor"].tolist()))==1:
            print(o,2)
            bfactor.append(np.nan)
        else:
            bfactor.append(bfile.iloc[found]["X-ray crystallographic B factor"]/np.mean(bfile["X-ray crystallographic B factor"]))

change_core["x_ray_bfactor"]=bfactor

os.chdir(r"location")

mypath=r"location"     
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "bfactor" in x] 
pdbfiles=[x for x in pdbfiles if "txt" in x] 

                
#read all files and store them in a dataframe
def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


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

change_core["x_ray_bfactor_total"]=bfactor
change_core["x_ray_bfactor_total"]=change_core.x_ray_bfactor_total.apply(np.nanmean)


os.chdir(r"location")

phospho=pd.read_csv("All_cores_secondary.csv")

phospho["uniprot"]=phospho.id.apply(lambda x: x.split("_")[1])

change_core=change_core.merge(phospho, how="left", left_on=["uniprot","position"], right_on=["uniprot","n"])


os.chdir(r"location")

phospho=pd.read_csv("Phosphorylation_site_dataset",sep="\t" ,skiprows=[0,1])
phospho=phospho.loc[phospho.ORGANISM=="human"]
phospho["position"]=phospho.MOD_RSD.apply(lambda x: int(x.split("-")[0][1::]))

known_phos=change_core.merge(phospho, how="left", left_on=["uniprot","position"], right_on=["ACC_ID","position"])

known_phos["DOMAIN_INFO"]=[1 if type(i)==str else 0 for i in known_phos.DOMAIN]


change_core_origin = pd.get_dummies(known_phos.DOMAIN,dtype=int)
change_core_res = pd.get_dummies(change_core.residue,dtype=int)
change_core_res.drop("S",inplace=True, axis=1)


#mutations
os.chdir(r"location\effects")
from os import listdir
from os.path import isfile, join
mypath=r"location\effects"       
pdbfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

mutation={}
for o,row in result.iterrows():
    file=[i for i in pdbfiles if i.split("_")[0]==row.uniprot]
    if len(file)==0:
        pass
    else:
        case=pd.read_csv(file[0], sep=";")
        print(file[0])
        mutation[row.uniprot+"_"+str(row.position)]=case.loc[case.pos==row.position]
        
mutation = {k:v for k,v in mutation.items() if len(v)>0}

mutation_2 = {k:i.loc[i.prediction_epistatic==min(i.prediction_epistatic)].iloc[0] for k,i in mutation.items()}
mutation_3 = {k:i.loc[i.prediction_epistatic==max(i.prediction_epistatic)].iloc[0] for k,i in mutation.items()}
for name,value in mutation_2.items():
    change_core.loc[(change_core.uniprot==name.split("_")[0]) & (change_core.position==value.pos),"min_mutation_score"]=value.prediction_epistatic
for name,value in mutation_3.items():
    change_core.loc[(change_core.uniprot==name.split("_")[0]) & (change_core.position==value.pos),"max_mutation_score"]=value.prediction_epistatic
    
    
#%% organization of features
machine=change_core[change_core.columns[-32:]]
machine=pd.concat([machine, change_core_origin], axis=1, sort=False)
machine=pd.concat([machine, change_core_res], axis=1, sort=False)
machine.drop(columns=["Final Position","intact change",\
                     "freesasa_mean","id","seq","n","rsa","asa","q3","q8","p[q3_H]","p[q3_E]","p[q3_E]"],inplace=True)


machine['just change']=machine['just change'].replace({'Dynamic\nCore':1,'Static\nCore':0})  

machine["function"]=[1 if i>=0.5 else 0 for i in machine.functional_score]
machine.drop(columns=["functional_score"],inplace=True)

Y = machine['function']
X = machine.drop('function', axis=1)

#%% run machine learning system with forest (other methods are implemented beforehand to choose best one)

from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.feature_selection import SelectKBest, chi2, SelectFromModel
from sklearn.svm import LinearSVC
from sklearn.feature_selection import RFE
from sklearn.feature_selection import RFECV
from sklearn.tree import DecisionTreeClassifier
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
import fancyimpute



model_tree = RandomForestClassifier(random_state=100, n_estimators=50)


steps = [('imputation', fancyimpute.IterativeImputer(verbose=0)),
         ('scaler', StandardScaler()),
         ("over", SMOTE(random_state=42)),
         ("under", RandomUnderSampler()),
         ('tree', RandomForestClassifier())]

# Create the pipeline: pipeline 
pipeline = Pipeline(steps)

#best result with selected parameters
# Number of trees in random forest
n_estimators = [1400]
# Number of features to consider at every split
max_features = ['auto']
# Maximum number of levels in tree
max_depth = [130]
# Minimum number of samples required to split a node
min_samples_split = [8]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1]
# Method of selecting samples for training each tree
bootstrap = [False]

"""
#can be tried with these parameters for validity
# Number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 100, stop = 2000, num = 20)]
# Number of features to consider at every split
max_features = ['auto', 'sqrt']
# Maximum number of levels in tree
max_depth = [int(x) for x in np.linspace(10, 150, num = 15)]
max_depth.append(None)
# Minimum number of samples required to split a node
min_samples_split = [int(x) for x in np.linspace(2, 15, num = 14)]
# Minimum number of samples required at each leaf node
min_samples_leaf = [int(x) for x in np.linspace(1, 10, num = 10)]
# Method of selecting samples for training each tree
bootstrap = [True, False]

"over__k_neighbors": [1, 2, 3, 4, 5, 6, 7]
"""

# Create the random grid
paramaters = {'tree__n_estimators': n_estimators,
               'tree__max_features': max_features,
               'tree__max_depth': max_depth,
               'tree__min_samples_split': min_samples_split,
               'tree__min_samples_leaf': min_samples_leaf,
               'tree__bootstrap': bootstrap,
               "over__k_neighbors": [6]}

# Create train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=42)

# Create the GridSearchCV object: gm_cv
gm_cv = GridSearchCV(pipeline, param_grid = paramaters, cv = 10, verbose=2, n_jobs = -1)
# Fit to the training set
gm_cv.fit(X_train, y_train)

y_pred = gm_cv.predict(X_test)
print(classification_report(y_test, y_pred))

y_fit_pred = gm_cv.predict(X_train)
print(classification_report(y_train, y_fit_pred))

# Compute and print the metrics
r2 = gm_cv.score(X_test, y_test)
print("Tuned ElasticNet Alpha: {}".format(gm_cv.best_params_))
print("Tuned ElasticNet R squared: {}".format(r2))


model_tree = RandomForestClassifier(random_state=100, n_estimators=50)


steps = [('imputation', fancyimpute.IterativeImputer(verbose=0)),
         ('scaler', StandardScaler()),
         ("over", SMOTE(random_state=42, k_neighbors=6)),
         ("under", RandomUnderSampler()),
         ('elasticnet', RandomForestClassifier(n_estimators=1400,max_features="auto",max_depth=130,\
                                               min_samples_leaf=1,min_samples_split=8,bootstrap=False))]

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=42)

# Create the pipeline: pipeline 
pipeline = Pipeline(steps)
pipeline.fit(X_train, y_train)

importance=pipeline.steps[4][1].feature_importances_

importance_df=pd.DataFrame(list(zip(X.columns,importance)),columns=["Feature","Value"])
importance_df.sort_values(by="Value",ascending=False,inplace=True)
importance_df.reset_index(drop=True,inplace=True)
importance_df=importance_df.loc[0:25]

#%% feature selection plot


import matplotlib.patches as mpatches


fig = plt.figure(figsize=(10,8))
ax=plt.hlines(y=range(0,len(importance_df.index)), xmin=0, xmax=importance_df['Value'], color='#ffa600',linewidth=3, alpha =0.8)
ax=sns.stripplot(x="Value", y="Feature", data=importance_df, color="#d62728", s=8, alpha=1)
plt.grid()
labels = [item.get_text() for item in ax.get_yticklabels()]
labels[0] = 'RSA change amount'
labels[6] = 'Core status'
labels[3] = 'Disorder Score'

ax.set_yticklabels(labels,fontsize=16, color="darkviolet")
ax.tick_params(axis='x', labelsize=18)
ax.get_yticklabels()[0].set_color("firebrick")
ax.get_yticklabels()[0].set_fontweight("semibold")
ax.get_yticklabels()[0].set_fontsize(16)
ax.get_yticklabels()[6].set_color("firebrick")
ax.get_yticklabels()[6].set_fontweight("semibold")
ax.get_yticklabels()[6].set_fontsize(16)

ax.get_yticklabels()[2].set_color("seagreen")
ax.get_yticklabels()[5].set_color("seagreen")
ax.get_yticklabels()[8].set_color("seagreen")
ax.get_yticklabels()[10].set_color("seagreen")
ax.get_yticklabels()[13].set_color("seagreen")
ax.get_yticklabels()[15].set_color("seagreen")
ax.get_yticklabels()[16].set_color("seagreen")
ax.get_yticklabels()[17].set_color("seagreen")
ax.get_yticklabels()[18].set_color("seagreen")

ax.get_yticklabels()[11].set_color("cornflowerblue")
ax.get_yticklabels()[22].set_color("cornflowerblue")
ax.get_yticklabels()[23].set_color("cornflowerblue")
ax.get_yticklabels()[24].set_color("cornflowerblue")
ax.get_yticklabels()[25].set_color("cornflowerblue")

ax.get_yticklabels()[7].set_color("deeppink")
ax.get_yticklabels()[21].set_color("deeppink")

ax.get_yticklabels()[12].set_color("olive")
ax.get_yticklabels()[19].set_color("olive")


patch2 = mpatches.Patch(color='darkviolet', label='Structural Information')
patch1 = mpatches.Patch(color='firebrick', label='Core Dynamicity')
patch3 = mpatches.Patch(color='seagreen', label='Secondary Structure')
patch4 = mpatches.Patch(color='cornflowerblue', label='Domain')
patch5 = mpatches.Patch(color='deeppink', label='Residue')
patch6 = mpatches.Patch(color='olive', label='Mutation')

ax.legend(handles=[patch1, patch2, patch3, patch4, patch5, patch6],handlelength=0.9,fontsize=14,loc='lower right')

ax.set_ylabel('Features',fontsize=24, labelpad=15)
ax.set_xlabel('Value',fontsize=24, labelpad=15)
ax.set_title('Feature Importance',fontsize=28)

plt.savefig("feature_importance.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

