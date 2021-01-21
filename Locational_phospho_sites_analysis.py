
import numpy as np
import pandas as pd
from scipy import stats
import os
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp as multi
from gprofiler import GProfiler
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import math
import urllib.parse
import urllib.request
from gprofiler import GProfiler
import seaborn as sns
import re
from os import listdir
from os.path import isfile, join

#%% [read files]
os.chdir(r"location")

paper_data=pd.read_excel("all_phospho_paper_paper_data_data.xlsx")


#%% [seperation of phospho sites]

cores=paper_data.loc[paper_data["Final Position"]=="core"]
others=paper_data.loc[(paper_data["Final Position"]!="core") & paper_data["Final Position"].notnull()]
interface=paper_data.loc[paper_data["Final Position"]=="Interface"]
intermediate=paper_data.loc[paper_data["Final Position"]=="intermediate"]
surface=paper_data.loc[paper_data["Final Position"]=="surface"]
paper_data=paper_data.loc[paper_data["Final Position"].notnull()]

#%% [statistical analysis]

# checked distribution
stats.kstest(paper_data.functional_score.dropna(),"norm")
stats.kstest(paper_data.functional_score.dropna().apply(np.log2),"norm")

# one-way anova
stats.f_oneway(cores["functional_score"].apply(np.log2),interface["functional_score"].apply(np.log2),\
               intermediate["functional_score"].apply(np.log2),surface["functional_score"].apply(np.log2))

#tukey hsd test
mc1=multi.MultiComparison(paper_data["functional_score"],paper_data["Final Position"])
res1 = mc1.tukeyhsd()
print(res1.summary())

# functional score plot

sns.set(style="whitegrid", palette="deep",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(6, 9))
ax = sns.boxplot(x=paper_data["Final Position"], y=paper_data.functional_score, fliersize=4, notch=True, width=0.5)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7, labelpad=15)
ax.set_ylabel("Functional Score",size = 20,alpha=0.7, labelpad=15)
ax.set(xticklabels=["Intermediate","Surface","Interface","Core"])


#%% [B-factor analysis]

# get list of all pdb ids

pdbs=list(set(paper_data["PDB_ID"].apply(lambda x: x.split(":")[0]).tolist()+paper_data["model PDB_ID"].apply(lambda x: x.split(":")[0]).tolist()))

paper_data=paper_data.copy()
# fill blanks with model pdbs
paper_data["PDB_ID"]=paper_data["PDB_ID"].replace("[]",np.nan)
paper_data["model PDB_ID"]=paper_data["model PDB_ID"].replace("[]",np.nan)

paper_data["PDB_ID"].fillna(paper_data["model PDB_ID"], inplace=True)
# write it out
with open("pdbs.txt", "w+") as f:
    for item in pdbs:
        f.write("%s\n" % item)

#%% [b-factor datas from ignm]

pdbs=pd.read_csv("pdbs.txt",sep=" ",header=None)[0].tolist()

mypath=r"location\Bfactor"
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
# get b-factor datas from ignm site
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


#%% [read b-factor files]

mypath=r"location\Bfactor"
pdbfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
pdbfiles=[x for x in pdbfiles if "bfactor" in x]
pdbfiles=[x for x in pdbfiles if "txt" in x]


# number control function
def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

#read all files and store them in a dictionary

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

# looped through each phospho-site and found their location in bfactor data and obtain x-ray b-factor
bfactor=[]
for o,pdb in paper_data.iterrows():
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

paper_data["x_ray_bfactor"]=bfactor

# looped through each phospho-site and found their location in bfactor data and obtain predicted b-factor

bfactor=[]
for o,pdb in paper_data.iterrows():
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
        if bfile.iloc[found]["GNMbpredicted B- factor"]==0:
            print(o)
            bfactor.append(np.nan)
        elif len(set(bfile["GNMbpredicted B- factor"].tolist()))==1:
            print(o,2)
            bfactor.append(np.nan)
        else:
            bfactor.append(bfile.iloc[found]["GNMbpredicted B- factor"]/np.mean(bfile["GNMbpredicted B- factor"]))

paper_data["gnm_bfactor"]=bfactor

#%% [z-score elimination]
# took z-score for each b-factor and eliminated outliers
z = np.abs(stats.zscore(paper_data["x_ray_bfactor"], nan_policy="omit"))
print(z)
print(np.where(z > 3))
x_ray_wo = paper_data[(z < 3)]

z = np.abs(stats.zscore(paper_data["gnm_bfactor"], nan_policy="omit"))
print(z)
print(np.where(z > 3))
gnm_wo = paper_data[(z < 3)]


x_ray_wo=x_ray_wo.dropna(subset=["x_ray_bfactor"])
gnm_wo=gnm_wo.dropna(subset=["gnm_bfactor"])

#%% [z-score elimination]

#%% [statistical analysis]
x_ray_wo["x_ray_bfactor"].hist()
gnm_wo["x_ray_bfactor"].hist()
x_ray_wo["x_ray_bfactor"].apply(np.log2).hist()
gnm_wo["x_ray_bfactor"].apply(np.log2).hist()

# we continued with x-ray due to their correlated paper_data
df = [np.array(x["x_ray_bfactor"].apply(np.log2)) for _, x in x_ray_wo.groupby('Final Position')]

stats.f_oneway(df[0],df[1],df[2],df[3])

#tukeyhsd
mc1=multi.MultiComparison(x_ray_wo["x_ray_bfactor"].apply(np.log2),x_ray_wo["Final Position"])
res1 = mc1.tukeyhsd()
print(res1.summary())

sns.set(style="whitegrid", palette="deep",font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(6, 9))
ax = sns.boxplot(x=x_ray_wo["Final Position"], y=x_ray_wo.x_ray_bfactor.apply(np.log2),\
                 fliersize=4, notch=True, width=0.5)
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7,labelpad=15)
ax.set_ylabel("log2(Normalized B-factor)",size = 20,alpha=0.7,labelpad=15)
ax.set(xticklabels=["Intermediate","Surface","Interface","Core"])


plt.savefig("mean_x_ray_bfactor_scores_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [cellular location of phospho sites]

#obtain genename from uniprot id by usng uniprot web server
a=paper_data.uniprot

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

paper_data=paper_data.merge(rawData, left_on="uniprot", right_on="From", how="left")

#%%

os.chdir(r"location")

localization=pd.read_csv("proteinatlas_c85499c9.tsv", sep="\t")

places=paper_data.merge(localization[["Uniprot","Subcellular main location"]], left_on="uniprot", right_on="Uniprot", how="left")

places.groupby("Subcellular main location")["Final Position"].count()

cores=places.loc[places["Final Position"]=="core"]
interface=places.loc[places["Final Position"]=="Interface"]
intermediate=places.loc[places["Final Position"]=="intermediate"]
surface=places.loc[places["Final Position"]=="surface"]
# top 10 location in each phospho-site
loc_top_10=cores["Subcellular main location"].value_counts()[0:10].index.tolist()
loc_top_10.extend(interface["Subcellular main location"].value_counts()[0:10].index.tolist())
loc_top_10.extend(intermediate["Subcellular main location"].value_counts()[0:10].index.tolist())
loc_top_10.extend(surface["Subcellular main location"].value_counts()[0:10].index.tolist())

# fill out based on number of each phosphorylation (percentage)
loc_top_10=list(set(loc_top_10))
local=pd.DataFrame(columns=["core","interface","intermediate","surface"], index=loc_top_10)
local["core"]=cores["Subcellular main location"].value_counts()\
    .loc[cores["Subcellular main location"].value_counts().index.intersection(loc_top_10)]\
        /sum(~cores["Subcellular main location"].isnull())*100
local["interface"]=interface["Subcellular main location"].value_counts().loc[interface["Subcellular main location"].value_counts().index.intersection(loc_top_10)]/sum(~interface["Subcellular main location"].isnull())*100
local["intermediate"]=intermediate["Subcellular main location"].value_counts().loc[intermediate["Subcellular main location"].value_counts().index.intersection(loc_top_10)]/sum(~intermediate["Subcellular main location"].isnull())*100
local["surface"]=surface["Subcellular main location"].value_counts().loc[surface["Subcellular main location"].value_counts().index.intersection(loc_top_10)]/sum(~surface["Subcellular main location"].isnull())*100

loc_percen=pd.melt(local.reset_index(), id_vars=["index"])
# bar plot
sns.set(style="whitegrid", palette="deep",font_scale = 1.5,font="Arial")
f, ax = plt.subplots(figsize=(15, 12))
ax=sns.barplot(y="index", x="value", hue="variable", data=loc_percen)
ax.set_xlabel("Percentage",size = 20,alpha=0.7)
ax.set_ylabel("",size = 20,alpha=0.7)
plt.legend(title='', fontsize=18)
plt.savefig("localization.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

# statistical analysis of location
stat=pd.DataFrame(columns=["core","interface","intermediate","surface"], index=loc_top_10)
stat["core"]=cores["Subcellular main location"].value_counts().loc[cores["Subcellular main location"].value_counts().index.intersection(loc_top_10)]
stat["interface"]=interface["Subcellular main location"].value_counts().loc[interface["Subcellular main location"].value_counts().index.intersection(loc_top_10)]
stat["intermediate"]=intermediate["Subcellular main location"].value_counts().loc[intermediate["Subcellular main location"].value_counts().index.intersection(loc_top_10)]
stat["surface"]=surface["Subcellular main location"].value_counts().loc[surface["Subcellular main location"].value_counts().index.intersection(loc_top_10)]


# chi square test for each location
for i in loc_top_10:
    oddsratio, pvalue, x, y = stats.chi2_contingency([np.array(stat.loc[i]),\
                        np.array([sum(~cores["Subcellular main location"].isnull()),\
                                  sum(~interface["Subcellular main location"].isnull()),\
                                   sum(~intermediate["Subcellular main location"].isnull()),\
                                    sum(~surface["Subcellular main location"].isnull())])-np.array(stat.loc[i])])
    if pvalue<0.05:
        print("OddsR: ", oddsratio, "p-Value:", pvalue, "Name:", i)

# relative residue distribution of each phosphorylation types
# residue distribution
x=paper_data.groupby(["residue","Final Position"])["Unnamed: 0"].count()
pd.melt(x.reset_index(), id_vars=["residue"])

#transform data for statistical analysis (chi-square)
x=x.reset_index().pivot(index='residue', columns='Final Position')["Unnamed: 0"]
x2=np.array(x)

stats.chi2_contingency(np.array(x))
#create a relative residue distribution table
relative_table=(x2/(x2.sum(axis=0)))
relative_table.mean(axis=1)
STY=pd.DataFrame(relative_table-np.vstack(relative_table.mean(axis=1)),columns=x.columns, index=x.index).T
STY.to_excel("STY_distribution.xlsx")


# relative bar plot

sns.set_context("paper")
sns.set(font_scale = 2)
fig, ax = plt.subplots()
fig.set_size_inches(11.7, 8.27)
sns.heatmap(STY.set_index("Final Position"), annot=True, annot_kws={"size":18}, linewidths=.6,cmap="YlGnBu", ax=ax)\
    .set_title('Relative Abundance', size=24)
ax.set_ylabel('')

plt.savefig("Figure_2_sty_distribution.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%% [EV mutation score]

#%%

os.chdir(r"location\effects")
# find all mutation tables
mypath=r"location\effects"
mutationfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

# load proteins with phosphorylation sites
mutation={}
for o,row in paper_data.iterrows():
    file=[i for i in mutationfiles if i.split("_")[0]==row.uniprot]
    if len(file)==0:
        pass
    else:
        case=pd.read_csv(file[0], sep=";")
        print(file[0])
        mutation[row.uniprot+"_"+str(row.position)]=case.loc[case.pos==row.position]
# delete empty mutation files
mutation = {k:v for k,v in mutation.items() if len(v)>0}
#get minimum values
mutation = {k:i.loc[i.prediction_epistatic==min(i.prediction_epistatic)].iloc[0] for k,i in mutation.items()}
# match utation scores with phospho-sites
for name,value in mutation.items():
    paper_data.loc[(paper_data.uniprot==name.split("_")[0]) & (paper_data.position==value.pos),"mutation_score"]=value.prediction_epistatic

# distribution check
stats.kstest(paper_data.mutation_score.dropna(),"norm")

# anova
df = [np.array(x["mutation_score"]) for _, x in paper_data.dropna(subset=["mutation_score"]).groupby('Final Position')]

stats.f_oneway(df[0],df[1],df[2],df[3])

# tukey
data_3=paper_data.dropna(subset=["mutation_score"])
mc2=multi.MultiComparison(data_3["mutation_score"],data_3['Final Position'])
res2 = mc2.tukeyhsd()
print(res2.summary())

# sort values based on location
data=paper_data.dropna(subset=["mutation_score"]).sort_values("Final Position")

# violin plot with counts
sns.set(style="whitegrid", palette=sns.color_palette("Set2", 4),font_scale = 1.4,font="Arial")
f, ax = plt.subplots(figsize=(10, 10))
ax = sns.violinplot(x=data["Final Position"], y=data["mutation_score"], inner="quartile")
sns.despine(left=True)
ax.set_xlabel("",size = 20,alpha=0.7, labelpad=15)
ax.set_ylabel("min(EV Mutation Score)",size = 20,alpha=0.7, labelpad=15)
# add counts as text
def add_n_obs(df,group_col,y):
    medians_dict = {grp[0]:grp[1][y].median() for grp in df.groupby(group_col)}
    xticklabels = [x.get_text() for x in plt.gca().get_xticklabels()]
    n_obs = df.groupby(group_col)[y].size().values
    for (x, xticklabel), n_ob in zip(enumerate(xticklabels), n_obs):
        plt.text(x, medians_dict[xticklabel]*1.040, " #obs :\n"+str(n_ob), horizontalalignment='center', fontdict={'size':14}, color='white')

add_n_obs(data,group_col='Final Position',y='mutation_score')
ax.set(xticklabels=["Interface","Core","Intermediate","Surface"])


plt.savefig("EV_mutation_scores_violinplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)

#%%
#protein length
paper_data["protein_length"]=paper_data["Protein Sequence"].apply(len)
cores=paper_data.loc[paper_data["Final Position"]=="core"]
others=paper_data.loc[(paper_data["Final Position"]!="core") & paper_data["Final Position"].notnull()]
interface=paper_data.loc[paper_data["Final Position"]=="Interface"]
intermediate=paper_data.loc[paper_data["Final Position"]=="intermediate"]
surface=paper_data.loc[paper_data["Final Position"]=="surface"]
alls=paper_data.loc[paper_data["Final Position"].notnull()]

no_dublicates=alls.drop_duplicates(["uniprot","Final Position"])

paper_data.groupby("Final Position")["protein_length"].apply(np.mean)

a=cores.drop_duplicates("uniprot")["protein_length"]
b=interface.drop_duplicates("uniprot")["protein_length"]
c=intermediate.drop_duplicates("uniprot")["protein_length"]
d=surface.drop_duplicates("uniprot")["protein_length"]
e=pd.concat([cores.drop_duplicates("uniprot"),interface.drop_duplicates("uniprot")\
             ,intermediate.drop_duplicates("uniprot"),surface.drop_duplicates("uniprot")])

stats.shapiro(no_dublicates.protein_length.apply(np.log2))
stats.kstest(alls.protein_length.dropna(),"norm")
stats.kstest(change_core.functional_score.dropna().apply(np.log2),"norm")

paper_data.protein_length.hist()

stats.f_oneway(np.log2(a),np.log2(b),np.log2(c),np.log2(d))

stats.kruskal(a,b,c,d)

mc1=multi.MultiComparison(np.log2(e["protein_length"]),e["Final Position"])
res1 = mc1.tukeyhsd()
print(res1.summary())

sns.set(style="whitegrid", palette="deep", font_scale=1.2)
f, ax = plt.subplots(figsize=(5, 8))
ax = sns.boxplot(x=e["Final Position"], y=np.log2(e["protein_length"]), fliersize=4, notch=True, width=0.6)
sns.despine(left=True)
ax.set_xlabel("Positions",size = 18,alpha=0.7)
ax.set_ylabel("Log2(Protein Length)",size = 18,alpha=0.7)
ax.set(xticklabels=["Core","Interface","Intermediate","Surface"])

plt.savefig("log2protein_legth_boxplot.pdf",dpi=2000,quality=100,bbox_inches="tight", transparent=True)
