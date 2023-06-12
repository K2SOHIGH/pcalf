import os
import glob
import pandas as pd
import tqdm 

def _parse_biosample_attribute(lk,lv):
    desc = {}
    l = ["isolate source",
            "isolationsource",
            "source of isolate",
            "isolation_source",
            "env_broad_scale",
            "geo_loc_name"
        ]
    for k,v in zip(lk,lv):        
        if k in l:
            desc[k] = v 
    return desc

def parse_ncbi_report(df):
    columns = [
        "Organism Name",
        "Assembly Level",
        "Assembly Name",
        "Assembly Notes",
        "is_mag",
        "ANI Best ANI match ANI",
        "ANI Best ANI match Assembly",
        "ANI Best ANI match Assembly Coverage",
        "ANI Best ANI match Organism",
        "Annotation Release Date",
        "Assembly BioProject Lineage Accession",
        "Assembly BioProject Lineage Title",
        "Assembly BioSample Accession",
        "Assembly BioSample Attribute Name",
        "Assembly BioSample Attribute Value",
        "Assembly BioSample Description Organism Taxonomic ID",
        "Assembly BioSample Owner Name",
        "warning"
    ]
    
    if df.index.name != "Assembly Accession" :
        df.set_index("Assembly Accession",inplace=True)
    datas = {}
    for _ , subdf in df.groupby("Assembly Accession"):
        lk = list(subdf["Assembly BioSample Attribute Name"])
        lv = list(subdf["Assembly BioSample Attribute Value"])
        subdf.drop("Assembly BioSample Attribute Name",axis=1,inplace=True)
        subdf.drop("Assembly BioSample Attribute Value",axis=1,inplace=True)        
        desc = {}
        for c in columns:
            if c in subdf.columns:
                desc[c] = list(subdf[c])[0]
        desc.update(_parse_biosample_attribute(lk,lv))
        datas[_]=desc
    
    return pd.DataFrame(datas).T

def report2df(f):
    return pd.read_csv(
        str(f),sep="\t",header=0,index_col=0,low_memory=False)
# for assembly in str(snakemake.input.assembly_dir):
#     if os.path.isdir(assembly) and assembly.startswith("GC"):


# There is one row per <assembly>.<latest> assembly
# We make a mapping dict by generalizing <assembly>.<latest> to <assembly>
df = report2df(str(snakemake.input.report))
df = parse_ncbi_report(df)
df.index.name = "latest version"
df.reset_index(inplace=True)
df.index = df.apply(lambda x : x["latest version"].split(".")[0],axis=1)
mapper = df.T.to_dict()

datas = {}
rootdir = str(snakemake.params.assembly_dir)
for it in os.scandir(rootdir):
    if it.is_dir():
        acc = os.path.basename(it) #<assembly>.<version>
        gid = acc.split(".")[0]    #<assembly>
        if acc not in datas:
            datas[acc] = {}
        if gid in mapper:
            datas[acc] = mapper[gid] 

unreachable_datas = {}
for acc in df["latest version"].unique():
    if acc not in datas:
        unreachable_datas[acc] = mapper[acc.split(".")[0]]

pd.DataFrame(datas).T.to_csv(
    str(snakemake.output.report) , sep="\t", header=True,index=True
)
pd.DataFrame(unreachable_datas).T.to_csv(
    str(snakemake.output.fail) , sep="\t", header=True,index=True
)