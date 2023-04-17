import os
import pandas as pd
import requests
import tqdm
import json
import logging

logging.basicConfig(encoding='utf-8', level=logging.INFO)

def download(url,outfile):
    # url = t[0]
    # outfile = t[1]
    try:   
        data = requests.get(url)
        # Save file data to local copy
        if data.status_code==200:
            os.makedirs(os.path.dirname(outfile),exist_ok=True)
            with open(outfile, 'wb') as stream:
                stream.write(data.content)
        return data.status_code            
    except requests.exceptions.MissingSchema:
        logging.info("Skipping invalid url {} ... ".format(url))
        return 500

def report2df(f):
    return pd.read_csv(
        str(f),sep="\t",header=0,index_col=0,low_memory=False)

def make_url(basurl, suffix):
    basurl = basurl.replace("ftp://", "https://")
    return os.path.join(basurl, basurl.split("/")[-1] + "_" + suffix ) 


assembly_dir = str(snakemake.params.assembly_dir)
mapdict = {}
for d in os.listdir(assembly_dir):
    if os.path.isdir(d):
        if d.startswith("GC"):
            mapdict[os.path.basename(d)] = "NCBI Datasets CLI"

ftp_df = report2df(str(snakemake.input.ftp_report))
urls = {}
for acc, url in ftp_df["FTP Path"].items():
    if isinstance(acc,str):
        gid = acc.split(".")[0]
        if gid not in urls:
            urls[gid] = []
        urls[gid].append(url)

accessions = list(report2df(str(snakemake.input.cli_report)).index.unique())
with open(str(snakemake.output.fb),'w') as out:
    with tqdm.tqdm(total=len(accessions)) as pbar:
        while accessions: 
            acc = accessions.pop() # acc : <assembly>.<version>       
            gid = acc.split(".")[0] # gid : <assembly>
            if gid in urls:
                for url in urls[gid]: # get all urls related to <assembly> 
                    #url : <url>/<assembly>.<version>_<name>/
                    acc_from_url = "_".join(url.split("/")[-1].split("_")[0:2]) # <assembly>.<version>
                    # check if assembly dir exists
                    if os.path.isdir(os.path.join(assembly_dir,acc_from_url)):
                        out.write("{}\tNCBI Cli\n".format(acc))
                        continue # Entry have already been downloaded
                    else:
                        gf = os.path.join( assembly_dir , acc , 
                                os.path.basename(url) + "_genomic.fna.gz" )
                        genome_url = make_url(url, "genomic.fna.gz")
                        exc = download(genome_url,gf) # wget genome file
                        if exc > 400: # entry doesn't exists in NCBI ftp server             
                            out.write("{}\tnot available for download\n".format(acc))
                            continue  # continue with next url avalaible for <assembly>
                        else:
                            out.write("{}\tNCBI FTP server\n".format(acc))
                            cf = os.path.join( assembly_dir , acc , 
                                    "cds_from_genomic.fna.gz" )
                            cds_url = make_url(url, "cds_from_genomic.fna.gz") # wget cds if any...
                            exc = download(cds_url, cf)              
            pbar.update(1)  


# for acc, basurl in tqdm.tqdm(ftp_df["FTP Path"].items()):
#     if os.path.isdir( os.path.join(assembly_dir , acc)):
#         # assembly have been download with ncbi CLI dataset tool.
#         continue
#     else:
#         # assembly is missing.
#         # try with wget 
#         gf = os.path.join( assembly_dir , acc , 
#                             os.path.basename(basurl) + "_genomic.fna.gz" )
#         genome_url = make_url(basurl, "genomic.fna.gz")
#         exc = download(genome_url,gf)
#         if exc > 400:
#             mapdict[acc] = "not avalaible for download"
#             continue
#         else:
#             mapdict[acc] = "NCBI FTP server"
        
#         cf = os.path.join( assembly_dir , acc , 
#                             "cds_from_genomic.fna.gz" )
#         cds_url = make_url(basurl, "cds_from_genomic.fna.gz")
#         exc = download(cds_url, cf)    
        
# json.dump(mapdict , open(str(snakemake.output.fb),"w"))