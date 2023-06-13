#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import logging

import pandas as pd
try:
    log = snakemake.log
    if log:
        log = str(log)    
        os.makedirs(os.path.abspath(os.path.dirname(str(log))),exist_ok=True)
        logging.basicConfig(filename= log , filemode='a', encoding='utf-8', level=logging.DEBUG)
except NameError:
    log = None
    logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

logging.info("Retrieving NCBI reports ... ")



coltypes = {
    "Assembly Accession":"str",
    "bioproject":"str",
    "SubGroup":"str",
    "biosample":"str",
    "TaxID":"int",
    "#Organism/Name":"str",
    "Release Date":"str",
    "Modify Date":"str",
    "Status":"str",
    "FTP Path":"str",
    "Strain":"str", 
    "gbrs_paired_asm":"str",	
    "paired_asm_comp":"str" 
}

def download_prokaryotes_report(report="https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt",subgroup=None,taxid=None):
    """
        report = ncbi prokaryotic report
        subgroup = 
    """
    logging.info("Prokaryotes report to dataframe(%s) ...." % report)    

    df = pd.read_csv(report,sep="\t",header=0,low_memory=False)
    df = df[["Assembly Accession","BioProject Accession","BioSample Accession","TaxID","SubGroup","#Organism/Name","Release Date","Modify Date","Status","FTP Path","Strain"]]
    df.columns = ["Assembly Accession","bioproject","biosample","TaxID",'SubGroup',"#Organism/Name","Release Date","Modify Date","Status","FTP Path","Strain"]    

    df["TaxID"] = df.TaxID.fillna(0)
    #df["TaxID"] = df.TaxID.astype(int)
    df.set_index("Assembly Accession",inplace=True)
    logging.info("Download prokaryotes.txt report finished ...." )    

    return df

def download_refseq_report(report="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"):
    logging.info("RefSeq assembly report to dataframe (%s) ...." % report)    

    df = pd.read_csv(report,sep="\t",dtype='str',skiprows=0,header=1,low_memory=False)
    df = df[["#assembly_accession","bioproject","biosample","species_taxid","organism_name","seq_rel_date", "assembly_level","ftp_path", "infraspecific_name", "gbrs_paired_asm",	"paired_asm_comp"]]
    df.columns = ["Assembly Accession","bioproject","biosample","TaxID","#Organism/Name","Release Date","Status","FTP Path","Strain", "gbrs_paired_asm",	"paired_asm_comp"]    
    df["TaxID"] = df.TaxID.fillna(0)
    #df["TaxID_bis"] = df.TaxID_bis.astype(int)
    df.set_index("Assembly Accession",inplace=True)
    # logging.info("Download RefSeq assembly report finished ...." )    
    return df

def download_genbank_report(report="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"):
    
    logging.info("GenBank assembly report to dataframe (%s) ...." % report)    
    df = pd.read_csv(report,sep="\t",dtype='str',skiprows=0,header=1,low_memory=False)
    df = df[["#assembly_accession","bioproject","biosample","species_taxid","organism_name","seq_rel_date", "assembly_level","ftp_path", "infraspecific_name", "gbrs_paired_asm",	"paired_asm_comp"]]
    df.columns = ["Assembly Accession","bioproject","biosample","TaxID","#Organism/Name","Release Date","Status","FTP Path","Strain", "gbrs_paired_asm",	"paired_asm_comp"]      
    df["TaxID"] = df.TaxID.fillna(0)
    #df["TaxID"] = df.TaxID_bis.astype(int)
    df.set_index("Assembly Accession",inplace=True)
    # logging.info("Download GenBank assembly report finished ...." )    

    return df



def merge_report():
    logging.info("Download NCBI reports  ...." )    

    genbank = download_genbank_report()
    refseq = download_refseq_report()
    prok = download_prokaryotes_report()
    # 1) merge prok and genbank
    logging.info("Merge NCBI reports  ...." )    
    df = pd.concat([prok,genbank],axis=0,sort=False)
    df.groupby("Assembly Accession").first() # merge row with identical accession  - priority to prokaryotes.txt
    # 2) merge refseq
    df = pd.concat([df,refseq],axis=0,sort=False)
    df = df.reset_index().astype(coltypes)
    logging.info("Remove duplicated entries  ...." )    

    derep = df.sort_values("FTP Path",ascending=False).groupby('Assembly Accession').first()
    
    return derep
    
def get_snakargs():    
    return str(snakemake.output)

def get_args():
    if "snakemake" in globals():
        return get_snakargs()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('out', type=str,
			help='Output file')   
    return str(parser.parse_args().out)



    

def main():
    out = get_args()
    # get reports from NCBI ftp server
    df = merge_report()
    df.to_csv(out,sep="\t",header=True,index=True)
    logging.info("Download NCBI reports finished  ...." )    



    
if __name__ == "__main__":
	main()
