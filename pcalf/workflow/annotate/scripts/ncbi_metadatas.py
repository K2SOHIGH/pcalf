import sys
import os
import argparse

import logging

import numpy as np
import tqdm
import pandas as pd
import ncbi.datasets.metadata
import ncbi.datasets.metadata.genome as ncbi_genome_mtd

logging.basicConfig()
# see https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/ for a list of attributes

def parse_biosample(query):
    BIOSAMPLE_ATTRIBUTES = {
        "Isolation source" : ["isolate source",
                "isolationsource",
                "isolation_source",
                "source of isolate",
                "isolation_source"],
        "Environment (biome)":["env_broad_scale","biome",
                "broad scale environmental context",
                "env biome",
                "env broad range",
                "env broad scale",
                "environment (biome)",
                "environment biome"],
        "Geographic location": ["geo_loc_name","country",
                "geo loc name",
                "geographic location (country and/or sea region)",
                "geographic location (country and/or sea)",
                "geographic location (country and/or sea, region)",
                "geographic location (country and/or sea,region)",
                "geographic location (country)",
                "geographic location (country:region,area)",
                "geographic location (locality)",
                "geographic location country and or sea",
                "geographic locations",
                "geographic origin",
                "geographical location",
                "geographical location (country:region, location)",
                "geolocname"],
        "Culture collection":["culture_collection",
                "culturecollection"],
        "Collection date" : ["collection_date",
                "colection date",
                "collection date (yyyymmdd)",
                "collection year",
                "collectiondate",
                "date of collection",
                "date sample collected",
                "isolation year",
                "sample collection date",
                "sample date",
                "sampling date",
                "time of sample collection",
                "year isolated"
                ],
        "Sample type" : ["sample_type"],
    }
    
    datas = {}
    
    if "assembly" in query:
        query = query["assembly"]
    
    biosample_attributes = {}
    biosample_accession  = None
    
    if "biosample" in query:
            biosample_accession = query["biosample"]["accession"]
            biosample_attributes = {attr["name"]:attr["value"] for attr in query["biosample"]["attributes"]}
            # flatten biosample attributes

        # given a list of biosample attributes names :
    datas["Biosample"] = biosample_accession
    for label, names in BIOSAMPLE_ATTRIBUTES.items():
        for n in names:            
            if biosample_attributes and n in biosample_attributes:
                datas[label] = biosample_attributes[n]
                break
            else:
                datas[label] = None
    return datas
    
def parse_assembly(query):
    ASSEMBLY_ATTRIBUTES = {
        "Assembly name": ["display_name"],
        "Submitter": ["submitter"],
        "Submission date": ["submission_date"],
    }
    
    datas = {}
    
    if "assembly" in query:
        query = query["assembly"]
    
    for label, names in ASSEMBLY_ATTRIBUTES.items():
        for n in names:            
            if query and n in query:
                datas[label] = query[n]
                break
            else:
                datas[label] = None
    return datas

    
def parse_org(query):
    ORGS_ATTRIBUTES = {
        "Isolate": ["isolate"],
        "TaxID": ["tax_id"],
        "Organism": ["title","sci_name"],
    }
    datas = {}
    
    if "assembly" in query:
        query = query["assembly"]
    
    org_attributes = {}

    if "org" in query:
        org_attributes = query["org"]

    for label, names in ORGS_ATTRIBUTES.items():
        for n in names:            
            if org_attributes and n in org_attributes:
                datas[label] = org_attributes[n]
                break
            else:
                datas[label] = None
    return datas
    
    
def fetch_ncbi_mtds(lacc:list):
    lacc = parse_accession(lacc)    
    
    assembly_datas = { i:{} for i in lacc }
    valid_for_query = [i for i in lacc if i.startswith("GCA_") or i.startswith("GCF_")]
    queries = {}
    
    for query in ncbi_genome_mtd.get_assembly_metadata_by_asm_accessions( valid_for_query ) :            
        query = query["assembly"]        
        if isinstance(
                query,
                ncbi.datasets.openapi.model.v1_assembly_dataset_descriptor.V1AssemblyDatasetDescriptor):            
            query_accession  = query["assembly_accession"]
            queries[query_accession] = query
    for acc , user_acc in lacc.items() :
        query = {}
        if acc in queries:
            query = queries[acc] 
        
        biosample = parse_biosample(query)
        organism  = parse_org(query)
        assembly  = parse_assembly(query)

        assembly_datas[user_acc].update(assembly) 
        assembly_datas[user_acc].update(organism)
        assembly_datas[user_acc].update(biosample)     
        
    return assembly_datas
        



def batch_generator(items, batch_size):
    count = 1
    chunk = []
    for item in items:
        if count % batch_size:
            # fill batch
            chunk.append(item)
        else:
            chunk.append(item)            
            # Batch is full, yield filled batch and clear batch
            yield chunk
            chunk.clear()
        count += 1
    
    # Check that last chunk is not empty and yield it
    if len(chunk):
        yield chunk


def get_args():
    parser = argparse.ArgumentParser(
                    prog='ncbi_metadatas',
                    description='get NCBI metadatas (subset) using the NCBI datasets API'
                    )
    parser.add_argument('accession',help="Either one or more assembly accession separated by a comma (e.g GCF_903937765.1) or a file in one assembly accession per line")   
    parser.add_argument('-o', '--out',default=sys.stdout,help="Output file, default is stdin")     
    args = parser.parse_args()
    return args

def get_snakargs():
    args = argparse.Namespace()
    args.accession = snakemake.input.accession
    args.out = snakemake.output.out
    return args
    
def parse_args():    
    if "snakemake" in globals():
        return get_snakargs()
    else:
        return get_args()


def parse_accession(accs):
    raccs = {}
    for a in accs:
        if a.startswith("GCA_") or a.startswith("GCF_"):
            raccs[ "_".join(a.split("_")[0:2]) ] = a
        else:
            raccs[a] = a
    return raccs

def fetch(accs:dict):
    cpt=0
    while cpt != 10:
        try:
            assembly_datas = {}        
            progesstotal = round(len(accs)/250) 
            if progesstotal == 0:
                progesstotal = 1
            for b in tqdm.tqdm(batch_generator(accs,250),total=progesstotal):                
                res = fetch_ncbi_mtds(b)
                assembly_datas.update(res)
            return assembly_datas
        except:
            logging.info("Error while recovering metadatas from NCBI , let's retry [{}/10]".format(cpt))
            cpt+=1            
    raise NameError("Can't retrieved metadatas from NCBI ")

if __name__ == "__main__":
    args = parse_args()
    accs = []
    if os.path.exists(args.accession):
        with open(args.accession) as fh:
            for line in fh.readlines():
                accs.append(line.strip())
    elif isinstance(args.accession, str):
        accs = [acc.strip() for acc in args.accession.split(",")]
    else:
        logging.error("Wrong type of input")
        get_args()
        exit(-1)
    
    
    assembly_datas = fetch(accs)
    
    df  = pd.DataFrame(assembly_datas).T
    df.fillna(np.NaN, inplace=True)
    df.index.name = "Accession"
    cols = list(df.columns)    
    df[sorted(cols)].to_csv(args.out, sep="\t", header=True, index=True)


