
import logging
from threading import local
import pandas as pd
import sys
import os



try:
    log = snakemake.log
    if log:
        log = str(log)    
        os.makedirs(os.path.abspath(os.path.dirname(str(log))),exist_ok=True)
        logging.basicConfig(filename= log , filemode='a', encoding='utf-8', level=logging.DEBUG)
except NameError:
    log = None
    logging.basicConfig(encoding='utf-8', level=logging.DEBUG)


logging.info("filtering NBCI reports based on taxid")

def main():
    logging.info("Loading NCBI reports ....")    
    ncbi = pd.read_csv(str(snakemake.input.reports),sep="\t",header=0,index_col=0,low_memory=False)
    # request_accession = []
    # with open(str(snakemake.input.accs),'r') as fh:
    #     logging.info(" ..........  Desired accession :")
    #     for l in fh.readlines():
    #         request_accession.append(l.strip().split()[0]) 
    #         logging.info(l.strip().split()[0])

    #print(request_accession)
    # logging.info("Taxonomy based filtering ................... ")
    # request_taxidname = []
    # logging.info(" ..........  Desired taxonomic group :")

    # with open(str(snakemake.input.names),'r') as fh:
    #     for l in fh.readlines():
    #         request_taxidname.append(l.strip().split()[0])
    #         logging.info(l.strip().split()[0])
    #print(request_taxidname)
    request_taxid = []
    with open(str(snakemake.input.taxids),'r') as fh:
        logging.info("TaxID based filtering ................... ")
        logging.info(" ..........  Desired TaxID :")
        for l in fh.readlines():
            request_taxid.append(int(l.strip()))
            logging.info(l.strip())
    
    logging.info("Filtering reports .... ")    
    subset = ncbi[ (ncbi.TaxID.isin(request_taxid)) ]
    with open(str(snakemake.output), 'w')  as stream:
        for acc in list(subset.index.unique()):
            stream.write("{}\n".format(acc))

    # all_request = pd.concat([subset,subset_gcf])
    # logging.info("Remove duplicated entries ...... ")    
    # request_derep = all_request.groupby('Assembly Accession').first()
    # request_derep.to_csv(str(snakemake.output),header=True,index=True,sep='\t')
    

if __name__ == "__main__":
    main()
