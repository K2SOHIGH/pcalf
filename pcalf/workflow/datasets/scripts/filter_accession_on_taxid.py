
import logging
import pandas as pd
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




def main():
    logging.info("Filter accession based on TaxID")
    logging.info("Loading NCBI reports ....")    
    ncbi = pd.read_csv(str(snakemake.input.report),header=None,sep="\t",low_memory=False)
    ncbi.columns = ["Accession","TaxID"]
    ncbi.set_index("Accession",inplace=True)
    request_taxid = []
    logging.info("Filtering .... ")
    with open(str(snakemake.input.taxids),'r') as fh:
        for l in fh.readlines():
            request_taxid.append(int(l.strip()))
            
    
    subset = ncbi[ (ncbi.TaxID.isin(request_taxid)) ]
    with open(str(snakemake.output), 'w')  as stream:
        for acc in list(subset.index.unique()):
            stream.write("{}\n".format(acc))
    logging.info("Done.")

if __name__ == "__main__":
    main()
