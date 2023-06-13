#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
from gimme_taxa import NCBITaxa, desc_taxa , taxon_info

try:
    log = snakemake.log
    if log:
        log = str(log)    
        os.makedirs(os.path.abspath(os.path.dirname(str(log))),exist_ok=True)
        logging.basicConfig(filename= log , filemode='a', encoding='utf-8', level=logging.DEBUG)
except NameError:
    log = None
    logging.basicConfig(encoding='utf-8', level=logging.DEBUG)





def get_taxid(taxidFH):
    t=[]    
    with open(taxidFH,'r') as fh:
        for l in fh.readlines():            
            t.append(l.strip().split()[1])
    return t


def main():
    # args
    taxid = str(snakemake.params.taxid)
    logging.info("Resolve TaxID children for %s ... " % str(taxid))

    taxidF = str(snakemake.output.taxid)
    taxidname = str(snakemake.output.taxidname)
    taxidFH = open(taxidF,'w')
    nameFH = open(taxidname,'w')
    
    # resolve child taxid / lineage
    ncbi = NCBITaxa(dbfile=None)
    
    taxon_info(taxid,ncbi,nameFH)
    desc_taxa(taxid, ncbi,taxidFH, False)

    taxidFH.close()
    nameFH.close()
    logging.info("End.")

    
if __name__ == "__main__":
	main()
