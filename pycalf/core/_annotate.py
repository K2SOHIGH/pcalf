import tempfile
import concurrent.futures
from importlib import resources
import time
import pandas as pd
from pycalf.core import cter 
from pycalf.core import nter 
from pycalf.utils import utils
import logging
import tqdm
import gzip
import itertools

logger = logging.getLogger()
def _annotate(fasta, glyx3_phmm , glyzip_hmms , nterdb_fa , nterdb_tab, args):
    logger.info("Load and sanitize input sequences" )
    # store fasta sequence in easlfasta format into list - https://pyhmmer.readthedocs.io/en/stable/api/easel.html#pyhmmer.easel.Sequence
    raw_sequences = utils.easelfasta(fasta)    
    sequences = {}
    # quick check - remove long sequences
    for i in range(len(raw_sequences)): 
        if len(raw_sequences[i]) > 10000:
            logger.warning("Sequence {} seems very long and will not be scan through hmmer".format(raw_sequences[i].name.decode("UTF-8")))            
            pass
        else:
            sequences[raw_sequences[i].name.decode("ascii")]  = raw_sequences[i]
    logger.info("Number of amino acid sequences :")
    logger.info("   raw file : {}".format(len(raw_sequences)))
    logger.info("   clean file : {}".format(len(sequences)))
   
    ####################################################
    #                   
    #                   SEARCH GLYX3
    #
    ####################################################
    logger.info("Search GlyX3 in sequence database ... ")
    # get hits against GlyX3 HMM profile
    # glyx3hits is a list of Hit object with the following property : seqid, domid, start, end, evalue, coverage
    all_hits_glyx3 = cter.findglyx3(
            sequences,
            glyx3_phmm, 
            glyx3evalue= args.gly3_evalue_threshold, 
            glyx3ievalue = args.gly3_i_evalue_threshold, 
            cpus = args.threads,
            domZ = args.domz
        )
    valid_glyx3_seqs = cter.filtersequences(
        sequences, all_hits_glyx3,
        args.gly3_coverage_threshold
    )
    logger.info("       {} hits against GlyX3 hmm profile.".format(len(valid_glyx3_seqs)))
    logger.info("Done.")
    summary_df = None
    # there is at least one hit against the GlyX3 HMM profile     
    if valid_glyx3_seqs:
        ####################################################
        #                   
        #                   SEARCH GLY1,GLY2,GLY3
        #
        ####################################################
        logger.info("Search glyzip in glyx3+ sequences")
        all_hits_glyzip = cter.findglyzips(
            sequences = valid_glyx3_seqs,
            glyzip_hmms = glyzip_hmms,
            glyzipevalue = args.glyzip_i_evalue, 
            cpus = args.threads,
            domZ = args.domz
        ) 
        logger.info("Done.")

        # write full length sequence with a GlyX3 into a fasta file for blast:     
        tmp = tempfile.NamedTemporaryFile(mode="w+")
        logger.debug("Store sequence with glyx3 in tempfile for blast search.")
        for seqid,seq in valid_glyx3_seqs.items():
            tmp.write(">{}\n{}\n".format(
                seqid,seq.textize().sequence
            ))        
        tmp.flush()    
        ####################################################
        #                   
        #                   SEARCH NTERs
        #
        ####################################################
        logger.info("Search similar N-ter in %s " % args.nterdb_fa )
        all_hits_nter = nter.blastp(                
            #res_dir + "/intermediates/sequence_with_glyX3.fasta",
            tmp.name,
            nterdb_fa,
            nterdb_tab,
            args.nter_evalue,
            blastpexec=args.blastp
        )
        logger.info("done.")
        #summary_df = utils.maketable(all_hits_glyx3 + all_hits_glyzip +  all_hits_nter)
        hits = all_hits_glyx3 + all_hits_glyzip +  all_hits_nter
        for h in hits:            
            h.extract_from_fasta(sequences[h.seqid])
        return hits, sequences 
    return [],[]
        

def annotate_st(*args,**kwargs): 
    return _annotate(*args) , kwargs

def _annotate_mt(p_input):
    #print(p_input)
    args,kwargs = p_input    
    return _annotate(*args) , kwargs

def annotate_mt(*args):
    hits = []
    logger.handlers = [
        h for h in logger.handlers if not isinstance(h, logging.StreamHandler)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        hits += list(tqdm.tqdm(executor.map(_annotate_mt, *args), total=len(*args)))
    
    #finish = time.perf_counter()
    #print(finish)
    return hits#list(itertools.chain.from_iterable(hits))

    
