import re
import os
import logging
import tempfile
import concurrent
import concurrent.futures
import multiprocessing

import pandas as pd
from shutil import which

import pyhmmer
import tqdm

from Bio.SeqFeature import SeqFeature, FeatureLocation

# from utils import log
from .biohmm import Hmm
from .bioseq  import Sequences , Hit

from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s",
    level=logging.INFO,
    handlers=[]
)

def _search_calcyanin(fasta:str, 
                    src:str,
                    glyx3_hmm:pyhmmer.plan7.HMM, 
                    gly1_hmm:pyhmmer.plan7.HMM, 
                    gly2_hmm:pyhmmer.plan7.HMM, 
                    gly3_hmm:pyhmmer.plan7.HMM,
                    nter_fa:str, 
                    domz:int=10000,
                    glyx3_evalue_threshold:float=1e-30,
                    glyx3_coverage_threshold:float=0.5,
                    glyzip_evalue_threshold:float=3.6e-4,
                    glyzip_coverage_threshold:float=0.7,
                    nter_coverage_threshold:float = 0.80, 
                    nter_evalue_threshold:float = 1e-07,
                    tmp_dir = None,
                    ):
    """Find calcyanin within a fasta file using several HMM profiles and 'database' of known N-ter.
    
    Args:
        fasta (str): Path to a fasta file - required.
        src (str): will be attached to calcyanin. Useful when searching across multiple files.  
        glyx3_hmm (pyhmmer.plan7.HMM):  Glyx3 HMM profile - required
        gly1_hmm (pyhmmer.plan7.HMM): Gly1 HMM profile - required
        gly2_hmm (pyhmmer.plan7.HMM): Gly2 HMM profile - required
        gly3_hmm (pyhmmer.plan7.HMM): Gly3 HMM profile - required
        nter_fa (str): Path to fasta file containing known N-ter sequences - required
        glyx3_evalue_threshold (float): E-value threshold for glyx3.
        glyx3_coverage_threshold (float): Coverage theshold for glyx3.
        glyzip_evalue_threshold (float): E-value threshold for glycine zipper 3.6e-4.
        glyzip_coverage_threshold (float): Coverage threshold for glycine zipper 0.7.
        nter_coverage_threshold (float): Coverage threshold for blastp. 
        nter_evalue_threshold (float): E-value threshold for blastp.
        tmp_dir (str): path to a directory where raw results (i.e , before filtering steps) will be stored.

    Raises:
        TypeError if one sequence is not an instance of pyhmmer.easel.DigitalSequence.

    Return:
        Hmm object
    """
    if tmp_dir:
        tmp_dir = os.path.join(tmp_dir, src )
        os.makedirs(tmp_dir, exist_ok=True)

    # Setup sequences - Parse fasta into Sequences object 
    sequences = Sequences(fasta,src=src)
    
        
    # Search sequences against GlyX3 profile
    glyx3_sequences = sequences.hmmsearch([glyx3_hmm])
    glyx3_sequences.to_hits_table().to_csv(os.path.join(tmp_dir,"glyx3.hits.tsv"),sep="\t",header=True,index=True)
    for seq in glyx3_sequences.sequences:
        # We filter hits based on their e-value and the total coverage against the GlyX3 HMM profile.

        hit_below_evalue_threshold = [ h for h in seq.hits if h.score < glyx3_evalue_threshold]     
        if hit_below_evalue_threshold:          
            start_glyx3 = min([h.start_location for h in hit_below_evalue_threshold])
            end_glyx3 = max([h.stop_location for h in hit_below_evalue_threshold])
            scores = [h.score for h in hit_below_evalue_threshold]
            hmm_profile_len = [h.target_len for h in hit_below_evalue_threshold][0]
            coverage_on_glyx3 = (end_glyx3-start_glyx3)/hmm_profile_len #glyx3_easel_hmm[0].M
            mean_identity = sum([h.identity for h in hit_below_evalue_threshold])/len(hit_below_evalue_threshold)

            # COVERAGE THESHOLD 
            if coverage_on_glyx3 > glyx3_coverage_threshold:

                f = SeqFeature(
                    FeatureLocation(
                        start_glyx3, 
                        end_glyx3
                        ),
                    id="GlyX3",
                    type="domain", 
                    qualifiers={
                        "n_hits":len(hit_below_evalue_threshold),
                        "score":min(scores),
                        "src":"GlyX3_HMM_Profile",
                        "src_len":hmm_profile_len,
                        "identity":mean_identity
                    })
                # Clean sequence hits list for next step.
                seq.hits=[]
                seq.features.append( f )

    # we free some memory.
    del sequences 
    # keep only sequence with a hit against the glyx3 profile with a E-value below and a coverage above E-value and coverage thresholds respectively.
    glyx3_sequences.sequences = [s for s in glyx3_sequences.sequences if s.features] 
    # We can then use those sequences for the next steps aka Glycine zipper annotation and N-ter comparison.
    
    # Search GlyZip profiles
    glzips_easel_hmm = [gly1_hmm,gly2_hmm,gly3_hmm]
    glyx3_sequences = glyx3_sequences.hmmsearch(glzips_easel_hmm)
    glyx3_sequences.to_hits_table().to_csv(os.path.join(tmp_dir,"glyzips.hits.tsv"),sep="\t",header=True,index=True)
    for seq in glyx3_sequences.sequences:
        # Hit filtering and features creation.
        seq.per_residue_annotation()
        current = Hit()
        features = []
        features_len = []
        fraglen = 0
        for res,hit in seq.res.items():   
            if hit:
                # E-VALUE THRESHOLD
                if hit[0].score < glyzip_evalue_threshold:        
                    if current.hid == hit[0].hid:
                        fraglen += 1
                        continue    
                    else:
                        if current.hid:                            
                            features.append(current)
                            features_len.append(fraglen)
                        current = hit[0]           
                        fraglen=0                               

            else:
                if current.hid:
                    features.append(current)
                    features_len.append(fraglen)
                current = Hit()
                fraglen=0
        
        # COVERAGE THRESHOLD 
        feat = [ f for f,fl in zip(features,features_len)  if fl/f.target_len > glyzip_coverage_threshold ]
        # We can then remove all hits that do not respect the coverage threshold.
        seq.hits = feat
        # For each hit we create a feature 
        for valid_glyzip_hit in seq.hits:
            f = SeqFeature(
                FeatureLocation(
                    valid_glyzip_hit.start_location, 
                    valid_glyzip_hit.stop_location
                    ),
                id=valid_glyzip_hit.hid,
                type="domain", 
                qualifiers={
                    "n_hits":1,
                    "score":valid_glyzip_hit.score,
                    "src":"{}_HMM_Profile".format(valid_glyzip_hit.hid),
                    "src_len":valid_glyzip_hit.target_len,
                    "identity":valid_glyzip_hit.identity
                })

            seq.features.append( f )
        # Again, we clean the hits attributes for the next step.
        seq.hits = []

    # Search N-ter
    blastp_raw_results = glyx3_sequences.blastp(nter_fa)
    if  blastp_raw_results is not None:
        if not blastp_raw_results.empty:
            glyx3_sequences.to_hits_table().to_csv(os.path.join(tmp_dir,"nter.hits.tsv"),sep="\t",header=True,index=True)

    for seq in glyx3_sequences.sequences:
        for valid_nter_hit in seq.hits:
            # E-value and coverage filtering.
            if valid_nter_hit.score < nter_evalue_threshold and \
                valid_nter_hit.coverage >= nter_coverage_threshold and \
                valid_nter_hit.hid.split("||")[-1] != seq.id:
                f = SeqFeature(
                    FeatureLocation(
                        valid_nter_hit.start_location, 
                        valid_nter_hit.stop_location
                        ),
                    id="N-ter",
                    type="domain", 
                    qualifiers={
                        "n_hits":len(seq.hits),
                        "score":valid_nter_hit.score,
                        "src":valid_nter_hit.hid,
                        "src_len":valid_nter_hit.target_len,
                        "identity":valid_nter_hit.identity
                    })                
                seq.features.append( f )

                # we break the loop because we keep only the best hit - aka the nearest neighboor.
                break
        seq.hits = []
    # we return the calcyanin as a Sequences object with all their features.
    return glyx3_sequences

def _search_calcyanin_mt(args):    
    """helper function to run search_calcyanin() from concurrent futures."""    
    # mem = (psutil.Process().memory_info().rss / (1024 * 1024))/1000
    # logging.debug("memory usage : {}".format(mem))
    return _search_calcyanin(*args)

def _search_calcyanin_concurrent(fastas:list,  
                    srcs:list,                                
                    glyx3_hmm:Hmm, 
                    gly1_hmm:Hmm, 
                    gly2_hmm:Hmm, 
                    gly3_hmm:Hmm,
                    nter_fa:str, 
                    domz:int=10000,
                    glyx3_evalue_threshold:float=1e-30,
                    glyx3_coverage_threshold:float=0.5,
                    glyzip_evalue_threshold:float=3.6e-4,
                    glyzip_coverage_threshold:float=0.7,
                    nter_coverage_threshold:float = 0.80, 
                    nter_evalue_threshold:float = 1e-07,
                    tmp_dir = None):
    """Wrapper to run multiple file at once with the function _search_calcyanin.
    
    Args:
        fastas (list): list of fasta file - required.
        srcs (list): list of flag that will be attached to calcyanin. Useful when searching across multiple files.  
        glyx3_hmm (pyhmmer.plan7.HMM):  Glyx3 HMM profile - required
        gly1_hmm (pyhmmer.plan7.HMM): Gly1 HMM profile - required
        gly2_hmm (pyhmmer.plan7.HMM): Gly2 HMM profile - required
        gly3_hmm (pyhmmer.plan7.HMM): Gly3 HMM profile - required
        nter_fa (str): Path to fasta file containing known N-ter sequences - required

        glyx3_evalue_threshold (float): E-value threshold for glyx3 
        glyx3_coverage_threshold (float): Coverage theshold for glyx3

    Raises:
        AssertionError if one hmm is not an instance of Hmm.

    Return:
        Hmm object
    """                 
                       
    assert isinstance(glyx3_hmm,Hmm)
    assert isinstance(gly1_hmm, Hmm)
    assert isinstance(gly2_hmm, Hmm)
    assert isinstance(gly3_hmm, Hmm)
    assert len(srcs)==len(fastas)

    pools = []
    l = []
    for fasta_file,src in zip(fastas,srcs):
        pools.append(
            (fasta_file,
            src,
            glyx3_hmm.hmm,
            gly1_hmm.hmm,
            gly2_hmm.hmm,
            gly3_hmm.hmm,
            nter_fa,
            domz,
            glyx3_evalue_threshold,
            glyx3_coverage_threshold,
            glyzip_evalue_threshold,
            glyzip_coverage_threshold,
            nter_coverage_threshold, 
            nter_evalue_threshold,
            tmp_dir
            )
        )
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()-1) as executor:
        
        l+=list(tqdm.tqdm(executor.map(_search_calcyanin_mt, pools), total=len(pools) , colour="yellow") )
        
    calcyanin_sequences = Sequences()
    for e in l:
        calcyanin_sequences.sequences += e.sequences
    return calcyanin_sequences

def __modifyline__(string,val):
    """Modify the weight of Glycine residue"""
    liste = string.split()
    if len(liste) >= 21:
        liste[6] = str(val*float(liste[6]))
    elif len(liste) == 20:
        liste[5] = str(val*float(liste[5]))
    else : 
        return string
    return("  ".join(liste))

def increase_glycine_weight(hmm,pc):
    """Increase weight of Glycine within a hmm profile.
    
    Args:
        hmm (Hmm): Hmm object.
        pc (float): weigth factor.
    Raises:
        None
    Return:
        weighted pyhmmer.plan7.HMM
    """
    # dump hmm as tmpfile
    tmp = tempfile.NamedTemporaryFile(mode="wb")
    hmm.hmm.write(tmp,binary=False)
    tmp.flush()
    # update lines 
    start =  True
    val = 1 - pc
    l = []
    with open(tmp.name,"r") as fin:
        for line in fin:
            if start:
                res = re.match(".*COMPO.*",line)
                if res:
                    start = False
                    newline = __modifyline__(line,val)+"\n"
                else:
                    newline = line
            else:
                    newline = __modifyline__(line,val)+"\n"
            l.append(newline)    
    tmp = tempfile.NamedTemporaryFile(mode="w")
    with open(tmp.name,'wb') as fout:
        for _ in l:
            tmp.write(_)
    tmp.flush()
    
    with pyhmmer.plan7.HMMFile( tmp.name ) as hmm_file:
        hmm = hmm_file.read()
    return hmm

def update_hmm(l_seq,hmm,is_update_iterative=True):
    """Update Hmm by aligning sequences against it
    
    Args:
        l_seq (list): List of Seq object.
        hmm (Hmm): Hmm object.
    Raises:
        None
    Return:
        Updated Hmm object.
    """    

    l_digital_seq = [ s.digitize() for s in l_seq ]
    if is_update_iterative:
        new_hmm = hmm.iterative_align(l_digital_seq)
    else:        
        new_hmm = hmm.hmmalign( l_digital_seq )
    
    new_hmm.hmm = increase_glycine_weight( new_hmm,0.2 )
    return new_hmm

def parse_nterdb(nterdb):   
    """parse the nterdb file to dictionnary"""
    nter_dict={} 
    with open(nterdb) as db:
        for line in db.readlines():
            nter, strain, seq = line.strip().split()
            nter_dict[strain] = (nter,seq)
            
    return nter_dict

def search(
        out:str,
        fastas:list, 
        srcs:list,
        glyx3:Hmm, 
        gly1:Hmm, 
        gly2:Hmm, 
        gly3:Hmm,
        nterdb:str, 
        domz:int=10000,
        glyx3_evalue_threshold:float=1e-30,  
        glyx3_coverage_threshold:float=0.60,  
        glyzip_evalue_threshold:float=3.6e-4,
        glyzip_coverage_threshold:float=0.7,
        nter_coverage_threshold:float = 0.80, 
        nter_evalue_threshold:float = 1e-07,
        max_iteration:int=1,
        is_update_iterative:bool=False): 
    """Search for calcyanin within one or more fasta files using several HMM profiles and 'database' of known N-ter.

    If is_iterative is set , then HMMs profiles and N-ter database is updated at each iteration with new calcyanin if any.
    it's first search for calcyanin based on the Glycine zipper repeat (aka Glyx3). Then each glycine zipper is searched using 
    a specific HMM profile. Finally, sequences with a Glyx3 are searched against a database of known N-ter using blastp. 
    
    Args:
        fastas (list): list of fasta file - required.
        srcs (list): list of flag that will be attached to calcyanin. Useful when searching across multiple files. 
        glyx3_hmm (pyhmmer.plan7.HMM):  Glyx3 HMM profile - required
        gly1_hmm (pyhmmer.plan7.HMM): Gly1 HMM profile - required
        gly2_hmm (pyhmmer.plan7.HMM): Gly2 HMM profile - required
        gly3_hmm (pyhmmer.plan7.HMM): Gly3 HMM profile - required
        nter_fa (str): Path to fasta file containing known N-ter sequences - required
        glyx3_evalue_threshold (float): E-value threshold for glyx3 
        glyx3_coverage_threshold (float): Coverage theshold for glyx3
        glyzip_evalue_threshold (float): E-value theshold for glycine zipper
        glyzip_coverage_threshold (float): Coverage thrshold for glycine zipper
        max_iteration (int): Number of iteration 
        is_update_iterative (bool): Update hmm sequence by sequence.
    Raises:
        TypeError if one sequence is not an instance of pyhmmer.easel.DigitalSequence.
    Return:
        A tuple containing calcyanins as Sequences object, Glyx3, Gly1, Gly2 and Gly3 as Hmm objects.
    """

    
    not_converged = True
    logging.info("Max Iteration : {}".format(max_iteration))
    
    ite=0    
    seq_per_iteration = {}
    logging.info("Number of sequences in the N-ter database at start : {}".format(len(nterdb)))

    while not_converged:
        tmp_dir = os.path.join(out,"tmp","iteration_{}".format(ite))
        logging.info("Dump Nter DB as fasta.")
        nterfa = tempfile.NamedTemporaryFile(mode="w+")
        for sid, nter in nterdb.items():
            nterfa.write(">{}||{}\n{}\n".format(nter[0],sid,nter[1]))
        nterfa.flush()
        # search calcyanin        
        logging.info("[iteration {}] |  Search calcyanin".format(ite))
        calcyanin_sequences = _search_calcyanin_concurrent(fastas, 
                srcs,
                glyx3,
                gly1,
                gly2,
                gly3,
                nterfa.name,
                domz,
                glyx3_evalue_threshold,
                glyx3_coverage_threshold,
                glyzip_evalue_threshold,
                glyzip_coverage_threshold,
                nter_coverage_threshold, 
                nter_evalue_threshold,
                tmp_dir
            )
        not_converged=False       
        new_calc = 0
        new_seq = 0
        valid_calcyanin = Sequences() 
        
        if calcyanin_sequences.sequences:            
            for seq in calcyanin_sequences.sequences:     
                if seq.id not in seq_per_iteration:                                              
                    seq_per_iteration[seq.id] = ite
                    new_seq+=1                    
                    nter = ",".join([f.features[0].qualifiers["src"].split("|")[0] for f in  seq.get_feature("N-ter")])
                    cter = ",".join([f.features[0].id for f in seq.get_feature("Gly1" )+
                                    seq.get_feature("Gly2" )+
                                    seq.get_feature("Gly3" )])
                    flag = decision_tree(nter,cter)
                    if flag == "Calcyanin with known N-ter" or flag == "Calcyanin with new N-ter":
                        new_calc += 1 
                        valid_calcyanin.sequences.append(seq)
                        not_converged=True                                        

        logging.info("[iteration {}] | New sequences with a match against GlyX3 [+{}] ! :)".format(ite , new_seq))                    
        logging.info("[iteration {}] | New calcyanin [+{}] ! :)".format(ite , new_calc))                            
        if not_converged:
            logging.info("Convergence not reached. Updating HMM profiles and N-ter DB.")
            # logging.info("[iteration {}] | New calcyanin found [+{}] ! :)".format(ite , new_calc))                    
            glyx3features = valid_calcyanin.get_feature("GlyX3")
            if glyx3features:  
                logging.info("[iteration {}] |  Updating GlyX3".format(ite))
                glyx3 = update_hmm(glyx3features,glyx3,is_update_iterative)

            
            gly1features = valid_calcyanin.get_feature("Gly1")        
            if gly1features:
                logging.info("[iteration {}] |  Updating Gly1".format(ite))
                gly1 = update_hmm(gly1features,gly1,is_update_iterative)
  
            
            gly2features = valid_calcyanin.get_feature("Gly2")        
            if gly2features:
                logging.info("[iteration {}] |  Updating Gly2".format(ite))
                gly2 = update_hmm(gly2features,gly2,is_update_iterative)

            
            gly3features = valid_calcyanin.get_feature("Gly3")        
            if gly3features:
                logging.info("[iteration {}] |  Updating Gly3".format(ite))
                gly3 = update_hmm(gly1features,gly3,is_update_iterative)

            
            
            nterfeatures = valid_calcyanin.get_feature("N-ter")
            if nterfeatures:                
                logging.info("[iteration {}] |  Updating N-ter".format(ite))
                for nf in tqdm.tqdm(nterfeatures,colour="yellow"):
                    ntype,src = nf.features[0].qualifiers["src"].split("||")
                    idt = nf.features[0].qualifiers["identity"]                    
                    if idt<100:                        
                        nterdb.update(
                            {nf.id:( ntype , str(nf.seq) )}
                        )
                    else:
                        logging.warning("Sequence already in DB ...")                     
            logging.info("Number of sequences in the N-ter database after update : {}".format(len(nterdb)))                    
            ite += 1
            if ite == max_iteration:
                logging.warning("Max iteration reached.")
                not_converged=False
        else:
            logging.info("[iteration {}] | No new calcyanin detected.".format(ite))
            logging.info("[iteration {}] | Stop.".format(ite))
    
    return calcyanin_sequences, glyx3 , gly1, gly2, gly3 , nterdb , seq_per_iteration

            
def decision_tree(nter , cter):
    """give a flag to a calcyanin based on its C-ter modular organization and its N-ter."""
    if re.search("Gly1,Gly2,Gly3",cter):
        if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
            flag = "Calcyanin with known N-ter"
        else:
            flag = "Calcyanin with new N-ter"
    elif re.search("Gly1,Gly3",cter) and nter == "Y-type":
        flag = "Calcyanin with known N-ter"
    elif re.search ("Gly(1|2|3)",cter):
        if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
            flag = "Atypical Gly region with known N-ter"
        else:
            flag = "Atypical Gly region with new N-ter"
    else:
        flag="Ancestral gly containing protein"
    return flag

def valid_calcyanin(sequences):
    """Subset 
    
    Args:
        l_seq (list): List of Seq object.
        hmm (Hmm): Hmm object.
    Raises:
        None
    Return:
        Updated Hmm object.
    """    


def parse_nterdb(nterdb):   
    """parse the nterdb file to dictionnary"""
    nter_dict={} 
    with open(nterdb) as db:
        for line in db.readlines():
            nter, strain, seq = line.strip().split()
            nter_dict[strain] = (nter,seq)
            
    return nter_dict
