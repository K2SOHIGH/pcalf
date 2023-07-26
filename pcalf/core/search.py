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

from .biohmm import Hmm
from .bioseq  import Sequences , Hit, Seq


def glyzips_to_features(seq,
                        glyzip_evalue_threshold,
                        glyzip_coverage_threshold,
                        ):

    seq.per_residue_annotation() # for each residue we create an empty list
    
    current = None
    glyzip_hits = []
    glyzip_frag_len = []        
    fraglen = 0
    for _,hits in seq.res.items():
        # return a tuple (pos,res)
        # none, one or multiple hit exists for this residue, best hit is at the the top of the list.        
        if hits:                
            # if hit(s) exists , keep only glyzip                
            hits = [h for h in hits if h.hid in ["Gly1","Gly2","Gly3"]]
            if hits:
                # E-VALUE THRESHOLD          
                if hits[0].score <= glyzip_evalue_threshold:        
                    #The hit is below the E-value threshold, we keep it.
                    if current and current.hid == hits[0].hid: 
                        # The hit is the same as before, so we continue
                        fraglen += 1
                        continue    
                    else: 
                        # The hit is different from before, 
                        # therefore we store residues that share the 
                        # same annotation as a new feature.
                        if current:           
                            fraglen+=1                                              
                            glyzip_hits.append(current)
                            glyzip_frag_len.append(fraglen)                    
                        # feature have changed.
                        current = hits[0]           
                        fraglen=0
            # this residue do not have any hit. 
            else:
                if current:
                    fraglen+=1
                    glyzip_hits.append(current)
                    glyzip_frag_len.append(fraglen)
                current = None
    # finally we add the final feature.
    if current:
        fraglen+=1                                              
        glyzip_hits.append(current)
        glyzip_frag_len.append(fraglen)
    

    # COVERAGE THRESHOLD
    valid_glyzip_hit = [ f for f,fl in zip(glyzip_hits,glyzip_frag_len)  if  fl/f.target_len >= glyzip_coverage_threshold ]        


    # For each hit we create a feature 
    for hit in valid_glyzip_hit:  
        if hit.score:
            f = SeqFeature(
            FeatureLocation(
                hit.start_location, 
                hit.stop_location
                ),
            id=hit.hid,
            type="domain", 
            qualifiers={
                "n_hits"  :1,
                "score"   :  "{:.2e}".format(hit.score),
                "src":"{}_HMM_Profile".format(hit.hid),
                "src_len" : hit.target_len,
                "coverage": hit.coverage,
                "identity": hit.identity
            })
            seq.features.append( f )
    return seq


def glyx3_to_feature(
    seq:Seq,
    glyx3_evalue_threshold:float=1e-30,
    glyx3_coverage_threshold:float=0.3,
    ):
    # We filter hits based on their e-value and the total coverage against the GlyX3 HMM profile.    
    is_glyx3_feature_exist = False
    hit_below_evalue_threshold = [ h for h in seq.hits if h.score <= glyx3_evalue_threshold and h.hid =="Glyx3"]             
    if hit_below_evalue_threshold:          
        start_glyx3 = min([h.start_location for h in hit_below_evalue_threshold])
        end_glyx3 = max([h.stop_location for h in hit_below_evalue_threshold])
        scores = min([str("{:.2e}".format(h.score)) for h in hit_below_evalue_threshold])            
        coverage_on_glyx3 = (end_glyx3-start_glyx3)/hit_below_evalue_threshold[0].target_len
        identity = min([str(h.identity) for h in hit_below_evalue_threshold])
        range = ",".join([str(h.start_location)+":"+str(h.stop_location) for h in hit_below_evalue_threshold])
        # COVERAGE THESHOLD 
        if coverage_on_glyx3 >= glyx3_coverage_threshold:
            f = SeqFeature(
                FeatureLocation(
                    start_glyx3, 
                    end_glyx3
                    ),
                id="GlyX3",
                type="domain", 
                qualifiers={
                    "n_hits":len(hit_below_evalue_threshold),
                    "range":range,
                    "score":scores,
                    "src":"GlyX3_HMM_Profile",
                    "src_len":hit_below_evalue_threshold[0].target_len,
                    "coverage":coverage_on_glyx3,             
                    "identity":identity
                })
            seq.features.append( f )   
            is_glyx3_feature_exist = True 
    return seq

def nter_to_feature(
    seq,
    nter_evalue_threshold:float = 1e-07,
    nter_coverage_threshold:float = 0.80,
    ):
    self_match = None
    valid_match = []
    nter_hits = [h for h in seq.hits if h.method=="blastp"]
    if nter_hits:
        for hit in nter_hits:
            # E-value and coverage filtering.
            if hit.score < nter_evalue_threshold and \
                hit.coverage >= nter_coverage_threshold:

                f = SeqFeature(
                    FeatureLocation(
                        hit.start_location, 
                        hit.stop_location
                        ),
                    id="N-ter",
                    type="domain", 
                    qualifiers={
                        "n_hits":len(seq.hits),
                        "score":hit.score,
                        "src":hit.hid,
                        "src_len":hit.target_len,
                        "coverage":hit.coverage,
                        "identity":hit.identity
                    })     
                if hit.hid.split("||")[-1] == seq.id:
                    self_match = f
                else:
                    valid_match.append(f)                                                        
        if valid_match:
            seq.features.append(valid_match[0]) # because we keep only the best hit - aka the nearest neighboor.
        elif self_match:
            seq.features.append(self_match) # because there is no nearest neighboor but a self match exist.
        else:
            pass # because there is not hit at all.
    return seq

def hits_2_features(
            sequences,
            glyx3_evalue_threshold:float=1e-30,
            glyx3_coverage_threshold:float=0.3,
            glyzip_evalue_threshold:float=1e-5,
            glyzip_coverage_threshold:float=0.7,
            nter_coverage_threshold:float = 0.80, 
            nter_evalue_threshold:float = 1e-07):
    """
    Args:
        glyx3_evalue_threshold (float): E-value threshold for glyx3.
        glyx3_coverage_threshold (float): Coverage theshold for glyx3.
        glyzip_evalue_threshold (float): E-value threshold for glycine zipper 3.6e-4.
        glyzip_coverage_threshold (float): Coverage threshold for glycine zipper 0.7.
        nter_coverage_threshold (float): Coverage threshold for blastp. 
        nter_evalue_threshold (float): E-value threshold for blastp.
    """

    #logging.info("Filtering GlyX3 hits and convert to feature.")
    for seq in sequences.sequences:
        glyx3_to_feature(
            seq,
            glyx3_evalue_threshold,
            glyx3_coverage_threshold
        )
        #logging.debug("Filtering glyzip for {}.".format(seq.id))
        # Hit filtering and features creation.
        # init per_residue_annotation 
        glyzips_to_features(
            seq,
            glyzip_evalue_threshold,
            glyzip_coverage_threshold
        )

        #logging.debug("Resolve N-ter nearest neighboor for {}.".format(seq.id))        
        nter_to_feature(
            seq,
            nter_evalue_threshold,
            nter_coverage_threshold
        )
    return sequences

def search_calcyanin(fasta:str, 
            src:str,
            glyx3_hmm:pyhmmer.plan7.HMM, 
            gly1_hmm:pyhmmer.plan7.HMM, 
            gly2_hmm:pyhmmer.plan7.HMM, 
            gly3_hmm:pyhmmer.plan7.HMM,
            nter_fa:str='', 
            Z:int=10000,
            domZ:int=10000,
            ):
    """Find calcyanin within a fasta file using several HMM 
    profiles and 'database' of known N-ter.
    
    Args:
        fasta (str): Path to a fasta file - required.
        src (str): will be attached to calcyanin. Useful when searching across multiple files.  
        glyx3_hmm (pyhmmer.plan7.HMM):  Glyx3 HMM profile - required
        gly1_hmm (pyhmmer.plan7.HMM): Gly1 HMM profile - required
        gly2_hmm (pyhmmer.plan7.HMM): Gly2 HMM profile - required
        gly3_hmm (pyhmmer.plan7.HMM): Gly3 HMM profile - required
        nter_fa (str): Path to fasta file containing known N-ter sequences - required
    Raises:
        TypeError if one sequence is not an instance of pyhmmer.easel.DigitalSequence.

    Return:
        Hmm object
    """

    # Setup sequences - Parse fasta into Sequences object 
    sequences = Sequences(fasta,src=src)
    
    # 1) Search sequences against GlyX3 profile
    logging.debug("1) Search {} against GlyX3".format(fasta))
    logging.debug("[{}] Number of seq : {}".format(glyx3_hmm.name.decode(), glyx3_hmm.nseq ))
    logging.debug("[{}] Number of nodes : {}".format(glyx3_hmm.name.decode(), glyx3_hmm.M ))
    
    # .hmmsearch() return only sequences with at least one hit.
    sequences_out_glyx3_search = sequences.hmmsearch([glyx3_hmm],Z=Z, domZ=domZ )
    logging.debug("Number of sequences out glyx3 search : {}".format(
        len(sequences_out_glyx3_search.sequences)
    ))
    # We will only keep sequence with a glyX3 hit, thus we can free some memory.
    del sequences 
    sequences_out_glyx3_search.sequences = [s for s in sequences_out_glyx3_search.sequences if s.hits]    
    logging.debug("[GlyX3] number of sequence with a hit against GlyX3  : {}".format(len(sequences_out_glyx3_search.sequences)))

    # We can then use those sequences for the next steps 
    # i.e, Glycine zipper annotation and N-ter comparison.
    # 2) Search GlyZip profiles
    logging.debug("2) Search {} sequences against Glyzip HMM profiles".format(len(sequences_out_glyx3_search.sequences)))
    glzips_easel_hmm = [gly1_hmm,gly2_hmm,gly3_hmm]
    for _ in glzips_easel_hmm:
        logging.debug("[{}] Number of seq : {}".format(_.name.decode(), _.nseq ))
        logging.debug("[{}] Number of nodes : {}".format(_.name.decode(), _.M ))    
    
    sequences_out_glyx3_search.hmmsearch(glzips_easel_hmm, Z=Z ,domZ=domZ)

    # 3) Compare sequence with known N-ter.
    if nter_fa:
        logging.debug("3) blast {} sequences against known N-ters".format(len(sequences_out_glyx3_search.sequences)))
        sequences_out_glyx3_search.blastp(nter_fa)
    else:
        logging.debug("3) Skip n-ter detection as no fasta file is provided.")
    return sequences_out_glyx3_search

def _search_calcyanin_mt(args):    
    """helper function to run search_calcyanin() from concurrent futures."""    
    # mem = (psutil.Process().memory_info().rss / (1024 * 1024))/1000
    # logging.debug("memory usage : {}".format(mem))
    return search_calcyanin(*args) # return a Sequences object # bioseq.py

def search_calcyanin_concurrent(fastas:list,  
                    srcs:list,                                
                    glyx3_hmm:Hmm, 
                    gly1_hmm:Hmm, 
                    gly2_hmm:Hmm, 
                    gly3_hmm:Hmm,
                    nter_fa:str, 
                    Z:int=10000,
                    domZ:int=10000
            ):
    """Wrapper to run multiple file at once with the function _search_calcyanin.
    
    Args:
        fastas (list): list of fasta file - required.
        srcs (list): list of flag that will be attached to calcyanin. Useful when searching across multiple files.  
        glyx3_hmm (pyhmmer.plan7.HMM):  Glyx3 HMM profile - required
        gly1_hmm (pyhmmer.plan7.HMM): Gly1 HMM profile - required
        gly2_hmm (pyhmmer.plan7.HMM): Gly2 HMM profile - required
        gly3_hmm (pyhmmer.plan7.HMM): Gly3 HMM profile - required
        nter_fa (str): Path to fasta file containing known N-ter sequences - required

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
            Z,
            domZ
            )
        )
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()-1) as executor:
        
        l+=list(tqdm.tqdm(executor.map(_search_calcyanin_mt, pools), total=len(pools) , colour="CYAN") )
    
    # concatenate Sequences objects
    sequences = Sequences()
    for e in l:
        sequences.sequences += e.sequences
    return sequences

def _modifyline(string,val):
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
                    newline = _modifyline(line,val)+"\n"
                else:
                    newline = line
            else:
                    newline = _modifyline(line,val)+"\n"
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

    l_digital_seq = [ s.digitize() for s in l_seq]# if s.id not in sequences_in_msa]
    
    logging.info("Number of sequence provided to update_hmm() : {} ".format(
        len(l_digital_seq)))
    logging.info("Number of sequences and nodes in HMM : {}, {}".format(
        hmm.hmm.nseq,hmm.hmm.M))
    new_hmm = hmm.hmmalign(l_digital_seq , iterative = is_update_iterative)
    logging.info("Number of sequences and nodes in updated HMM : {}, {}".format(
        new_hmm.hmm.nseq,new_hmm.hmm.M))
    new_hmm.hmm = increase_glycine_weight( new_hmm , 0.2 )

    return new_hmm

def parse_nterdb(nterdb):   
    """parse the nterdb file to dictionnary"""
    nter_dict={} 
    with open(nterdb) as db:
        for line in db.readlines():
            nter, strain, seq = line.strip().split()
            nter_dict[strain] = (nter,seq)
            
    return nter_dict

def get_coverage_and_evalue_threshold(features_table,feature_src):
    max_e_value = 0
    min_coverage = 1
    for _,row in features_table[features_table.feature_src.isin(feature_src)].iterrows():
        e_value = max([float(i) for i in str(row['e-value']).split(';')])
        coverage  = min([float(i) for i in str(row['coverage']).split(';')])
        if coverage < min_coverage:
            min_coverage=coverage
        if e_value > max_e_value:
            max_e_value = e_value   
    return float('{:3e}'.format(max_e_value*10)),round(min_coverage-0.1,1)

def auto_thresholds(reference,glyx3,gly1,gly2,gly3,Z,domZ):    
    sequences = search_calcyanin(
        reference,
        '-',     
        glyx3.hmm,
        gly1.hmm,
        gly2.hmm,
        gly3.hmm,
        Z=Z,
        domZ=domZ,    
    )    
    hits_2_features(
        sequences,
        glyx3_evalue_threshold=1e-3,
        glyx3_coverage_threshold=0.5,
        glyzip_evalue_threshold=1e-3,
        glyzip_coverage_threshold=0.5,
        nter_coverage_threshold = 0, 
        nter_evalue_threshold = 1
    )
    features_table = sequences.to_feature_table()
    glyx3_evalue_threshold, glyx3_coverage_threshold = get_coverage_and_evalue_threshold(
        features_table,['GlyX3_HMM_Profile'])
    glyzip_evalue_threshold, glyzip_coverage_threshold = get_coverage_and_evalue_threshold(
        features_table,['Gly1_HMM_Profile','Gly2_HMM_Profile','Gly3_HMM_Profile'])

    return (glyx3_evalue_threshold, glyx3_coverage_threshold, glyzip_evalue_threshold, glyzip_coverage_threshold)


def pcalf(        
        fastas:list, 
        srcs:list,
        glyx3:Hmm, 
        gly1:Hmm, 
        gly2:Hmm, 
        gly3:Hmm,
        nterdb:str, 
        Z:int=10000,
        domZ:int=10000,
        glyx3_evalue_threshold:float=None,  
        glyx3_coverage_threshold:float=None,  
        glyzip_evalue_threshold:float=None,
        glyzip_coverage_threshold:float=None,
        nter_coverage_threshold:float = 0.80, 
        nter_evalue_threshold:float = 1e-4,
        is_update_iterative:bool=False,
        **kwargs
        ): 
    """Search for calcyanin within one or more fasta files using several HMM profiles and 'database' of known N-ter.
    
    it's first search for calcyanin based on the Glycine zipper repeat (aka Glyx3). Then each glycine zipper is searched using 
    a specific HMM profile. Finally, sequences with a Glyx3 are searched against a database of known N-ter using blastp. 
    
    Args:
        fastas (list): list of fasta file - required.
        srcs (list): list of flag that will be attached to calcyanin. Useful when searching across multiple files. 
        glyx3 (BioHMM.Hmm): Glyx3 HMM profile - required
        gly1  (BioHMM.Hmm): Gly1 HMM profile - required
        gly2  (BioHMM.Hmm): Gly2 HMM profile - required
        gly3  (BioHMM.Hmm): Gly3 HMM profile - required
        nterdb (str): Path to fasta file containing known N-ter sequences - required
        glyx3_evalue_threshold (float): E-value threshold for glyx3 
        glyx3_coverage_threshold (float): Coverage theshold for glyx3
        glyzip_evalue_threshold (float): E-value theshold for glycine zipper        
        glyzip_coverage_threshold (float): Coverage theshold for glycine zipper        
        max_iteration (int): Number of iteration 
        is_update_iterative (bool): Update hmm sequence by sequence.
    Raises:
        TypeError if one sequence is not an instance of pyhmmer.easel.DigitalSequence.
    Return:
        A tuple containing calcyanins as Sequences object, Glyx3, Gly1, Gly2 and Gly3 as Hmm objects.
    """

    logging.info('Get E-value and coverage thresholds based on input HMM profiles.')
    reference_file = tempfile.NamedTemporaryFile(mode="w")    
    glyx3.dump(reference_file.name)    
    reference_file.flush()
    auto_glyx3_e_value,auto_glyx3_cov,auto_glyzip_e_value,auto_glyzip_cov =  auto_thresholds(
        reference_file.name,
        glyx3,
        gly1,
        gly2,
        gly3,
        Z,
        domZ
    )

    glyx3_evalue_threshold    = glyx3_evalue_threshold    if glyx3_evalue_threshold    else auto_glyx3_e_value  
    glyx3_coverage_threshold  = glyx3_coverage_threshold  if glyx3_coverage_threshold  else auto_glyx3_cov
    glyzip_evalue_threshold   = glyzip_evalue_threshold   if glyzip_evalue_threshold   else auto_glyzip_e_value  
    glyzip_coverage_threshold = glyzip_coverage_threshold if glyzip_coverage_threshold else auto_glyzip_cov

    logging.info("""\nConfiguration:
    GlyX3 E-value threshold: {} [auto: {}]
    GlyX3 coverage threshold: {} [auto: {}]
    GlyZip E-value threshold: {} [auto: {}]
    GlyZip coverage threshold: {} [auto: {}]
    blast E-value threshold for N-ter extremity: {}
    blast coverage threshold for N-ter extremity: {}
    """.format(
        glyx3_evalue_threshold, auto_glyx3_e_value,
        glyx3_coverage_threshold, auto_glyx3_cov,
        glyzip_evalue_threshold,  auto_glyzip_e_value,
        glyzip_coverage_threshold, auto_glyzip_cov,
        nter_evalue_threshold, nter_coverage_threshold
    ))

    for profile in [glyx3,gly1,gly2,gly3]:
        logging.info('{}:  Number of sequences = {}.'.format(
            profile.hmm.name.decode(),profile.hmm.nseq))
    
    logging.info("Dump Nter DB as fasta.")
    logging.debug("Number of sequences in the N-ter database at start : {}".format(len(nterdb)))
    nterfa = tempfile.NamedTemporaryFile(mode="w+")
    for sid, nter in nterdb.items():
        nterfa.write(">{}||{}\n{}\n".format(nter[0],sid,nter[1]))
    nterfa.flush()
    # search calcyanin        
    
    logging.info("Start search for {} files.".format(len(fastas)))    
    sequences = search_calcyanin_concurrent(fastas, 
            srcs,
            glyx3,
            gly1,
            gly2,
            gly3,
            nterfa.name,
            Z,
            domZ
        )
    logging.info("Total number of sequences with at least one GlyX3 hit: {}".format(len(sequences.sequences)))    

    logging.info("Start features creation.")
    sequences = hits_2_features(
        sequences,
        glyx3_evalue_threshold,  
        glyx3_coverage_threshold,  
        glyzip_evalue_threshold,
        glyzip_coverage_threshold,
        nter_coverage_threshold, 
        nter_evalue_threshold) 
    

   
    new_calc = 0
    new_seq = 0
    valid_calcyanin = Sequences() 
    if sequences.sequences:            
        for seq in sequences.sequences:        
            if seq.get_feature("GlyX3" ):
                new_seq += 1           
                nter = ",".join([f.features[0].qualifiers["src"].split("|")[0] for f in  seq.get_feature("N-ter")])
                cter = ",".join([f.features[0].id for f in seq.get_feature("Gly1" )+
                                    seq.get_feature("Gly2" )+
                                    seq.get_feature("Gly3" )])
                flag = decision_tree(nter,cter)
                if flag == "Calcyanin with known N-ter" or flag == "Calcyanin with new N-ter":
                    new_calc += 1 
                    valid_calcyanin.sequences.append(seq)
            

        logging.info("New sequences with a match against GlyX3 [+{}]".format(new_seq))                    
        logging.info("New calcyanin [+{}] ! :)".format(new_calc))                            
        
        logging.info("Updating HMM profiles and N-ter DB.")
        
        glyx3features = valid_calcyanin.get_feature("GlyX3")
        if glyx3features:  
            logging.info("Updating GlyX3")
            glyx3 = update_hmm(glyx3features,glyx3,is_update_iterative)

            
        gly1features = valid_calcyanin.get_feature("Gly1")        
        if gly1features:    
            logging.info("Updating Gly1")
            gly1 = update_hmm(gly1features,gly1,is_update_iterative)

            
        gly2features = valid_calcyanin.get_feature("Gly2")        
        if gly2features:
            logging.info("Updating Gly2")
            gly2 = update_hmm(gly2features,gly2,is_update_iterative)

            
        gly3features = valid_calcyanin.get_feature("Gly3")        
        if gly3features:
            logging.info("Updating Gly3")
            gly3 = update_hmm(gly1features,gly3,is_update_iterative)

        nterfeatures = valid_calcyanin.get_feature("N-ter")
        if nterfeatures:                
            logging.info("Updating N-ter")
            for nf in nterfeatures:
                ntype,src = nf.features[0].qualifiers["src"].split("||")
                idt = nf.features[0].qualifiers["identity"]                    
                if idt<100:                        
                    nterdb.update(
                        {nf.id:( ntype , str(nf.seq) )}
                    )
                else:
                    logging.debug("Already in DB ... {} has 100% identity with {} which is already in the N-ter DB.".format(                            
                        nf.id,
                        src,
                    ))                     
            logging.info("Number of sequences in the N-ter database after update : {}".format(len(nterdb)))                    
    return sequences, glyx3 , gly1, gly2, gly3 , nterdb 
          


def iterative_pcalf(*args,**kwargs):
    ite = 0
    max_ite = kwargs.get('max_iteration')
    logging.info('Iterative search enable [max-iteration: {}]'.format(max_ite))
    keep_going = True
    fastas,names,glyx3,gly1,gly2,gly3,nterdb = args
    glyx3_initial_size = glyx3.hmm.nseq
    while keep_going:
        keep_going = False
        # do 
        new_sequences, glyx3 , gly1, gly2, gly3 , nterdb  = pcalf(
            fastas,names,glyx3,gly1,gly2,gly3,nterdb,**kwargs)
        if glyx3.hmm.nseq > glyx3_initial_size:
            keep_going = True
            glyx3_initial_size = glyx3.hmm.nseq

    ite+=1
    logging.info('{}/{} iterations done.'.format(ite,max_ite))
    if ite==max_ite:
        keep_going = False
    logging.info('Converged (i.e , no new sequence detected) in {} iterations [max-iteration: {}]'.format(ite,max_ite))
    return new_sequences, glyx3 , gly1, gly2, gly3 , nterdb



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

def parse_nterdb(nterdb):   
    """parse the nterdb file to dictionnary"""
    nter_dict={} 
    with open(nterdb) as db:
        for line in db.readlines():
            nter, strain, seq = line.strip().split()
            nter_dict[strain] = (nter,seq)
            
    return nter_dict
