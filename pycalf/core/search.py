
import sys
import logging
import concurrent
import concurrent.futures
from shutil import which

if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

import pyhmmer
import tqdm

from Bio.SeqFeature import SeqFeature, FeatureLocation

from pycalf.utils import log
from pycalf.core.Hmm import Hmm
from pycalf.core.Sequences  import Sequences, Seq , Hit




logging.basicConfig(
    format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s",
    level=logging.INFO,
    handlers=[]
)



def _search_calcyanin(fasta:str, 
                    src:str,
                    glyx3_easel_hmm:pyhmmer.plan7.HMM, 
                    gly1_hmm:pyhmmer.plan7.HMM, 
                    gly2_hmm:pyhmmer.plan7.HMM, 
                    gly3_hmm:pyhmmer.plan7.HMM,
                    nter_fa:str, 
                    glyx3_evalue_threshold:float=1e-30,
                    glyx3_coverage_threshold:float=0.5,
                    ):

    # Setup sequences
    sequences = Sequences(fasta,src=src)

    # Search GlyX3
    glyx3_sequences = sequences.hmmsearch([glyx3_easel_hmm])
    for seq in glyx3_sequences.sequences:
        hit_below_evalue_threshold = [ h for h in seq.hits if h.score < glyx3_evalue_threshold]     
        if hit_below_evalue_threshold:          
            start_glyx3 = min([h.start_location for h in hit_below_evalue_threshold])
            end_glyx3 = max([h.stop_location for h in hit_below_evalue_threshold])
            scores = [h.score for h in hit_below_evalue_threshold]
            hmm_profile_len = [h.target_len for h in hit_below_evalue_threshold][0]
            coverage_on_glyx3 = (end_glyx3-start_glyx3)/hmm_profile_len #glyx3_easel_hmm[0].M
            mean_identity = sum([h.identity for h in hit_below_evalue_threshold])/len(hit_below_evalue_threshold)
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
    
    glyx3_sequences.sequences = [s for s in glyx3_sequences.sequences if s.features] 
    
    # Search GlyZip

    glzips_easel_hmm = [gly1_hmm,gly2_hmm,gly3_hmm]
    glyx3_sequences = glyx3_sequences.hmmsearch(glzips_easel_hmm)
    
    for seq in glyx3_sequences.sequences:
        seq.per_residue_annotation()
        current = Hit()
        features = []
        features_len = []
        fraglen = 0
        for res,hit in seq.res.items():   
            if hit:
                if hit[0].score < 1:        
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
        
        glyzip_coverage_threshold = 0.8
        feat = [ f for f,fl in zip(features,features_len)  if fl/f.target_len > glyzip_coverage_threshold ]
        seq.hits = feat
    
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
        seq.hits = []


    # Search N-ter
    nter_evalue_t = 1e-4
    nter_cov_t = 0.7
    glyx3_sequences.blastp(nter_fa)
    for seq in glyx3_sequences.sequences:
        for valid_nter_hit in seq.hits:
            if valid_nter_hit.score < nter_evalue_t and valid_nter_hit.coverage >= nter_cov_t:
                f = SeqFeature(
                    FeatureLocation(
                        valid_nter_hit.start_location, 
                        valid_nter_hit.stop_location
                        ),
                    id=valid_nter_hit.hid,
                    type="domain", 
                    qualifiers={
                        "n_hits":len(seq.hits),
                        "score":valid_nter_hit.score,
                        "src":"blastp",
                        "src_len":valid_nter_hit.target_len,
                        "identity":valid_nter_hit.identity
                    })
                seq.features.append( f )
                # because we keep only the best one - aka the nearest neighboor.
                break
        seq.hits = []
    return glyx3_sequences

def _search_calcyanin_mt(args):    
    return _search_calcyanin(*args)

def search_calcyanin(fastas:list,  
                    srcs:list,                                
                    glyx3_hmm:pyhmmer.plan7.HMM, 
                    gly1_hmm:pyhmmer.plan7.HMM, 
                    gly2_hmm:pyhmmer.plan7.HMM, 
                    gly3_hmm:pyhmmer.plan7.HMM,
                    nter_fa:str, 
                    glyx3_evalue_threshold:float=1e-30,
                    glyx3_coverage_threshold:float=0.5,   
                    ):
    assert len(srcs)==len(fastas)
    pools = []
    l = []
    for fasta_file,src in zip(fastas,srcs):
        pools.append(
            (fasta_file,
            src,
            glyx3_hmm,
            gly1_hmm,
            gly2_hmm,
            gly3_hmm,
            nter_fa,
            glyx3_evalue_threshold,
            glyx3_coverage_threshold,
            )
        )
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        l+=list(tqdm.tqdm(executor.map(_search_calcyanin_mt, pools), total=len(pools)))
        
    calcyanin_sequences = Sequences()
    for e in l:
        calcyanin_sequences.sequences += e.sequences

    return calcyanin_sequences

def iterative_search(fastas:list, 
                    srcs:list,
                    glyx3:Hmm, 
                    gly1:Hmm, 
                    gly2:Hmm, 
                    gly3:Hmm,
                    nter_fa:str, 
                    glyx3_evalue_threshold:float=1e-30,
                    glyx3_coverage_threshold:float=0.5,                
                    ): 
    logging.info("Iterative : True")
    not_converged = True
    lseqs = []
    max_iteration = 10
    logging.info("Max Iterative : {}".format(max_iteration))
    ite=0
    while not_converged:
        # search calcyanin
        
        ite += 1
        if ite == max_iteration:
            logging.warning("Max iteration reached.")
            break
        logging.info(" [iteration {}] |  Search calcyanin".format(ite))
        calcyanin_sequences = search_calcyanin(fastas, 
                srcs,
                glyx3.hmm,
                gly1.hmm,
                gly2.hmm,
                gly3.hmm,
                nter_fa,
                glyx3_evalue_threshold,
                glyx3_coverage_threshold 
            )
        
        if calcyanin_sequences.sequences:
            
            for seq in calcyanin_sequences.sequences:
                not_converged=False
                if seq.id not in lseqs:
                    logging.info(" [iteration {}] |  New calcyanin detected.".format(ite))
                    lseqs.append(seq.id)
                    not_converged=True
        
            logging.info(" [iteration {}] |  Updating GlyX3".format(ite))
            glyx3features = calcyanin_sequences.get_feature("GlyX3")
            if glyx3features:        
                glyx3_hmm,glyx3_msa = glyx3.hmmalign(glyx3features)
                glyx3.hmm = glyx3_hmm
                glyx3.msa = glyx3_msa
            
            logging.info(" [iteration {}] |  Updating Gly1".format(ite))
            gly1features = calcyanin_sequences.get_feature("Gly1")        
            if gly1features:
                gly1_hmm,gly1_msa = gly1.hmmalign(gly1features)
                gly1.hmm = gly1_hmm
                gly1.msa = gly1_msa

            logging.info(" [iteration {}] |  Updating Gly2".format(ite))
            gly2features= calcyanin_sequences.get_feature("Gly2")   
            if gly2features:     
                gly2_hmm,gly2_msa = gly2.hmmalign(gly2features)
                gly2.hmm = gly2_hmm
                gly2.msa = gly2_msa

            logging.info(" [iteration {}] |  Updating Gly3".format(ite))
            gly3features = calcyanin_sequences.get_feature("Gly3")        
            if gly3features:
                gly3_hmm,gly3_msa = gly3.hmmalign(gly3features)
                gly3.hmm = gly3_hmm
                gly3.msa = gly3_msa
    return calcyanin_sequences, glyx3 , gly1, gly2, gly3    

            
def to_fasta(filehandle, sequences):
    assert not isinstance(filehandle,str), "Expect file like object not str"
    SeqIO.write(sequences, filehandle,format="fasta")






if __name__ == "__main__":
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(log.CustomFormatter())
    logging.getLogger('').addHandler(
        console)
    #f = "/Users/mmillet/harley_databases/tmp/datas/taxids/1126/ncbi_dataset/data/GCA_008757435.1/cds_from_genomic.faa.gz"
    glyx3 = Hmm("../datas/GlyX3.msa.fa")
    gly1 = Hmm("../datas/Gly1.msa.fa")
    gly2 = Hmm("../datas/Gly2.msa.fa")
    gly3 = Hmm("../datas/Gly3.msa.fa")
    nter = "../datas/nterdb.fasta"

    fastas = []
    names = []
    with open("test.txt") as stream:
        for line in stream.readlines():
            n,f = line.strip().split()
            fastas.append(f)
            names.append(n)
    logging.info("Search start.")
    fastas = ["../test/ncbi_dataset/data/GCF_000307995.1/protein.faa"]
    names = ["GCF_000307995"]
    search_calcyanin(
        fastas,#[0:5],
        names,#[0:5],
        glyx3,
        gly1,
        gly2,
        gly3,
        nter,        
        glyx3_evalue_threshold=1e-30)
    logging.info(" End. ")
    #to_fasta(sys.stdout,sequences = calc.sequences)
    #print(Sequence(s).digitize())
    

    #print(seqs.digitize()[0:110])