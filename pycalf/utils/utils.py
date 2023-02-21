
import os
import tempfile
import pandas as pd
import glob
import pyhmmer.easel
import yaml
from utils import log
import logging
import numpy as np

#logger = log.setlogger(__name__)
logger = logging.getLogger()
    
def parse_nterdb(nterdb):   
    """parse the nterdb file to dictionnary"""
    nter_dict={} 
    with open(nterdb) as db:
        for line in db.readlines():
            nter, strain, seq = line.strip().split()
            nter_dict[strain] = (nter,seq)
            
    return nter_dict


# def args2dict(args):
#     config = {}
#     for arg in vars(args):
#         config[arg] = getattr(args,arg)
#     return config


# def input_from_dir(input , extension):
#     if input:
#         if os.path.isdir(input):
#             #IDS, = glob_wildcards(input + "/{id}." + extension)
#             IDS = glob.glob(input + "/*" + extension)
#             return {
#                 os.path.basename(i).split(extension)[0] : i for i in IDS
#                 }

# def input_from_yaml(input):
#     if input:
#         if os.path.exists(input):
#             conff = open(input)
#             datas = yaml.load(conff,Loader=yaml.FullLoader)
#             return datas
#         else:
#             msg="""WORKFLOW INPUT : {} not found
#             """
#             logger.error(msg.format(input))
#             raise FileNotFoundError(msg.format(input))
#     else:
#         return None

# def input_is_fasta(input):
#     if input:
#         if os.path.exists(input):
#             n = os.path.basename(input)
#             return {n : os.path.abspath(input)}


# def parse_input(input , extension):
#     if input:
#         if os.path.isdir(input):
#             return input_from_dir(input, extension)
#         elif os.path.isfile(input):
#             if input.endswith(".yaml"):
#                 return input_from_yaml(input)
#             return input_is_fasta(input)
#         else:
#             return None
#     return None

# def module(snakefile,cprefix,args):

#     logger.info("Snakemake will install conda environment in %s" % cprefix)

#     configfile = tempfile.NamedTemporaryFile(mode="w+")
#     CONFIG = args2dict(args)
#     yaml.dump( CONFIG, configfile )
#     configfile.flush()

#     cmd = """
#         snakemake --snakefile  {snakefile} -j{threads} --rerun-triggers mtime --use-conda --configfile {config} --conda-prefix {cp} {snakargs}
#     """.format(
#         snakefile = snakefile ,
#         threads = args.threads,
#         cp = cprefix,
#         config = configfile.name,
#         snakargs = args.snakargs
#     )

#     logger.info("""running :
#         %s """ % cmd )
#     excode = os.system(cmd)

#     if excode != 0:
#         logger.error("Hum ... something went wrong while executing the workflow ... :( ")
#         exit(-1)
#     logger.info("Great , pycalf-collection finished without error ! :)")
#     return excode


# def easelhmm(f:str):
#     with pyhmmer.plan7.HMMFile(f) as hmm_file:
#         hmm = hmm_file.read()
#     return hmm


# def easelfasta(f,digital:bool=True):
#     if isinstance(f,str):
#         if not os.path.exists(f):
#             logger.error("{} file not found".format(f))
#             exit(-1)
#     else:
#         tmp = tempfile.NamedTemporaryFile(mode="w+")
#         tmp.write(f.read())                      
#         tmp.flush()
#         f = tmp.name                
#         logger.info("create temporary file {}".format(f))
#     try:
#         with pyhmmer.easel.SequenceFile(f, digital=digital, ignore_gaps=True) as seqs_file:
#             return list(seqs_file)
#     except ValueError:
#         # ValueError: Could not determine alphabet of file: /file
#         with pyhmmer.easel.SequenceFile(f, digital=False, ignore_gaps=True ) as seqs_file:
#             return [ s.digitize(pyhmmer.easel.Alphabet.amino())  for s in list(seqs_file)]


# def getseqbyname(name, sequences:list):
#     for i in sequences:
#         if i.name.decode("ascii") == name:
#             return i
#     raise ValueError("Sequence not found : %s" % name)


# def pyhmmsearch(sequences:list,hmms:list,cpus=4,**kwargs):
#     assert isinstance(hmms,list)
#     assert isinstance(sequences,list)    
    
#     for hmm in hmms:
#         assert isinstance(hmm, pyhmmer.plan7.HMM)
#     for seq in sequences:
#         assert isinstance(seq, pyhmmer.easel.DigitalSequence)

#     all_hits = list(pyhmmer.hmmsearch(hmms,sequences, cpus=4 , **kwargs))
#     return all_hits


# def targetview(tophits):
#     """
#         tophits : list of query hmm tophits
#     """
#     targets = {}
#     for queryhmmhits in tophits:
#         for hit in queryhmmhits:
#             print(hit.name,len(queryhmmhits))
#             if hit.name not in targets:
#                 targets[hit.name] = []
#             targets[hit.name].append(hit)
#     return targets


# def hitcoverage(
#     hmm_query_length:int,  
#     hit:pyhmmer.plan7.Hit
#     ):
#     """
#         per query (hmm in case of hmmsearch) cumulative coverage.
#     """    
#     poscovered = []
#     for domain in hit.domains:        
#         poscovered += list(
#             range(
#                 domain.alignment.hmm_from,
#                 domain.alignment.hmm_to,
#             )
#         )

#     cov = len(list(set(poscovered)))/hmm_query_length
#     return cov



# def targethitpos(hit , evalue = 1e-10):    
#     fstart = []
#     fend = []             
#     print(len(hit.domains))
#     for d in hit.domains:        
#         if d.i_evalue < evalue:
#             print("here")
#             fstart.append(d.alignment.target_from)
#             fend.append(d.alignment.target_to)    
#     print(fstart,fend)
#     return min(fstart) , max(fend)
    

# def writefasta(sequences:dict,file:str):
#     f = os.path.abspath(file)
#     if not os.path.isdir(os.path.dirname(f)):
#         os.makedirs(os.path.dirname(f) , exist_ok = True)
#     with open(f,'w') as stream:
#         for seqid,seq in sequences.items():
#             assert isinstance(seq,pyhmmer.easel.DigitalSequence)
#             stream.write(">{} {}\n".format(seqid , seq.description.decode('ascii') ))
#             stream.write(seq.textize().sequence + "\n")

# def filterfasta(sequences:list,identifier:list):
#     keep = [] 
#     for seq in sequences:
#         if seq.name.decode('ascii') in identifier:
#             keep.append(seq)
#     return keep


# def maketable(data:list):
#     d = [i.to_dict() for i in data]
#     print(d)
#     df = pd.DataFrame(d)
#     return df

        
# def summarize(nter , cter):
#     import re
#     if re.search("Gly1,Gly2,Gly3",cter):
#         if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
#             flag = "Calcyanin with known N-ter"
#         else:
#             flag = "Calcyanin with new N-ter"
#     elif re.search("Gly1,Gly3",cter) and nter == "Y-type":
#         flag = "Calcyanin with known N-ter"
#     elif re.search ("Gly(1|2|3)",cter):
#         if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
#             flag = "Atypical Gly region with known N-ter"
#         else:
#             flag = "Atypical Gly region with new N-ter"
#     else:
#         flag="Ancestral gly containing protein"
#     return flag

# def domain_to_fasta(hits,records):
#     dom_records = []
#     for h in hits:
#         print(h.domid,h.seqid)
#         dom_records.append(h.extract_from_fasta(records))
#     return dom_records