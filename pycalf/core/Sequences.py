import os
import io
import sys
import logging
import tempfile
import gzip
import subprocess
import multiprocessing

from shutil import which
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

import numpy as np
import pandas as pd
import pyhmmer.easel
from   Bio import SeqIO
from   Bio.SeqRecord import SeqRecord
from   Bio.SeqFeature import SeqFeature, FeatureLocation

#from pycalf.utils import log

class SequencesIO:
    def _openfile(self,filename):
        if filename.endswith('.gz'):
            return gzip.open(filename,"rt") 
        else:
            return open(filename, "r" )

    def _parse_fasta(self,fasta,src=""):        
        if isinstance(fasta,str):
            if not os.path.exists(fasta):
                logging.error("{} file not found".format(fasta))
                exit(-1)
            else:
                fhl = self._openfile(fasta)
        else:
            fhl = fasta
        
        assert isinstance(fhl,io.TextIOWrapper)

        sequences = []
        for record in SeqIO.parse(fhl,format="fasta"):
            if len(record.seq) < 10000:
                sequences.append(Seq(record,src=src))
            else:
                logging.warning("Sequence {} seems very long and can't \
                    be scan through hmmer".format(record.id))
        return sequences

    def to_fasta(self, filehandle, sequences):
        assert not isinstance(filehandle,str), "Expect file like object not str"
        SeqIO.write(sequences, filehandle,format="fasta")


    def load_msa(self,msafile):
        with pyhmmer.easel.MSAFile("data/msa/LuxC.sto") as msa_file:
            msa_file.set_digital(pyhmmer.easel.Alphabet.amino())
            msa = next(msa_file)

class Hit:
    def __init__(self,        
        hid="",
        target_len:int=0,
        start_location:int=0,
        stop_location:int=0,
        score:float=0,
        identity:float=0,
        coverage:float=0,
        method:str="",
        src:str="",
        ):
        
        
        self.hid = hid
        self.target_len=target_len
        self.start_location = start_location
        self.stop_location = stop_location
        self.score = score
        self.identity=identity
        self.coverage = round(coverage,2) if coverage else coverage
        self.method = method
        self.src = src


class Seq(SeqRecord):     
    def __init__(self,seqrecord,src=None):        
        SeqRecord.__init__(self,seq=seqrecord.seq)
        self.id = seqrecord.id
        self.name= seqrecord.name
        self.seq = self.sanitize_record(seqrecord.seq)
        self.description = seqrecord.description
        self.src = src       
        self.hits = []
        self.res = {}

    def sanitize_record(self,seq):
        return seq.replace("-","")

    def digitize(self,alphabet=pyhmmer.easel.Alphabet.amino()):
        return pyhmmer.easel.TextSequence(
            sequence=str(self.seq),
            name=self.id.encode(),
            description=self.description.encode()
            ).digitize(alphabet)

    def get_feature(self,feature_id=None):
        fseq = []
        for f in self.features:
            keep = True            
            if feature_id:                
                keep = feature_id == f.id
            if keep:
                fseq.append(
                    Seq(SeqRecord(seq=f.extract(self.seq),id=self.id,description=f.id))
                )
        return fseq
        
    def per_residue_annotation(self):
        self.res = {aa:None for aa in range(1,len(self.seq)+1 )} 
        #hits = sorted(self.hits, key=lambda x: x.score, reverse=True)
        self.hits.sort(key=lambda x: x.score, reverse=False)

        for hit in self.hits:
            for aa in range(hit.start_location,hit.stop_location):
                if not self.res[aa]:
                    self.res[aa] = []#                
                self.res[aa].append(hit)
                #self.res[aa].evalue.append(hit.score)   
    
    def _keep_best_annotation(self):
        b = []
        for _ , residue in self.res.items():
            if residue:
                b.append(residue.keep_best())
            else:
                b.append("-")
        return b 

    def addhit(self,hit):
        assert isinstance(hit,Hit)
        self.hits.append(hit)        


class Sequences(SequencesIO):
    def __init__(self,fasta=None,src="",**kwargs):        
        self._sequences = self._parse_fasta(fasta,src=src) if fasta else []

    @property
    def sequences(self):
        return self._sequences

    @sequences.setter
    def sequences(self, value:dict):
        if isinstance(value,list):     
            for seq in value:
                if not isinstance(seq,Seq):
                    raise TypeError("sequence is not a Seq object")
            self._sequences = value
        else:
            raise TypeError("value must be a list of Seq object")        

    def get_seq_by_id(self,sid:str):
        for seq in self._sequences:
            if seq.id == sid:
                return seq

    def get_feature(self,feature_id=None):
        features = []
        for seq in self.sequences:
            features+=seq.get_feature(feature_id=feature_id)
        return features

    def digitize(self,alphabet=pyhmmer.easel.Alphabet.amino()):
        digital_records = []
        for seq in self._sequences:
            digital_records.append(seq.digitize(alphabet))
        return digital_records

    def hamming_distance(self,str1, str2):
        assert len(str1) == len(str2)
        return sum(chr1 != chr2 for chr1, chr2 in zip(str1, str2))/len(str1)     

        
    def hmmsearch(self,hmms:list,cpus=multiprocessing.cpu_count()-1, **kwargs):
        assert isinstance(hmms,list)
        for hmm in hmms:
            print(type(hmm.hmm))
            assert isinstance(hmm.hmm, pyhmmer.plan7.HMM)
        hmms = [hmm.hmm for hmm in hmms]
        hits = pyhmmer.hmmsearch(hmms,self.digitize(),cpus=cpus)

        hmm_datas = [(h.name.decode("UTF-8"),h.M)  for h in hmms]
        #[print(h.name.decode("UTF-8"), h.nseq , h.M ) for h in hmms ]
        seq_with_hit = []
        for hmm_datas , hit_by_hmm in zip(hmm_datas,hits):
            hmm_id,hmm_len=hmm_datas
            for hit in hit_by_hmm:
                sequence_identifier = hit.name.decode("UTF-8")
                if sequence_identifier not in seq_with_hit:
                    seq_with_hit.append(sequence_identifier)
                seq = self.get_seq_by_id(sequence_identifier)
                for dom in hit.domains:
                    
                    dh = Hit(
                            hid=hmm_id,
                            target_len=hmm_len,
                            start_location=dom.alignment.target_from,
                            stop_location=dom.alignment.target_to,
                            score = dom.i_evalue, # np.format_float_scientific( val ,precision=4) ,
                            method="hmmsearch",
                            coverage=(dom.alignment.target_to-dom.alignment.target_from)/hmm_len,#len(seq.seq),
                            src = hmm_id,
                            identity = round(
                                self.hamming_distance(dom.alignment.target_sequence.lower(),dom.alignment.hmm_sequence.lower()),2)                            
                        )
                    seq.addhit(dh)
        
        sequences_with_hit = Sequences()
        sequences_with_hit.sequences = [self.get_seq_by_id(sid) for sid in seq_with_hit]
        return sequences_with_hit


    def blastp(self,blast_subject_fasta,evalue="1e-4",outfmt="10 std slen"):        
        blastpexec = which('blastp')
        if blastpexec is None:
            logging.critical("blastp command not found...")
            raise OSError("blastp command not found...")
        query = tempfile.NamedTemporaryFile(mode="w+")
        self.to_fasta(query,self.sequences)
        query.flush()
        if blastpexec:
            logging.info(blastpexec)
            command = [
                    blastpexec,
                    "-query" , query.name,
                    "-subject", blast_subject_fasta,
                    "-evalue" , str(evalue),
                    "-outfmt" , '"10 std slen"'                    
                ]

        logging.info("running : " + " ".join(command))       
        # using shell=true is not a problem here
        o = subprocess.run(" ".join(command) , capture_output=True , shell=True)         
        res = o.stdout.decode('ascii').strip()

        if res:
            df = pd.read_table( StringIO(res)  , sep=","  , header=None )          
            df.columns = "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore slen".split(" ")
            df.sort_values("pident",inplace=True,ascending=False)    
            df.set_index("qacc",inplace=True)
            df["coverage"] = df.apply(lambda x: (x.send-x.sstart) / x.slen * 100, axis=1 )
            for _ , blastp_hit in df.iterrows():
                dh = Hit(
                        hid=blastp_hit.sacc,
                        target_len=blastp_hit.slen,
                        start_location=blastp_hit.sstart,
                        stop_location=blastp_hit.send,
                        score = blastp_hit.evalue,
                        method="blastp",
                        coverage=(blastp_hit.send-blastp_hit.sstart)/blastp_hit.slen ,
                        src = blastp_hit.sacc,
                        identity = blastp_hit.pident
                    )
                self.get_seq_by_id(_).addhit(dh)

    def to_feature_table(self,add_sequence=True):
        features = []
        for seq in self.sequences:
            for feature in seq.features:
                f = [
                    seq.id,
                    seq.src,
                    feature.type,
                    feature.location.start,
                    feature.location.end,
                    feature.id,
                    feature.qualifiers["identity"],
                    feature.qualifiers["score"],
                    feature.qualifiers["src"],
                    feature.qualifiers["src_len"]
                ]
                if add_sequence:                    
                    f.append(
                        str(feature.extract(seq.seq))
                    )
                features.append(f)
        return pd.DataFrame(features,columns=[
            "sequence_id","sequence_src","feature_type","feature_start","feature_end",
            "feature_id","pident","e-value","feature_src","feature_target_len","feature_seq"
        ]).set_index("sequence_id")
