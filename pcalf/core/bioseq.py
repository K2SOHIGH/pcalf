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



class SequencesIO:
    """Class to handle IO operation for Sequences object"""
    def _openfile(self,filename):
        """Open file in gzip or text mode.

        Args:
        filename (str): Path to file

        Returns:
        return Handle

        Raises:
        None
        """
        
        if filename.endswith('.gz'):
            return gzip.open(filename,"rt") 
        else:
            return open(filename, "r" )

    def _parse_fasta(self,fasta,src=""):   
        """Parse fasta file and return a list of Seq object.

        Args:
        fasta: Path to  file or file-like object.
        src (str): Name used to keep trace of the original file for each record. Useful when multiple fasta file will be mixed. 

        Returns:
        return a list of Seq.

        Raises:
        FileNotFoundError; If fasta doesn't exist. 
        """     
        if isinstance(fasta,str):
            if not os.path.exists(fasta):
                raise FileNotFoundError("{} not found".format(fasta))  
            else:
                fhl = self._openfile(fasta)
        else:
            fhl = fasta
        
        assert isinstance(fhl,io.TextIOWrapper)

        sequences = []
        ids = []
        #'ACDEFGHIKLMNPQRSTVWY-BJZOUX*~'
        for record in SeqIO.parse(fhl,format="fasta"):            
            if len(record.seq) < 10000:
                #suffix = 1
                if record.id in ids:
                    #record.id = record.id + "_pcalf{}".format(suffix)
                    logging.warning('Duplicate sequence identifier for {}. this one will be ignored...'.format(record.id))
                    continue
                    #suffix+=1
                ids.append(record.id)            
                record.seq = record.seq.strip()
                sequences.append(Seq(record,src=src))
            else:                
                logging.warning("Sequence {} seems very long and can't be scan through hmmer".format(record.id))
        return sequences

    def to_fasta(self , filehandle ):
        """Dump fasta records.

        Args:
        filehandle (File like object): Handle.
        
        Returns:
        None

        Raises:
        FileNotFoundError: If fasta can't be created.
        AssertionError: If filehandle is not a file-like object.
        """    
        assert not isinstance(filehandle,str), "Expect file like object not str"
        SeqIO.write(self.sequences, filehandle,format="fasta")


class Hit:
    """Class to store scores and location of features.
    
    Attributes:
        hid (str): Name of the hit.
        target_len (int): Length of the target.
        start_location (int): Start location.
        stop_location (int): Stop location.
        score (float): E-value.
        identity (float): Percentage of identity.
        coverage (float): Coverage.
        method (str): Method use to produced the hit.
        src (str): Additional information about the origin of the hit.
    """
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
        self.score = float('{:0.3e}'.format(score).strip("'")) if score else score
        self.identity= round(identity,3) if identity else identity
        self.coverage = round(coverage,2) if coverage else coverage
        self.method = method
        self.src = src


class Seq(SeqRecord):  
    """Seq object extend the SeqRecord class from biopython
    
    Attributes:
        hits (list): List of Hit object if any. 
        res (dict): Dictionnary with residue as key.
    """   

    def __init__(self,seqrecord,src=None):      
        SeqRecord.__init__(self,seq=seqrecord.seq , features = seqrecord.features)
        self.id = seqrecord.id
        self.name= seqrecord.name        
        self.description = seqrecord.description
        self.src = src       
        self.hits = []
        self.res = {}
        self.sanitize_record()
        

    def sanitize_record(self):
        """Remove gap character from sequence if any (inplace).
        
        Args:          
            None

        Raises:
            None

        Return:
            None
        """
        if self.seq:
            self.seq = self.seq.replace("-", "")

    def digitize(self,alphabet=pyhmmer.easel.Alphabet.amino()):
        """Convert Seq to pyhmmer.easel.DigitalSequence.

        Args:          
          alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.

        Raises:
            None

        Return:
            pyhmmer.easel.DigitalSequence.
        """
        easelseq = pyhmmer.easel.TextSequence(
            sequence=str(self.seq),
            name=self.id.encode(),
            description=self.description.encode()
            )
        return easelseq.digitize(alphabet)

    def get_feature(self,feature_id=None):
        """Extract features from Seq.

        Args:          
          feature_id (str): A feature identifier, if None all features will be extracted.

        Raises:
            None

        Return:
            list of Seq objects, one per feature.
        """
        fseq = []
        for f in self.features:
            keep = True            
            if feature_id:                
                keep = feature_id == f.id
            if keep:
                fseq.append(
                    Seq(SeqRecord(seq=f.extract(self.seq),id=self.id, description="",name=f.id, features = [f] )    )
                )
        return fseq
        
    def per_residue_annotation(self):
        """Put annotation hits under the residues (self.res) they cover [sorted by score].
        
        Args:
            None
        Raises:
            None
        Return:
            None
        """
        
        self.res = {aa:None for aa in range(1,len(self.seq)+1 )} 
        #hits = sorted(self.hits, key=lambda x: x.score, reverse=True)
        self.hits.sort(key=lambda x: x.score, reverse=False)

        for hit in self.hits:
            for aa in range(hit.start_location,hit.stop_location):
                if not self.res[aa]:
                    self.res[aa] = []#                
                self.res[aa].append(hit)


    def addhit(self,hit):
        """Add Hit object to Seq.
        
        Args:
            hit (Hit): Hit object.
        Raises:
            TypeError: If hit is not an instance of Hit.
        Return:
            None
        """
        if not isinstance(hit,Hit):
            raise TypeError("Except Hit object not {}".format(type(hit)))            
        self.hits.append(hit)        


class Sequences(SequencesIO):
    """Class to handle multiple Seq object.
    
    Attributes:
        sequences (list): List of Seq object if any.         
    """   
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
        """Get a single sequence by its identifier.
        
        Args:
            sid (str): sequence identifier.
        Raises:
            ValueError: If sequence identifier is not found.
        Return:
            Seq object
        """
        for seq in self._sequences:
            if seq.id == sid:
                return seq
        raise ValueError("can't find {}".format(sid))

    def get_feature(self,feature_id=None):
        """Get features from multiple sequences.
        
        Args:
            feature_id (str): A feature identifier.
        Raises:
            None
        Return:
            list of features as Seq objects.
        """        
        features = []
        for seq in self.sequences:
            features+=seq.get_feature(feature_id=feature_id)
        return features

    def digitize(self,alphabet=pyhmmer.easel.Alphabet.amino()):
        """Convert all sequences to pyhmmer.easel.DigitalSequence.
        
        Args:
            alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.
        Raises:
            None
        Return:
            list of pyhmmer.easel.DigitalSequence.
        """
        digital_records = []
        for seq in self._sequences:     
            digital_records.append(seq.digitize(alphabet))
        return digital_records

    def hamming_distance(self,str1, str2):
        """Compute the hamming distance between two string divided by the length of the alignment.
        
        Args:
            str1 (str): First string.
            str2 (str): Second string.
        Raises:
            ValueError: If str1 and str2 do not have the same length.
        Return:
            float
        """
        if len(str1) != len(str2):
            raise ValueError("Str1 and Str2 should have the same length.")
        return sum(chr1 != chr2 for chr1, chr2 in zip(str1, str2))/len(str1)     

    def hmmsearch(self,
                    hmms:list,
                    alphabet=pyhmmer.easel.Alphabet.amino(), 
                    cpus=multiprocessing.cpu_count()-1, 
                    **kwargs):
        """Search sequences against  HMM profile(s).
        
        Args:
            hmms (list): List of pyhmmer.plan7.HMM profile(s).
            alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.
            cpus (int): Number of cpu to use.
        Raises:
            AssertionError: If hmms is not a list or of an HMM is not a .
            TypeError: If at least one profile within hmms is not a pyhmmer.plan7.HMM.
        Return:
            Sequences object -> a new Sequences object containing sequence (Seq) with at least one hit against one or more profiles.
        """                    
        assert isinstance(hmms,list)
        for hmm in hmms:            
            if not isinstance(hmm, pyhmmer.plan7.HMM):
                raise TypeError("Except pyhmmer.plan7.HMM but received {}".format(type(hmm)))

        hits = pyhmmer.hmmsearch(hmms,self.digitize(alphabet),cpus=cpus,**kwargs)
        hmm_datas = [(h.name.decode("UTF-8"),h.M)  for h in hmms]
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

        return self



    def blastp(self,blast_subject_fasta,evalue=1e-4,outfmt="10 std slen"):   
        """Perform a sequence vs sequence blastp search and append \ 
        resulting hits to the hits attribute of their respective sequence.
        
        Args:
            blast_subject_fasta (str): path to another fasta file.
            evalue (float): E-value threshold.
            outfmt (str): blast output format.
        Raises:
            OSError: If blastp command is not found.
            TypeError: If at least one profile within hmms is not a pyhmmer.plan7.HMM.
        Return:
            None
        """  
        blastpexec = which('blastp')
        if blastpexec is None:
            #logging.critical("blastp command not found...")
            raise OSError("blastp command not found...")
        query = tempfile.NamedTemporaryFile(mode="w+")
        self.to_fasta(query)
        query.flush()
        if blastpexec:            
            command = [
                    blastpexec,
                    "-query" , query.name,
                    "-subject", blast_subject_fasta,
                    "-evalue" , str(evalue),
                    "-outfmt" , '"10 std slen"'                    
                ]

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
        return self        

    def to_feature_table(self,add_sequence=True,feature_id=""):
        """Generate a feature table.
        
        Args:
            add_sequence (bool): If set, the feature sequence (str) will be append to the table.
            feature_id (str): If set, only feature with the same identifier will be considered.            
        Raises:
            None
        Return:
            pd.DataFrame
        """          
        features = []
        for seq in self.sequences:         
            for feature in seq.features:
                keep=True
                if feature_id:
                    if feature_id!=feature.id:
                        keep=False
                if keep:
                    f = [
                        seq.id,
                        seq.src,
                        feature.type,
                        feature.location.start,
                        feature.location.end,
                        feature.id,
                        feature.qualifiers["identity"],
                        feature.qualifiers["coverage"],
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
            "feature_id","pident","coverage","e-value","feature_src","feature_target_len","feature_seq"
        ]).set_index("sequence_id")

    def to_hits_table(self):
        """Generate a hits table.
        
        Args:
            self
        Raises:
            None
        Return:
            pd.DataFrame
        """          
        hits_datas = []
        for seq in self.sequences:
            for hit in seq.hits:
                f = [
                    seq.id,
                    seq.src,
                    hit.target_len,                    
                    hit.start_location,
                    hit.stop_location,
                    hit.score,
                    hit.identity,
                    hit.coverage,
                    hit.method,
                    hit.src
                ]
                hits_datas.append(f)
        return pd.DataFrame(hits_datas,columns=[
            "sequence_id","sequence_src","hit_target_len","hit_start","hit_stop",
            "hit_e_value","hit_pident","hit_coverage","hit_method","hit_src"
        ]).set_index("sequence_id")
  