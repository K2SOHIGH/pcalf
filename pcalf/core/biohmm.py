import os
import tqdm
import pyhmmer
import logging

#logger = logging.getLogger()

class HmmIO:
    """Class to handle IO operation for Hmm object"""
    
    def load_msa(self,msafile:str,msa_name:str,alphabet=pyhmmer.easel.Alphabet.amino(),format=None):
        """Load a fasta file as pyhmmmer.easel.MSA object.

        Args:
            msafile (str): Path to fasta file
            msa_name (str): Name of your MSA
            alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use
            format (str): pyhmmer.easel.MSAFile format, default is None for auto detection. \
                See pyhmmer.easel.MSAFile documentation for a complete list of format.

        Returns: 
        pyhmmmer.easel.MSA

        Raises:
        ValueError: When format is not a valid MSA format.
        FileNotFound: If msafile doesn't exist.
        """
        if os.path.exists(msafile):
            with pyhmmer.easel.MSAFile(msafile,digital=True,alphabet=alphabet) as msa_file:                
                msa = next(msa_file)
            msa.name=msa_name.encode("UTF-8")
            return msa
        else:
            raise FileNotFoundError("{} doesn't exist".format(msafile))


    def dump(self,filename):
        msa = self._msa
        with open(filename,'w') as fh:
            if isinstance(msa ,pyhmmer.easel.DigitalMSA):
                msa = msa.textize()
            for s in msa.sequences:                
                r = ">{}\n{}\n".format(s.name.decode(),s.sequence)
                fh.write(r)


            


class Hmm(HmmIO):
    """Class to handle both a MSA and its corresponding HMM profile.
    
    Attributes:
        msa (pyhmmer.easel.MSA): Multiple Sequence Alignment.
        hmm (pyhmmer.plan7.HMM): HMM profile.
 
    """

    def __init__(self,msa_name,msafile=None, **kwargs) :        
        """Initializes the instance.

        Args:
            msa_name (str): Name of the msa/hmm object - required
            msafile (str): Path to the alignment file - optional
        """
        self._msa = self.load_msa(msafile,msa_name) if msafile else None        
        self._hmm = self.hmmbuild(self.msa, **kwargs) if msafile else None
        self._msa_file = msafile        
    
    @property
    def msa(self):
        return self._msa

    @property
    def msa_file(self):
        return self._msa_file

    @property
    def hmm(self):
        return self._hmm

    @msa.setter
    def msa(self, value):        
        if isinstance(value,pyhmmer.easel.MSA):                 
            self._msa = value
        else:
            raise TypeError("value must be pyhmmer.easel.MSA")        

    @msa_file.setter
    def msa_file(self, value):        
        self._msa_file = value

    @hmm.setter
    def hmm(self, value):        
        if isinstance(value,pyhmmer.plan7.HMM):                 
            self._hmm = value
        else:
            raise TypeError("value must be pyhmmer.plan7.HMM")        

    def _align(self,
                sequences:list,
                hmm:pyhmmer.plan7.HMM,
                trim:bool=True,
                alphabet:pyhmmer.easel.Alphabet=pyhmmer.easel.Alphabet.amino()):
        """Helper function to align sequences against temporary hmm 

        Args:
            sequences (list): List of  pyhmmer.easel.DigitalSequences - required.
            hmm (pyhmmer.plan7.HMM): Hmm profile.
            trim (bool): N-ter/C-ter trimming, see pyhmmer.hmmalign documentation for details
            alphabet: Kind of alphabet to use.  
        
        Raises:
            TypeError: If one sequence is not an instance of pyhmmer.easel.DigitalSequence.

        Return:
            Tuple of pyhmmer.easel.MSA and pyhmmer.plan7.HMM
        """
        for seq in sequences:
            if not isinstance(seq,pyhmmer.easel.DigitalSequence):
                raise TypeError("TypeError, get {} while pyhmmer.easel.DigitalSequence was expected".format(
                    type(seq)
                ))
        msa = pyhmmer.hmmalign(hmm, 
                               sequences, 
                               digitize=True,
                               trim=trim)     
        msa.name = hmm.name
        hmm = self.hmmbuild( msa )
        return msa,hmm

    def _count_kmers(self,sequence, k_size=5):
        data = {}
        kmers = []
        size = len(sequence)
        for i in range(size - k_size + 1):
            kmer = sequence[i: i + k_size]
            kmers.append(kmer)
            if kmer not in data:
                data[kmer] = 1
                continue        
            data[kmer] += 1 

        return data , kmers 

    def _jaccard_similarity(self, a, b):
        a = set(a)
        b = set(b)
        intersection = len(a.intersection(b))
        union = len(a.union(b))
        return intersection / union


    def hmmalign(self,sequences:list,
                 self_include=True,
                 iterative = False,
                 trim=False,
                 alphabet=pyhmmer.easel.Alphabet.amino()):
        """Align sequences to the hmm 

        Args:
          sequence (list): List of pyhmmer.easel.DigitalSequence - required
          self_include (bool): If set, sequences already present within the hmm will be added to the new hmm produced
          trim (bool): N-ter/C-ter trimming, see pyhmmer.hmmalign documentation for details
          alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.

        Raises:
            TypeError: If one sequence is not an instance of pyhmmer.easel.DigitalSequence.

        Return:
            Hmm -> A new Hmm object including the hmm profile and the multiple sequence alignment.
        """    

        accs = {} # keep sequence name in a list 
        for seq in sequences:            
            if isinstance(seq,pyhmmer.easel.DigitalSequence): 
                accs[seq.name] = seq.textize().sequence
                continue
            raise TypeError("TypeError, get {} while pyhmmer.easel.DigitalSequence was expected".format(
                type(seq)
            ))            
                
        if self_include: # Current MSA is kept.
            for seq in self.msa.sequences:
                if seq.name not in accs:
                    sequences.append(seq)
                else: # Duplicates      
                    k_size = 5              
                    _ , kmer_A = self._count_kmers(seq.textize().sequence,k_size= k_size)
                    _ , kmer_B = self._count_kmers(accs[seq.name],k_size= k_size)
                    jaccard = self._jaccard_similarity(kmer_A,kmer_B)
                    logging.debug('{} already in HMM. Skip.'.format(seq.name.decode()))
                    logging.debug("""Duplicate identifier : {}\nHMM seqlength: {};\nNEW seqlength: {};\n Jaccard similarity: {} [kmer length: {}]""".format(
                        seq.name.decode(),
                        len(accs[seq.name]),
                        len(seq.textize().sequence),
                        jaccard,
                        k_size
                    ))
        
        hmm = self.hmm
        msa = self.msa
        no_seqs = len(sequences)
        logging.debug("Current No nodes within hmm : {} ".format(hmm.nseq) )
        logging.debug("# Sequences to align : {} ".format(no_seqs) )
        query = []
        while sequences:        
            if not iterative:
                query = [sequences.pop() for _ in range(no_seqs)]
            else:                
                query.append(sequences.pop())
                                     
            msa, hmm = self._align(
                        query,
                        hmm,
                        trim=trim,
                        alphabet=alphabet)
            
            assert isinstance(hmm,pyhmmer.plan7.HMM)
            assert isinstance(msa,pyhmmer.easel.MSA)
        
        new_hmm = Hmm(hmm.name)
        new_hmm.hmm = hmm
        new_hmm.msa = msa     
        logging.debug("No of sequence in updated hmm/msa : {}/{}".format(
            hmm.nseq,len(msa.sequences)))
        return new_hmm
        

    def hmmbuild(self,
                msa,
                alphabet=pyhmmer.easel.Alphabet.amino(),
                matrix = "BLOSUM90",
                **kwargs):    
        """Build HMM from MSA 

        Args:
            msa (pyhmmer.easel.MSA): MSA. [required]
            alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.
            matrix (str): Substitution matrix to use.
            **kwargs: Additional parameter for pyhmmer.plan7.builder()
          
        Raises:
            None

        Return:
            pyhmmer.plan7.HMM object
        """

        builder = pyhmmer.plan7.Builder(alphabet,**kwargs)
        builder.score_matrix = matrix
        background = pyhmmer.plan7.Background(alphabet)        
        hmm, _ , __ = builder.build_msa(msa, background)
        return hmm
