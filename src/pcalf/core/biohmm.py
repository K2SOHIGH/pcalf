import os
import tqdm
import pyhmmer


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
        self._hmm = self.hmmbuild(self.msa, **kwargs) if msafile else None #self.load_hmm(hmmfile,**kwargs)
    
    @property
    def msa(self):
        return self._msa

    @property
    def hmm(self):
        return self._hmm

    @msa.setter
    def msa(self, value):        
        if isinstance(value,pyhmmer.easel.MSA):                 
            self._msa = value
        else:
            raise TypeError("value must be pyhmmer.easel.MSA")        

    @hmm.setter
    def hmm(self, value):        
        if isinstance(value,pyhmmer.plan7.HMM):                 
            self._hmm = value
        else:
            raise TypeError("value must be pyhmmer.plan7.HMM")        


    def hmmalign(self,sequences:list,self_include=True,trim=True,alphabet=pyhmmer.easel.Alphabet.amino()):
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
        if self_include:
            sequences += self.msa.sequences

        for seq in sequences:            
            if not isinstance(seq,pyhmmer.easel.DigitalSequence): 
                raise TypeError("TypeError, get {} while pyhmmer.easel.DigitalSequence was expected".format(
                    type(seq)
                ))
        
        msa = pyhmmer.hmmalign( self.hmm, sequences, digitize=True , trim=trim )         
        if not msa.name:
            msa.name=self.hmm.name 
            
        hmm = self.hmmbuild(msa)
        
        assert isinstance(hmm,pyhmmer.plan7.HMM)
        assert isinstance(msa,pyhmmer.easel.MSA)
        
        new_hmm = Hmm(hmm.name)
        new_hmm.hmm = hmm
        new_hmm.msa = msa

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
        builder.score_matrix = "BLOSUM90"
        background = pyhmmer.plan7.Background(alphabet)        
        hmm, _ , __ = builder.build_msa(msa, background)
        return hmm

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
        msa = pyhmmer.hmmalign(hmm, sequences , digitize=True ,trim=trim)     
        msa.name = hmm.name
        hmm = self.hmmbuild( msa )
        return msa,hmm

    def iterative_align(self,
                sequences:list,
                self_include=True,
                trim=True,
                alphabet=pyhmmer.easel.Alphabet.amino()):                        

        """Align sequence to HMM iteratively

        Args:
            sequences (list): List of pyhmmer.easel.DigitalSequences - required.
            self_include (bool): If set, sequences already present within the hmm will be added to the new hmm produced
            trim (bool): N-ter/C-ter trimming, see pyhmmer.hmmalign documentation for details
            alphabet (pyhmmer.easel.Alphabet): Kind of alphabet to use.  
        
        Raises:
            TypeError: If one sequence is not an instance of pyhmmer.easel.DigitalSequence.

        Return:
            Tuple of pyhmmer.easel.MSA and pyhmmer.plan7.HMM
        """

                
        
        hmm = self.hmm
        msa = None
        if self_include:
            msa = self.msa
        #seq_iterator = iter(sequences)
        #while (seq := next(seq_iterator, None)) is not None:                                                
        for seq in tqdm.tqdm(sequences,colour="yellow"):
            if not isinstance(seq,pyhmmer.easel.DigitalSequence):
                raise TypeError("Wrong type - {} , expected pyhmmer.easel.DigitalSequence".format(type(seq)))            
            
            lseq = [seq] 
            if msa:
                lseq += msa.sequences
            
            
            msa, hmm = self._align(
                        lseq,
                        hmm,
                        trim=trim,
                        alphabet=alphabet)
                        
        assert isinstance(hmm,pyhmmer.plan7.HMM)
        assert isinstance(msa,pyhmmer.easel.MSA)
        
        new_hmm = Hmm(hmm.name)
        new_hmm.hmm = hmm
        new_hmm.msa = msa

        return new_hmm