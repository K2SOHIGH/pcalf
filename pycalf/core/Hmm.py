import pyhmmer

class HmmIO:
    def __init__():
        pass
    
    def load_msa(self,msafile:str,msa_name:str,digital:bool=True,alphabet=pyhmmer.easel.Alphabet.amino()):
        assert isinstance(msa_name,str)
        with pyhmmer.easel.MSAFile(msafile,digital=digital) as msa_file:
            #msa_file.set_digital(alphabet)
            msa = next(msa_file)
        msa.name=msa_name.encode("UTF-8")
        return msa


class Hmm(HmmIO):
    def __init__(self,msafile=None,msa_name="", **kwargs) :        
        self._msa = self.load_msa(msafile,msa_name) if msafile else None        
        self._hmm = self._hmmbuild(self.msa, **kwargs) if msafile else None #self.load_hmm(hmmfile,**kwargs)

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
            self._msa = value
        else:
            raise TypeError("value must be pyhmmer.plan7.HMM")        


    def hmmalign(self,sequences:list,self_include=True,digital=False,trim=True,alphabet=pyhmmer.easel.Alphabet.amino()):                
        if self_include:
            sequences += self.msa.sequences
        digitalseqs = [s.digitize(alphabet) if not isinstance(s,pyhmmer.easel.DigitalSequence) else s for s in sequences   ]
        msa = pyhmmer.hmmalign(self.hmm, digitalseqs , digitize=digital ,trim=trim) 
        
        if not msa.name:
            msa.name=self.hmm.name 
            
        hmm = self._hmmbuild(msa.digitize(
            pyhmmer.easel.Alphabet.amino()
        ))

        return hmm,msa

    def _hmmbuild(self,msa,**kwargs):
        alphabet = pyhmmer.easel.Alphabet.amino()
        builder = pyhmmer.plan7.Builder(alphabet,**kwargs)
        background = pyhmmer.plan7.Background(alphabet)
        builder.score_matrix = "BLOSUM90"
        hmm, _ , __ = builder.build_msa(msa, background)
        return hmm


if __name__ == "__main__":
    import Sequences
    aln1 = "../test/test1.fa.msa" 
    fa2 = Sequences.Sequences("../test/test2.fa.msa")
    hmm = Hmm(aln1)
    print( hmm.msa,
            hmm.hmm.M , hmm.hmm.nseq )
    hmm,msa = hmm.hmmalign(fa2.sequences)
    print(
        hmm.M , hmm.nseq )
