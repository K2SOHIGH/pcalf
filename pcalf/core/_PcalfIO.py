import os
import pickle

class HarleyIOError(Exception):
    pass
    

class HarleyIO(HarleyIOError):
    def is_file(self,f):
        if os.path.exists(f):
            return os.path.abspath(f)
        else:
            raise FileNotFoundError

    def is_dir(self,f):
        if os.path.isdir(f):
            return os.path.abspath(f)
        else:
            raise IsADirectoryError
                    

    def is_valid_msa(self,f):
        seqlen = None
        seq = None
        current_len = 0

        with open(f,"r") as stream:
            line = stream.readline()
            while line:
                if line.startswith(">"):
                    if seqlen is None:
                        if current_len > 0:
                            seqlen = current_len
                            current_len=0
                    else:
                        if seqlen != current_len:                            
                            raise HarleyIOError("Sequence in a multiple sequence alignment should have the same length.")
                        current_len=0
                else:                
                    current_len+=len(line.strip())     
                line = stream.readline()
        return True


    def dump(self,file):
        dirname = os.path.abspath(os.path.dirname(file))
        os.makedirs(dirname,exist_ok=True)
        with open(file, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls,file):
        if os.path.exists(file):
            with open(file,'rb') as f:
                return pickle.load(f)
    