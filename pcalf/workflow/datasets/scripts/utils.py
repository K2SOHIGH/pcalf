import os
import yaml
import logging
import glob

def bins2batches(bins,batch_size):    
    batches = {"batch-{}".format(x//batch_size):list(bins.keys())[x:x+batch_size] for x in range(0,len(bins.keys()),batch_size) if x <= len(bins.keys())}
    return batches

def input_from_dir(input , extension):
    if input:
        if os.path.isdir(input):
            #IDS, = glob_wildcards(input + "/{id}." + extension)
            IDS = glob.glob(input + "/*" + extension)
            return { 
                os.path.basename(i).split(extension)[0] : i for i in IDS
                }

def input_from_yaml(input):  
    if input:
        if os.path.exists(input):
            conff = open(input)
            datas = yaml.load(conff,Loader=yaml.FullLoader)    
            return datas
        else:
            msg="""WORKFLOW INPUT : {} not found
            """
            logging.error(msg.format(input))
            raise FileNotFoundError(msg.format(input))
    else:        
        return None

def input_is_fasta(input):
    if input:
        if os.path.exists(input):
            n = os.path.basename(input)
            return {n : os.path.abspath(input)}


def parse_input(input , extension):
    if input:
        if os.path.isdir(input):
            return input_from_dir(input, extension)
        elif os.path.isfile(input):
            if input.endswith(".yaml"):
                return input_from_yaml(input)
            return input_is_fasta(input)
        else:
            return None
    return None