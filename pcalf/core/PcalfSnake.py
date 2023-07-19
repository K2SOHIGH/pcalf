import os
import logging
import re
import yaml
import multiprocessing


logger = logging.getLogger("pcalf-workflow")


class Snakemake:
    def __init__(self,workflow):
        self.workflow = workflow
        self.configfile = self.get_configfile()
        self.config = self.load_default_config()
        # self.default_snakargs = ["--rerun-triggers", "mtime", "--use-conda" ,"--jobs" , ]

    def get_snakefile(self):            
        return os.path.join(self.workflow,"Snakefile")
        
    def get_configfile(self):        
        return os.path.join(self.workflow,"config","config.yaml")

    def load_default_config(self):
        return yaml.load(open(self.get_configfile()),Loader=yaml.SafeLoader)
            
    def dump_config(self,file):
        dirname = os.path.abspath(os.path.dirname(file))
        os.makedirs(dirname,exist_ok=True)
        yaml.dump(self.config, open(file,'w') )
        return file

    def run(self, snakargs = []):
        args =  ['--snakefile' , self.get_snakefile()]
        if "configfile" not in snakargs:
            args += ["--configfile", self.configfile]
        if "--use-conda" not in snakargs:
            args += ["--use-conda"]
        if "--jobs" not in snakargs and not re.search("-j[0-9]+"," ".join(snakargs)):
            args += ["--jobs",str(multiprocessing.cpu_count())]

        args += snakargs        
    
        logger.info("""running : 
            snakemake %s """ % " ".join(args))
        
        # snakex = snakemake.main(args)
        # Running snakemake through main API seems to exit totaly 
        # when the workflow is completed (on success and on error).
        # Therefore , snakemake is run from subprocess.
        # (+ deploying the workflow on a cluster is easier this way)
        cmd = "snakemake {}".format( " ".join(args))
        snakex = os.system(cmd)
        if snakex != 0:
            logger.error("Hum ... something went wrong while executing the module ... :( ")
            exit(-1)            
        logger.info("Great , workflow finished without error :)" )       
        if "-n" in args or "--dryrun" in args or "--dry-run" in args:
            logger.warning("ps : it was a dryrun !")
            exit(0)
        return snakex    
