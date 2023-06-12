import os
import logging

class CustomFormatter(logging.Formatter):    
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    black = "\u001b[30m"    
    green = "\u001b[32m"
    blue = "\u001b[34m"
    magenta = "\u001b[35m"
    cyan = "\u001b[36m"
    white = "\u001b[37m"


    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: reset + format + reset,
        logging.INFO: magenta + format + reset,
        logging.WARNING: cyan + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)



def setlogger(name="__main__"):    
    logger = logging.Logger(name)          
    return logger 

def file_handler(logfile,level):    
    dirname = os.path.abspath(
        os.path.dirname(logfile)
    )
    os.makedirs(dirname,exist_ok=True)
    handler = logging.FileHandler( logfile , mode = 'a' )    
    handler.setFormatter( CustomFormatter() )    
    handler.setLevel(level)
    return handler
    
def stream_handler(level):    
    handler = logging.StreamHandler( )    
    handler.setFormatter( CustomFormatter() )    
    handler.setLevel(level)
    return handler
