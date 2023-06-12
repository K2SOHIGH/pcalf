import os
import logging
import inspect 


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
    lightpurple = "\x1b[38;5;69m" #[38;5;141m" #"\x1b[38;5;84m"
    pink90 = "\x1b[38;5;162m" #"\x1b[38;5;214m"
    
    format = "PCALF | [%(levelname)s] %(message)s"
    
    FORMATS = {
        logging.DEBUG: green + format + reset,
        logging.INFO: lightpurple + format + reset,
        logging.WARNING: pink90 + format + reset,
        logging.ERROR: bold_red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)                
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

