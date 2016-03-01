import logging
import colorlog
from logging.config import fileConfig


class mylogger():
    def __init__(self, sdict, logfn):
        fileConfig('./logging_config.ini', defaults={'logfilename': logfn})
        self.logger = logging.getLogger()
        self.sdict = sdict

    def __call__(self, type, key, value, msg):
        if key:
            self.sdict[key] = value
        if type==1:
            self.logger.info(msg)
        elif type==2:
            self.logger.warning(msg)
        elif type==3:
            self.logger.error(msg)
        elif type==4:
            self.logger.exception(msg)

