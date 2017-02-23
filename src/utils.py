import logging
import os
import time

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, filename="run-%s.log"%time.strftime("%Y%m%d-%H%M%S"), filemode="w")

logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

def create_dirs(dirlist):
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating directory %s" % (dirname))
            os.makedirs(dirname)

