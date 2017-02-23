############################################################################
#This script is modified from the original code
#obtained from https://github.com/bioinform/metasv/blob/master/metasv/external_cmd.py
#Copyright (c) 2014, Bina Technologies inc.
############################################################################


import time
import shlex
import subprocess
from threading import Timer
import unittest
import os
from utils import *

class TimedExternalCmd:
    def __init__(self, cmd, logger, raise_exception=False, env_dict={}):
        self.cmd = shlex.split(cmd)
        self.p = None
        self.did_timeout = False
        self.logger = logger
        self.raise_exception = raise_exception
        self.env_dict = env_dict
    def enforce_timeout(self):
        self.p.terminate()
        self.did_timeout = True
    def run(self, cmd_log_fd_out=None, cmd_log_fd_err=None,  cmd_log="", msg="", timeout=None):
        self.logger.info("Task: %s " % (msg))
        self.logger.info("Running \"%s\" " % (" ".join(self.cmd)))
        cmd_log_fd_err = cmd_log_fd_err or cmd_log_fd_out
        if self.env_dict:
            my_env = os.environ.copy()
            for k,v in self.env_dict.iteritems():
                my_env[k] = v
            self.p = subprocess.Popen(self.cmd, stderr=cmd_log_fd_err, stdout=cmd_log_fd_out, env=my_env)
        else:
            self.p = subprocess.Popen(self.cmd, stderr=cmd_log_fd_err, stdout=cmd_log_fd_out)
        
        start_time = time.time()
        if timeout:
            t = Timer(timeout, self.enforce_timeout)
            t.start()
        self.p.wait()
        if timeout:
            t.cancel()
            if self.did_timeout:
                if not self.raise_exception:
                    self.logger.error("Timed out after %d seconds.", timeout)
                    return None
                else:
                    self.logger.error("Aborting!")
                    raise Exception("Timed out after %d seconds."%timeout)
        retcode = self.p.returncode
        if retcode == 0:
            self.logger.info("Done %s " % msg)
        else:
            if self.raise_exception:
                self.logger.info("Returned code %d (%g seconds)" % (retcode, time.time() - start_time))            
                self.logger.error("Aborting!")
                if cmd_log:
                    raise Exception("Failed %s. Log file: %s" % (msg,cmd_log))
                else:
                    raise Exception(msg)
        self.logger.info("Returned code %d (%g seconds)" % (retcode, time.time() - start_time))
        return retcode


class TestTimedExternalCmd(unittest.TestCase):
    def test_run_complete(self):
        cmd = TimedExternalCmd("sleep 1", self.logger)
        self.assertEqual(cmd.run(timeout = 2), 0)
        self.assertFalse(cmd.did_timeout)
        return

    def test_run_timeout(self):
        start_tick = time.time()
        cmd = TimedExternalCmd("sleep 2", self.logger)
        cmd.run(timeout = 1)
        run_time = time.time() - start_tick
        self.assertTrue(cmd.did_timeout)
        self.assertAlmostEqual(run_time, 1, delta=0.2)
        return

    def test_run_no_timeout(self):
        cmd = TimedExternalCmd("sleep 1", self.logger)
        retcode = cmd.run()
        self.assertEqual(cmd.run(), 0)
        self.assertFalse(cmd.did_timeout)
        return

    def test_run_fail(self):
        cmd = TimedExternalCmd("sleep 1 2 3", self.logger)
        retcode = cmd.run(timeout = 1)
        self.assertIsNotNone(retcode)
        self.assertIsNot(retcode, 0)
        return

    logger = None


if __name__ == '__main__':
    TestTimedExternalCmd.logger = logging.getLogger(__name__)
    unittest.main()
