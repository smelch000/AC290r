#!/usr/bin/env python

import sys, os
import logging
from logging import handlers
import tempfile
import io
import errno

from contextlib import contextmanager
# from wurlitzer import pipes, sys_pipes

@contextmanager
def stdout_redirector(stream):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()
    # original_stdout_fd = 1

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        # M.MUPHY.fflush(0)
        # libf90.fflush(0)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)

class Mlogging(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, debug):

        # logging defs starts here
        logging.basicConfig(filename='.MagicLog', 
                        filemode='w', 
                        format='[%(asctime)s] [%(levelname)s] (%(threadName)-10s) %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',)

        self.logger = logging.getLogger('')

        if debug == 'yes':
            self.log_level = logging.DEBUG
        else:
            self.log_level = logging.INFO

        self.linebuf = ''

        self.logger.setLevel(self.log_level)

        #mail_handler = logging.handlers.SMTPHandler(mailhost='127.0.0.1', 
        #                                        fromaddr='Magic@universe.com', 
        #                                        toaddrs='simone.melchionna@gmail.com', 
        #                                        subject='Magic Info')
        #mail_handler.setLevel(logging.INFO)

        # console handler with a higher log level
        consolle_handler = logging.StreamHandler()
        consolle_handler.setLevel(logging.ERROR)

        # self.logger.addHandler(mail_handler)
        self.logger.addHandler(consolle_handler)
        logging.info('Starting...')

        """
        stdout_logger = logging.getLogger('STDOUT')
        #sl = StreamToLogger(stdout_logger, logging.INFO)
        #sys.stdout = sl
        sys.stdout = self

        stderr_logger = logging.getLogger('STDERR')
        #sl = StreamToLogger(stderr_logger, logging.ERROR)
        #sys.stderr = sl
        sys.stderr = self
        """

        # sys.stdout = stdout_redirector
        # sys.stdout.write = self.write
        # sys.stderr = stderr_redirector

        # print("Test to standard out")
        #raise Exception('Test to standard error')
        #sys.exit(1)

    @staticmethod
    def info(buf):
        logging.info(buf)

    @staticmethod
    def debug(buf):
        logging.debug(buf)

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())


