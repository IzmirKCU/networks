# encoding: utf-8
'''
fairpype.virusPipe --

A script to generate a set of blast databases using the blastdb command
this requires the NCBI Blast+ Toolkit to be installed
available at https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html

@author:     Xuan Jiang/Simon Rayner

@copyright:  2024 Oslo University Hospital. All rights reserved.

@license:    license

@contact:    xuan.jiang@medisin.uio.no/simon.rayner@medisin.uio.no
@deffield    updated: Updated
'''

import sys
import os
from datetime import datetime
import hashlib
import logging
import yaml

import argparse
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import glob


__all__ = []
__version__ = 0.1
__date__ = '2024-09-24'
__updated__ = '2024-09-24'

DEBUG = 1
TESTRUN = 0
PROFILE = 0





class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


def initLogger(projectFile):
    ''' setup log file based on project name'''
    md5string = hashlib.md5(projectFile.encode('utf-8')).hexdigest()
    projectBaseName = os.path.splitext(os.path.basename(projectFile))[0]
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string + ".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.INFO)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)  # set the new handler
    log.addHandler(handler)
    logging.info("+" + "*" * 78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*" * 78 + "+")
    logging.debug("debug mode is on")


def parseArgs(argv):
    '''parse out Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2025 Oslo University Hospital. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % (program_shortdesc, str(__date__))
    projectFile = ""
    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-p", "--projectfile", dest="projectfile", action="store",
                            help="project file in JSON format [default: %(default)s]")

        # Process arguments
        args = parser.parse_args()


        projectFile = args.projectfile
        if not projectFile:
            raise argparse.ArgumentTypeError(" you have to specify a project file")

        print(projectFile)
        return projectFile

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        if DEBUG or TESTRUN:
            raise (e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

def loadConfigData(projectFile):
    with open(projectFile, 'r') as file:

        projectData = yaml.safe_load(file)
    return projectData

import os, glob


def main(argv=None):  # IGNORE:C0111

    if argv is None:
        argv = sys.argv



    projectFile = parseArgs(argv)
    initLogger(projectFile)

    logging.info("project file is <" + projectFile + ">")
    configData = loadConfigData(projectFile)
    logging.info(configData)
    logging.info("done")




if __name__ == "__main__":
    if DEBUG:
        pass
        #sys.argv.append("-h")
        #sys.argv.append("-v")

    if TESTRUN:
        import doctest
        doctest.testmod()

    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'fairpype.virusPipe_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
