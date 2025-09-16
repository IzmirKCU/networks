import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import hashlib
from datetime import datetime

import logging
import os
from pathlib import Path
import sys

global df_miraw_results
global df_miraw_filtered_results
global probabilityMin, mfeMin
__all__ = []
__version__ = 0.1
__date__ = '2024-04-03'
__updated__ = '2024-04-03'


def initLogger(md5string, gffFile):

    ''' setup log file based on project name'''
    projectBaseName = Path(gffFile).stem

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    log.setLevel(logging.DEBUG)
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler
    log.addHandler(handler)
    logging.info("+" + "*"*78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*"*78 + "+")
    logging.debug("debug mode is on")


def parseArgs(argv):
    '''parse out Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    # program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2024 Oslo University Hospital. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % (program_name, str(__date__))
    try:
        # Setup argument parser
        parser = ArgumentParser(description="wiRAW Result File", formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-m", "--miraw_result_file", dest="mirrawresultfile", action="store",
                            help="gfflist file name [default: %(default)s]")
        parser.add_argument("-p", "--probability_min", dest="probmin", action="store",
                            help="only keep predictions with probability above this cutoff [default: %(default)s]")
        parser.add_argument("-e", "--mfe_min", dest="mfemin", action="store",
                            help="only keep predictions with mean free energy below this cutoff[default: %(default)s]")

        global probabilityMin, mfeMin
        # Process arguments
        args = parser.parse_args()
        mirResultFile = args.mirrawresultfile
        probabilityMin = args.probmin
        mfeMin = args.mfemin

        # check the user specified a fasta file, if not warn and and exit
        if mirResultFile:
            logging.info("miRAW result file is <" + mirResultFile + ">")
        else:
            logging.error("you must specify a miRAW result file")
            exit()

    except KeyboardInterrupt:
        # handle keyboard interrupt #
        return 0
    except Exception as e:
        print(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    return mirResultFile


def loadResultFile(miraw_result_file):
    global df_miraw_results

    try:
        logging.info("loading miRAW result file <" + miraw_result_file + ">")
        df_miraw_results = pd.read_csv(miraw_result_file, delimiter='\t', comment="#")
        logging.info("read <" + str(len(df_miraw_results)) + "> lines")

    except Exception as e:
        print(e)
        indent = len(miraw_result_file) * " "
        sys.stderr.write(miraw_result_file + ": " + repr(e) + "\n")
        return 2


def filterMiRAWResults():
    global df_miraw_results
    global df_miraw_filtered_results
    global mfeMin, probabilityMin
    df_miraw_filtered_results = df_miraw_results[(df_miraw_results['Prediction']>float(probabilityMin)) & (df_miraw_results['MFE']<float(mfeMin))]
    logging.info("filtering removed <" + str(len(df_miraw_results)-len(df_miraw_filtered_results)) + ">")

def summariseMiRAWResults():
    logging.info("miRNA occurrences are as follows:")
    logging.info(df_miraw_filtered_results['miRNA'].value_counts())

    logging.info("Probability distribution is as follows:")
    logging.info(df_miraw_filtered_results['Prediction'].value_counts())

    logging.info(df_miraw_filtered_results['MFE'].value_counts())


def writeFilteredResults(miraw_result_file):
    global df_miraw_filtered_results
    global probabilityMin, mfeMin
    mfeString =""
    if float(mfeMin) < 0:
        mfeString = 'm' + str(float(mfeMin) * 1.0)
    else:
        mfeString = mfeMin
    filteredFilename = Path(miraw_result_file).stem + "__p" + str(probabilityMin)+ "__e" + mfeString + ".tsv"
    filteredFilename = Path(Path(miraw_result_file).parent, filteredFilename)
    df_miraw_filtered_results.to_csv(filteredFilename, sep='\t')

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # parse_args to get filename
    mirResultFile = parseArgs(argv)
    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    initLogger(md5String, mirResultFile)
    loadResultFile(mirResultFile)
    filterMiRAWResults()
    summariseMiRAWResults()
    writeFilteredResults(mirResultFile)
    logging.info("program completed")


if __name__ == '__main__':
    sys.exit(main())
