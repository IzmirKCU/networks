import sys
from Levenshtein import distance
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import numpy as np
import pandas as pd
import logging
from datetime import datetime

from pathlib import Path
import sys

import os

__all__ = []
__version__ = 0.1
__date__ = '2024-04-03'
__updated__ = '2024-04-03'
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
        
        parser = ArgumentParser(description="a program to calculate Fibonacci's number", formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--miraw_result_file", dest="mirawres", action="store", help="result file from running miraw% [default: %(default)s]")
        parser.add_argument("-m", "--max_dist", dest="maxdist", action="store", help="filter Levenshtein distances greater than this [default: %(default)s]")
    
        # Process arguments
        args = parser.parse_args()

        
        miraw_result_file = args.mirawres
        maxDist = args.maxdist

        # check the user specified a fasta file, if not warn and and exit
        if miraw_result_file:
            print("miraw result file is <" + miraw_result_file + ">")
        else:
            print("you must specify a miraw result file")
            exit
            
        if maxDist:                    
            print("max Levenshtein distance to keep is <" + str(maxDist) + ">")
        
    except KeyboardInterrupt:
        # handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        


def readMiRAWResults(miraw_result_file):
    '''
    load specified fasta file and store header and sequence as entries in two lists
    :param self:
    :return:
    '''

    logging.info("load sequences from miRAW results file <" + miraw_result_file + ">")

    global df_miraw

    # load the fasta lines into a list
    try:
        df_miraw = pd.read_csv(miraw_result_file, delimiter="\t", comment="#",  index_col=False)
        logging.info("read <" + str(len(df_miraw)) + "> lines")
    except Exception as e:
        raise(e)



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


def createDistMatrix():
    uniqueGenes = df_miraw["GeneName"].unique()
    uniqueMiRNAs = df_miraw["miRNA"].unique()
    dim = 0
    if len(uniqueGenes) > len(uniqueMiRNAs):
        dim = len(uniqueGenes)
    else:
        dim = len(uniqueMiRNAs)

    ldDistMatrix = np.zeros((dim, dim))

    for ind in df_miraw.index:
        qGene = df_miraw['GeneName'][ind]
        qMiR = df_miraw['miRNA'][ind]
        im = np.where(uniqueMiRNAs == qMiR)[0][0]
        ig = np.where(uniqueGenes == qGene)[0][0]
        ldDistMatrix[ig][im] = -df_miraw['MFE'][ind]

    return ldDistMatrix

def testPlot(distMatrix):
    
    import networkx as nx
    import pandas as pd
    import matplotlib.pyplot as plt
    #G = nx.from_pandas_edgelist(df_miraw, "miRNA", "GeneName", ["MFE"])
    G=nx.from_numpy_array(distMatrix)
    weights = pd.DataFrame(distMatrix).reset_index().melt('index').values.tolist()
    tuples = [tuple(w) for w in weights]
    filteredTuples = []
    for t in tuples:
        if t[2]!=0.0:
            filteredTuples.append(t)
    
    
    #G.add_weighted_edges_from(tuples)
    pos=nx.spring_layout(G)

    edge_weight = nx.get_edge_attributes(G,'weight')
    d = nx.degree(G)
    node_adjacencies = []
    node_text = []
    for adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))


    nx.draw(G, pos, with_labels=True, font_weight='bold', node_size=[v * 100 for v in dict(d).values()], node_color=node_adjacencies)
    
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_weight)

    plt.show()
    
    print("done")


def main(argv=None): 

    if argv is None:
        argv = sys.argv
        
    # parse_args to get filename
    parseArgs(argv)
    
    # load fasta file
    readMiRAWResults(miraw_result_file)
    distMatrix = createDistMatrix()

    # make a pretty plot
    testPlot(distMatrix)
    
if __name__ == '__main__':

    sys.exit(main())