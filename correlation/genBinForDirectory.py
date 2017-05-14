from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import argparse
import numpy as np
from translate import np_to_bin


def make_parser():
  parser = argparse.ArgumentParser(
    description='generate binaries from directory'
  )
  parser.add_argument('dataDir', type=str, help='path to input data')
  parser.add_argument('outDir', type=str, help='path to output data')
  args = parser.parse_args()
  return args


args = make_parser()
dataDir = os.listdir(args.dataDir)
for dataFile in dataDir:
  dataFileArr = dataFile.split('.')
  if dataFileArr[-1] == 'npy':
    dataFileArr[-1] = 'bin'
    newDataFile = '.'.join(dataFileArr)
    np_to_bin(
      os.path.join(args.dataDir, dataFile), 
      outputfile=os.path.join(args.outDir, newDataFile),
    )
