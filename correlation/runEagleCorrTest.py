from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import argparse
import numpy as np
from applyCorr import (
  corr_simu_test_suite,
  corr_simu_func_suite,
)


def make_parser():
  parser = argparse.ArgumentParser(
    description='run Correlation test on Eagle data'
  )
  parser.add_argument('dataDir', type=str, help='path to input data')
  parser.add_argument('outDir', type=str, help='path to output data')
  parser.add_argument('--norm', action='store_true')
  parser.add_argument('--return_func', action='store_true')
  args = parser.parse_args()
  return args


def run_corr_test(dataDir, outDir, norm=False):
  all_cdm_corr = np.zeros((4, 64))
  all_wdm_corr = np.zeros((4, 64))

  for i in np.arange(1, 5)[::-1]:
    print('Operating on set %d' % i)
    cdm_corr, wdm_corr = corr_simu_test_suite(
      os.path.join(dataDir, 'cdm_diags_%d.rds' % i),
      os.path.join(dataDir, 'wdm_diags_%d.rds' % i),
      normalize=norm,
    )
    all_cdm_corr[i-1, :i**3] = cdm_corr
    all_wdm_corr[i-1, :i**3] = wdm_corr

  np.save(
    os.path.join(outDir, 'cdm_corr_norm(%d).npy' % norm),
    all_cdm_corr,
  )
  np.save(
    os.path.join(outDir, 'wdm_corr_norm(%d).npy' % norm),
    all_wdm_corr,
  )


def run_corr_func_test(dataDir, outDir, norm=False):
  all_cdm_func = np.zeros((4, 64, 100))
  all_wdm_func = np.zeros((4, 64, 100))

  for i in np.arange(1, 5)[::-1]:
    print('Operating on set %d' % i)
    cdm_func, wdm_func = corr_simu_func_suite(
      os.path.join(dataDir, 'cdm_diags_%d.rds' % i),
      os.path.join(dataDir, 'wdm_diags_%d.rds' % i),
      normalize=norm,
    )
    all_cdm_func[i-1, :i**3, :] = cdm_func
    all_wdm_func[i-1, :i**3, :] = wdm_func

  np.save(
    os.path.join(outDir, 'cdm_corr_func_norm(%d).npy' % norm),
    all_cdm_func,
  )
  np.save(
    os.path.join(outDir, 'wdm_corr_func_norm(%d).npy' % norm),
    all_wdm_func,
  )


if __name__ == "__main__":
  args = make_parser()
  if args.return_func:
    run_corr_func_test(args.dataDir, args.outDir, norm=args.norm)
  else:
    run_corr_test(args.dataDir, args.outDir, norm=args.norm)

