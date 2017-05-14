from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from applyCorr import (
  corr_voronoi_test_suite,
  corr_voronoi_func_suite
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
  all_base_corr = np.zeros((100, 15))
  all_foam_corr = np.zeros((100, 9, 15))

  for i in range(1, 101):
    print('Operating on set %d' % i)
    base_corr, foam_corr = corr_voronoi_test_suite(
      os.path.join(dataDir, 'baseline%d.rds' % i),
      os.path.join(dataDir, 'foam%d.rds' % i),
      normalize=norm,
    )
    all_base_corr[i-1, :] = base_corr
    all_foam_corr[i-1, :, :] = foam_corr

  np.save(
    os.path.join(outDir, 'base_corr_norm(%d).npy' % norm),
    all_base_corr,
  )
  np.save(
    os.path.join(outDir, 'foam_corr_norm(%d).npy' % norm),
    all_foam_corr,
  )


def run_corr_func_test(dataDir, outDir, norm=False):
  all_base_func = np.zeros((100, 15, 100))
  all_foam_func = np.zeros((100, 9, 15, 100))

  for i in range(1, 101):
    print('Operating on set %d' % i)
    base_func, foam_func = corr_voronoi_func_suite(
      os.path.join(dataDir, 'baseline%d.rds' % i),
      os.path.join(dataDir, 'foam%d.rds' % i),
      normalize=norm,
    )
    all_base_func[i-1, :, :] = base_func
    all_foam_func[i-1, :, :, :] = foam_func

  np.save(
    os.path.join(outDir, 'base_corr_func_norm(%d).npy' % norm),
    all_base_func,
  )
  np.save(
    os.path.join(outDir, 'foam_corr_func_norm(%d).npy' % norm),
    all_foam_func,
  )

if __name__ == "__main__":
  args = make_parser()
  if args.return_func:
    run_corr_func_test(args.dataDir, args.outDir, norm=args.norm)
  else:
    run_corr_test(args.dataDir, args.outDir, norm=args.norm)

