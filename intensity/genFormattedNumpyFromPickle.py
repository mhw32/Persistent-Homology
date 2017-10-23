import os
import numpy as np
import cPickle
import argparse


def format_eagle(pickle_path, out_dir, testType='intensity', norm=False):
    with open(pickle_path) as fp:
        data = cPickle.load(fp)
    # TODO: no longer hardcoding
    cdm_stats = np.zeros([4] + list(data[0][0].shape))
    cdm_stats[-1, ...] = data[0][0]
    wdm_stats = np.zeros([4] + list(data[0][1].shape))
    wdm_stats[-1, ...] = data[0][1]

    normstr = 'norm' if norm else 'unnorm'
    np.save(os.path.join(out_dir, 'cdm_{}_stats_{}.npy'.format(testType, normstr)), cdm_stats)
    np.save(os.path.join(out_dir, 'wdm_{}_stats_{}.npy'.format(testType, normstr)), wdm_stats)


if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('picklePath', type=str, help='input folder')
    ap.add_argument('outDir', type=str, help='input folder')
    ap.add_argument('--testType', type=str, default='intensity', help='test type')
    ap.add_argument('--norm', action='store_true')

    # parse args and run stuff
    args = vars(ap.parse_args())
    format_eagle(args['picklePath'], args['outDir'],
                 args['testType'], args['norm'])