import subprocess
from os import listdir
from os.path import isfile, join, splitext, abspath

ORIGINAL_FIG_DIR = abspath('figs')
COMPRESSED_FIG_DIR = abspath('compressed_figs')

onlyfiles = [f for f in listdir(ORIGINAL_FIG_DIR) if isfile(join(ORIGINAL_FIG_DIR, f))]
ignorelist = ['.DS_Store']

for fname in onlyfiles:
    if fname in ignorelist:
        print('ignoring {}...'.format(fname))
        continue

    print('compressing {}...'.format(fname))
    _name, _ext = splitext(fname)

    if _ext == '.pdf':
        cmd = 'bash shrinkpdf.sh {}/{} {}/{}'.format(
            ORIGINAL_FIG_DIR, fname, COMPRESSED_FIG_DIR, fname)
        subprocess.call(cmd, shell=True)

    elif _ext == '.png':
        cmd = 'optipng {}/{} -out {}/{}'.format(
            ORIGINAL_FIG_DIR, fname, COMPRESSED_FIG_DIR, fname)
        subprocess.call(cmd, shell=True)
