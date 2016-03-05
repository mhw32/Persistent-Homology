# parse.py
# Given the outputs created by testlib.r, process them python arrays
# Include easy visualization with Seaborn.

import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
# sns settings.
sns.palplot(sns.color_palette('hls', 8))
sns.set_style('whitegrid')

def extract_name(line, header):
    s = line.find(header) + len(header) + 1
    e = line.find(":", s) - 1
    return line[s:e]

def remove_trash(content):
    bad_stuff = ['[1] ""\n', '[1] "--------------------------------"\n']
    new_content = [line.replace('\n', '') for line in content if not line in bad_stuff]
    best_content = []
    for line in new_content:
        while '[' in line:
            b = find_brackets(line)
            line = line[b+1:]
        best_content.append(line)
    return [line for line in best_content if not line == '']

def find_brackets(s):
    seen = False
    for i, l in enumerate(s):
        if l == '[':
            seen = True
        elif l == ']' and seen == True:
            break
    return(i)

# filepath: full path to output file.
# names: names of the tests stored.
def parse(filepath):
    # Open the file and read all lines.
    with open(filepath) as f:
        content = f.readlines()
    f.close()
    # Remove excess lines.
    content = remove_trash(content)
    # Look for the test headers.
    header = 'Test for'
    info = np.array([(i, extract_name(line, header)) for i, line in enumerate(content) if "Test for" in line])
    startidx, nametag = info[:, 0].astype('int'), info[:,1]
    # Add the last line.
    startidx = np.append(startidx, len(content))
    segments = zip(startidx[:-1], startidx[1:])
    # Loop through and pull out lines.
    parsed = {}
    for i,(j,k) in enumerate(segments):
        parsed[nametag[i]] = content[j+1:k]
    # Loop through the keys, and convert to arrays.
    for k in parsed:
        v = parsed[k]
        if k in ['all-silh', 'euler', 'all-euler', 'silh-euler']:
            c = list(itertools.chain(*[row.split(' ') for row in v]))
            parsed[k] = np.array(filter(lambda a: a != '', c)).astype('float')
        else: # This is for all the other types.
            parsed[k] = np.array([np.array(filter(lambda a: a != '', row.split(' '))).astype('float') for row in v])
    return parsed

def exponentiate(parsed, key):
    for r in parsed:
        r[key] = np.exp(r[key])
    return parsed

# Extra functions to plot things.
def prepare1d(resArr, key):
    obj = np.array([[r[key][i] for r in resArr] for i in range(len(resArr[0][key]))])[:-1]
    return obj

# For 2-d matrices, parse by dimension.
def prepare2d(resArr, key, dim):
    newArr = np.array([r[key][:,dim] for r in resArr])
    obj = np.array([[r[i] for r in newArr] for i in range(len(newArr[0]))])[:-1]
    return obj

def plottest(p, key, save=True):
    meanobj = [np.average(i) for i in p]
    # plot IQR instead of max min
    minobj = [np.percentile(i, 25) for i in p]
    maxobj = [np.percentile(i, 75) for i in p]
    percFils = np.arange(0.1, 1, 0.1)
    plt.figure()
    plt.xlabel('Percent Filled')
    plt.ylabel('Proba')
    plt.title(key + ' characteristic (20 iter)')
    plt.plot(percFils, meanobj, '-o')
    plt.fill_between(percFils, minobj, maxobj, facecolor='#da70d6', interpolate=True, alpha=0.2)
    plt.savefig('../saved_states/test_samples/'+key+'-test-base09.png')
    plt.show()
