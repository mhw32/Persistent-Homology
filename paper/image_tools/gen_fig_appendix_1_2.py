import sys, numpy as np
sys.path.append('/Users/mikewu/Desktop/Research/persist-homology/')
import sub_parse
reload(sub_parse)
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

# -----------------------------------------------------------------------------------------------

def safelog10(B, noise=1e-5):
    B[B == 0] += noise
    return np.log10(B)

import rpy2.robjects as robjects
from rpy2.robjects.numpy2ri import numpy2ri
import numpy as np

def rds_to_np(Rfile):
    ''' Convert .RData to be able to load 
        the array into Python.

        Rfile := (str) location of file to 
        translate to Python.
    ''' 
    raw = robjects.r['readRDS'](Rfile)
    return raw

def get_align_dict(keys, step=0.01):
    alignment_hash = {}
    mid_pt = len(allkeys) / 2
    base_val = mid_pt * -step
    for i, key in enumerate(allkeys):
        alignment_hash[key] = base_val
        base_val += step
    return alignment_hash

# -----------------------------------------------------------------------------------------------

name = '/Users/mikewu/Desktop/Research/Cisewski-Lab/saved_states/large_set_test/results-'
norm = 'False'
base = '0.1'
paths = [name+str(i)+'-'+base+'baseNorm'+norm+'.txt' for i in range(1, 101)]
singles = ['all-silh', 'euler', 'all-euler', 'silh-euler']
doubles = ['indiv_silh', 'indiv-euler', 'contour']

resArr = np.array([sub_parse.parse(f) for f in paths])

bighash_no_norm = {}
for characteristic in singles:
    bighash_no_norm[characteristic] = sub_parse.prepare1d(resArr, characteristic)
    
for characteristic in doubles:
    for dim in [0,1,2]:
        bighash_no_norm[characteristic+'-dim-'+str(dim)] = sub_parse.prepare2d(
            resArr, characteristic, dim)

if norm == 'False':
    bighash['contour-dim-0'][-1, :] = bighash['contour-dim-0'][-2, :]

corr_no_norm = np.array(rds_to_np('/Users/mikewu/Desktop/Research/persist-homology/correlation/output/voronoi_proba_norm(0).rds'))
pi_no_norm = np.load('/Users/mikewu/Desktop/Research/persist-homology/intensity/output/proba/nonorm/ans-voronoi-nodim.npy')
wkc_no_norm = np.load('/Users/mikewu/Desktop/Research/persist-homology/intensity/output/proba/nonorm/ans-voronoi-bydim.npy')
wkc_no_norm = np.log(1 - np.exp(wkc_no_norm))
pi_no_norm = np.log(1 - np.exp(pi_no_norm))
wkc_0_no_norm = wkc_no_norm[:, 0, :]
wkc_1_no_norm = wkc_no_norm[:, 1, :]
wkc_2_no_norm = wkc_no_norm[:, 2, :]

bighash_corr_no_norm = {'corr': np.array(corr_no_norm.T)}
bighash_wik_no_norm = {'wik_0': np.array(wkc_0_no_norm.T),
                       'wik_1': np.array(wkc_1_no_norm.T),
                       'wik_2': np.array(wkc_2_no_norm.T),
                       'pi': np.array(pi_no_norm.T) }

# -----------------------------------------------------------------------------------------------

name = '/Users/mikewu/Desktop/Research/Cisewski-Lab/saved_states/large_set_test/results-'
norm = 'True'
base = '0.1'
paths = [name+str(i)+'-'+base+'baseNorm'+norm+'.txt' for i in range(1, 101)]
singles = ['all-silh', 'euler', 'all-euler', 'silh-euler']
doubles = ['indiv_silh', 'indiv-euler', 'contour']

resArr = np.array([sub_parse.parse(f) for f in paths])

bighash_yes_norm = {}
for characteristic in singles:
    bighash_yes_norm[characteristic] = sub_parse.prepare1d(resArr, characteristic)
    
for characteristic in doubles:
    for dim in [0,1,2]:
        bighash_yes_norm[characteristic+'-dim-'+str(dim)] = sub_parse.prepare2d(
            resArr, characteristic, dim)

corr_yes_norm = np.array(rds_to_np('/Users/mikewu/Desktop/Research/persist-homology/correlation/output/voronoi_proba_norm(1).rds'))   
pi_yes_norm = np.load('/Users/mikewu/Desktop/Research/persist-homology/intensity/output/proba/yesnorm/ans-voronoi-nodim.npy')
wkc_yes_norm = np.load('/Users/mikewu/Desktop/Research/persist-homology/intensity/output/proba/yesnorm/ans-voronoi-bydim.npy')
wkc_yes_norm = np.log(1 - np.exp(wkc_yes_norm))
pi_yes_norm = np.log(1 - np.exp(pi_yes_norm))
wkc_0_yes_norm = wkc_yes_norm[:, 0, :]
wkc_1_yes_norm = wkc_yes_norm[:, 1, :]
wkc_2_yes_norm = wkc_yes_norm[:, 2, :]

bighash_corr_yes_norm = {'corr': np.array(corr_yes_norm.T)}
bighash_wik_yes_norm = {'wik_0': np.array(wkc_0_yes_norm.T),
                       'wik_1': np.array(wkc_1_yes_norm.T),
                       'wik_2': np.array(wkc_2_yes_norm.T),
                       'pi': np.array(pi_yes_norm.T) }

# -----------------------------------------------------------------------------------------------

def hard_line_plot(allkeys, 
                   allticks, 
                   allcolors,
                   bighash_no_norm,
                   bighash_yes_norm,
                   custom_ylim=None,
                   save_path=None):

    xvalues = np.arange(0.1, 1, 0.1)
    align_values = get_align_dict(allkeys, 0.012)
    matplotlib.rc('xtick', labelsize=27) 
    matplotlib.rc('ytick', labelsize=27) 
    fig, ax = plt.subplots(figsize=(25,5))
    plt.ylabel(r'$\frac{p}{\|p\|}$', fontsize=60)
    plt.xlabel('PercFil', fontsize=40)

    for it in range(9):
        for k, c, t in zip(allkeys, allcolors, allticks):
            xvalue = xvalues[it] + align_values[k]
            y_value_no_norm = np.exp(bighash_no_norm[k][it])
            y_value_yes_norm = np.exp(bighash_yes_norm[k][it])
            # y_value = np.divide(y_value_yes_norm, y_value_no_norm)
            y_value = np.divide(y_value_no_norm, y_value_yes_norm)
            
            lower_error = np.percentile(y_value, 25)
            upper_error = np.percentile(y_value, 75)
            median_y_value = np.percentile(y_value, 50)

            pos = [xvalue]
            ypt = [median_y_value]
            err = [[median_y_value - lower_error], [upper_error - median_y_value]]
            plt.errorbar(pos, 
                         ypt, 
                         yerr=err, 
                         lw=6, 
                         alpha=0.4,
                         color=c, 
                         capsize=20, 
                         capthick=6)

        for k, c, t in zip(allkeys, allcolors, allticks):
            xvalue = xvalues[it] + align_values[k]
            y_value_no_norm = np.exp(bighash_no_norm[k][it])
            y_value_yes_norm = np.exp(bighash_yes_norm[k][it])
            # y_value = np.divide(y_value_yes_norm, y_value_no_norm)
            y_value = np.divide(y_value_no_norm, y_value_yes_norm)
            median_y_value = np.percentile(y_value, 50)
            
            plt.scatter([xvalue],
                        [median_y_value],
                        color=c,
                        marker='o',
                        alpha=1,
                        s=300)

    # put it here so dots are the legend
    plt.legend(
        allticks, 
        fontsize=30, 
        loc='lower left'
    )    
    ax.xaxis.grid(False)
    plt.tick_params(
        axis='both', 
        which='major', 
        labelsize=35
    )
    plt.tight_layout()
    plt.xlim(-0.1, 1)
    plt.xticks(np.arange(0, 1, 0.1))
    ax.set_xticklabels(['', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])

    if not custom_ylim is None:
        plt.ylim(custom_ylim[0], custom_ylim[1])

    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

# -----------------------------------------------------------------------------------------------

allkeys = ['euler', 'all-euler', 'indiv-euler-dim-0', 'indiv-euler-dim-1', 'indiv-euler-dim-2']
allticks = ['EC', 'EC (0:2)', 'EC (0)', 'EC (1)', 'EC (2)']
allcolors = ['#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E']
hard_line_plot(
    allkeys, 
    allticks, 
    allcolors, 
    bighash_no_norm,
    bighash_yes_norm, 
    custom_ylim=[-25, 5], 
    save_path='figure_8_all_euler_group.pdf'
)

allkeys = ['silh-euler', 'all-silh', 'indiv_silh-dim-0', 'indiv_silh-dim-1', 'indiv_silh-dim-2']
allticks = ['SIL (EC)', 'SIL (0:2)', 'SIL (0)', 'SIL (1)', 'Sil (2)']
allcolors = ['#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E']
hard_line_plot(
    allkeys, 
    allticks, 
    allcolors, 
    bighash_no_norm,
    bighash_yes_norm, 
    custom_ylim=[-10, 2], 
    save_path='figure_8_all_silhouette_group.pdf'
)

allkeys = ['contour-dim-0', 'contour-dim-1', 'contour-dim-2']
allticks = ['IK (0)', 'IK (1)', 'IK (2)']
allcolors = ['#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E']
hard_line_plot(
    allkeys, 
    allticks, 
    allcolors, 
    bighash_no_norm,
    bighash_yes_norm, 
    custom_ylim=[-3.0, 0.0], 
    save_path='figure_8_all_contour_group.pdf'
)

allkeys = ['corr']
allticks = ['CORR']
allcolors = ['#44B3C2']
hard_line_plot(
    allkeys, 
    allticks, 
    allcolors, 
    bighash_corr_no_norm,
    bighash_corr_yes_norm, 
    custom_ylim=[-14, 2], 
    save_path='figure_8_all_correlation_group.pdf'
)

allkeys   = ['wik_0', 'wik_1', 'wik_2', 'pi']
allticks  = ['WIK (0)', 'WIK (1)', 'WIK (2)', 'PI']
allcolors = ['#44B3C2', '#F1A94E', '#E45641', '#5D4C46']
hard_line_plot(
    allkeys, 
    allticks, 
    allcolors, 
    bighash_wik_no_norm,
    bighash_wik_yes_norm, 
    custom_ylim=[-3.5, 0.5], 
    save_path='figure_8_all_weighted_contour_group.pdf'
)

