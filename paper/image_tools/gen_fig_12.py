# -- Figure 12 : heatmaps --
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

def prepare(data):
    return (data[0, :1], 
            data[1, :8], 
            data[2, :27], 
            data[3, :64])

data = np.load('intermediate/fig_12_dists.npy')
data_dim = (data[:, :, 0], data[:, :, 1], data[:, :, 2])

for i, curdata in enumerate(data_dim):
    plt.figure()
    s, d, t, q = prepare(curdata)
    vmin, vmax = np.min(curdata), np.max(curdata)
    fig, axn = plt.subplots(3, 1, sharey=True, figsize=(15,3))
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    cbar_ax.tick_params(labelsize=18) 
    sns.heatmap([s], 
                ax=axn.flat[0], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                cmap="coolwarm", 
                cbar=False)
    sns.heatmap([d], 
                ax=axn.flat[1], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                yticklabels=False, 
                cmap="coolwarm", 
                cbar=False)
    sns.heatmap([q], 
                ax=axn.flat[2], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                yticklabels=False, 
                cmap="coolwarm", 
                cbar_ax=cbar_ax)
    plt.savefig("fig_12_hmap_dim{}_nonorm.pdf".format(i), 
                bbox_inches='tight')
    # plt.show()

# -- do same for normalized

data = np.load('intermediate/fig_12_dists_norm.npy')
data_dim = (data[:, :, 0], data[:, :, 1], data[:, :, 2])

for i, curdata in enumerate(data_dim):
    plt.figure()
    s, d, t, q = prepare(curdata)
    vmin, vmax = np.min(curdata), np.max(curdata)
    fig, axn = plt.subplots(3, 1, sharey=True, figsize=(15,3))
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    cbar_ax.tick_params(labelsize=18) 
    sns.heatmap([s], 
                ax=axn.flat[0], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                cmap="coolwarm", 
                cbar=False)
    sns.heatmap([d], 
                ax=axn.flat[1], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                yticklabels=False, 
                cmap="coolwarm", 
                cbar=False)
    sns.heatmap([q], 
                ax=axn.flat[2], 
                vmin=vmin, 
                vmax=vmax, 
                xticklabels=False, 
                yticklabels=False, 
                cmap="coolwarm", 
                cbar_ax=cbar_ax)
    plt.savefig("fig_12_hmap_dim{}_yesnorm.pdf".format(i), 
                bbox_inches='tight')
    # plt.show()

