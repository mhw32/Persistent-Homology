# -- Figure 13 : heatmap analysis (correlation) --
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

cdm_corr_funcs = np.load('intermediate/fig_13_cdm_corr.npy')
wdm_corr_funcs = np.load('intermediate/fig_13_wdm_corr.npy')

# handle max margin first
max_idx = 13 # 36
my_cdm_corr = cdm_corr_funcs[-1, max_idx, :]
my_wdm_corr = wdm_corr_funcs[-1, max_idx, :]

plt.figure()
plt.plot(my_cdm_corr, 
         linewidth=3, 
         color='#F88971', 
         label='CDM')
plt.plot(my_wdm_corr, 
         linewidth=3, 
         color='#3FDFDA', 
         label='WDM')
plt.tick_params(labelsize=18)
plt.xlabel('Grid Sequence', fontsize=18)
plt.ylabel('Correlation Func.', fontsize=18)
plt.legend(loc='upper right', fontsize=18)
plt.tight_layout()
# plt.show()
plt.savefig('figure_13_max_margin_corr.pdf')

# handle min margin first
min_idx = 63
my_cdm_corr = cdm_corr_funcs[-1, min_idx, :]
my_wdm_corr = wdm_corr_funcs[-1, min_idx, :]

plt.figure()
plt.plot(my_cdm_corr, 
         linewidth=3, 
         color='#F88971', 
         label='CDM')
plt.plot(my_wdm_corr, 
         linewidth=3, 
         color='#3FDFDA', 
         label='WDM')
plt.tick_params(labelsize=18)
plt.xlabel('Grid Sequence', fontsize=18)
plt.ylabel('Correlation Func.', fontsize=18)
plt.legend(loc='upper right', fontsize=18)
plt.tight_layout()
# plt.show()
plt.savefig('figure_13_min_margin_corr.pdf')