# -- Figure 13 : heatmap analysis (correlation) --
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

cdm_corr_funcs = np.load('intermediate/fig_13_cdm_corr.npy')
wdm_corr_funcs = np.load('intermediate/fig_13_wdm_corr.npy')

# handle max margin first
max_idx = 13 # 36
max_cdm_corr = cdm_corr_funcs[-1, max_idx, :]
max_wdm_corr = wdm_corr_funcs[-1, max_idx, :]

# handle min margin first
min_idx = 63
min_cdm_corr = cdm_corr_funcs[-1, min_idx, :]
min_wdm_corr = wdm_corr_funcs[-1, min_idx, :]

plt.figure()
plt.plot(max_cdm_corr, 
         linewidth=2, 
         color='#F88971', 
         label='CDM')
plt.plot(max_wdm_corr, 
         linewidth=2, 
         color='#3FDFDA', 
         label='WDM')
# plt.plot(min_cdm_corr,
#          # linestyle='--', 
#          linewidth=2, 
#          color='#CA6D59', 
#          label='CDM (Low)')
# plt.plot(min_wdm_corr, 
#          # linestyle='--',
#          linewidth=2, 
#          color='#2EA5A2', 
#          label='WDM (Low)')

plt.tick_params(labelsize=28)
plt.xlabel('Distance', fontsize=28)
plt.ylabel('CORR', fontsize=28)
plt.legend(loc='upper right', fontsize=28)
plt.tight_layout()
# plt.show()
plt.savefig('figure_13_maxmin_margin_corr.pdf')
