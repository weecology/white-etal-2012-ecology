"""Reproduce the published results for the METE SADs project"""

import macroeco_distributions as md
import mete
import mete_sads
import macroecotools
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from mpl_toolkits.axes_grid.inset_locator import inset_axes

#workdir = raw_input('Enter the directory where the data files are located:\n')
workdir = '/home/ethan/Dropbox/Research/MaxEnt/src/projects/sads/data/'

datasets=['bbs', 'cbc', 'fia', 'gentry', 'mcdb', 'naba']

input_filenames = (workdir + 'bbs_obs_pred.csv',
                   workdir + 'cbc_obs_pred.csv',
                   workdir + 'fia_obs_pred.csv',
                   workdir + 'gentry_obs_pred.csv',
                   workdir + 'mcdb_obs_pred.csv',
                   workdir + 'naba_obs_pred.csv')

input_filenames1 = (workdir + 'bbs_dist_test.csv', 
                   workdir + 'cbc_dist_test.csv',
                   workdir + 'fia_dist_test.csv',
                   workdir + 'gentry_dist_test.csv',
                   workdir + 'mcdb_dist_test.csv',
                   workdir + 'naba_dist_test.csv')

colors = ['#87a4ef', '#0033cc', '#97ca82', '#339900','#ff6600', '#990000']

#Rare species prediction plot

obs_pred_data = mete_sads.get_combined_obs_pred_data(datasets)
mete_sads.plot_numsp_obs_pred(obs_pred_data['site'], obs_pred_data['obs'],
                              1, 10)
plt.loglog([0.8, 300], [0.8, 300], 'k-', linewidth=2)
plt.axis([0.8, 300, 0.8, 300])
plt.xlabel('Predicted Number of Rare Species', fontsize=22)
plt.ylabel('Observed Number of Rare Species', fontsize=22)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('fig4.png', dpi=400, bbox_inches = 'tight', pad_inches=0.1)