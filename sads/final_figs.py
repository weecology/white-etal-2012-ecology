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

def single_var_plot(input_filenames, radius=2):
    """Plot all obs-predicted data together in a single density colored plot
    
    WARNING: this takes quit some time to run given the cost of determining 
             neighbors for such a large number of points.
    
    """
    fig = plt.figure()
    site = np.array([])
    obs = np.array([])
    pred = np.array([])
    for i in range(0,len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,i8,f8,f8", 
                           names = ['site','year', 'obs','pred'], delimiter = ",")
        site = np.concatenate([site, ifile["site"]])
        obs = np.concatenate([obs, ifile["obs"]])   
        pred = np.concatenate([pred, ifile["pred"]]) 
    axis_min = 0.5 * min(obs)
    axis_max = 2 * max(obs)
    axis_scale = 1
    macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=axis_scale, 
                                   plot_obj=plt.subplot(3,2,i+1))        
    plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
    plt.xlim(axis_min, axis_max)
    plt.ylim(axis_min, axis_max)
    plt.savefig('fig2.png', dpi=400, bbox_inches = 'tight', pad_inches=0) 
    
single_var_plot(input_filenames, radius = 3)

#figure 4
for i in range(0,6):
    ifile = np.genfromtxt(input_filenames[i], dtype = "S15,i8,i8,i8", 
                       names = ['site', 'year', 'obs','pred'], delimiter = ",")
    ifile2 = np.genfromtxt(input_filenames1[i], dtype = "S15,i8,i8,i8,f8,f8", 
                           names = ['site','year','S','N','p','weight'],
                           delimiter = ",")
    pr = ifile2["p"]
    S = ifile2["S"]
    usites = list(set(ifile["site"]))
    sites = ifile["site"]
    sites_for_p = ifile2["site"]
    obs_ab = ifile["obs"]
    mete_sads.plot_avg_deviation_from_logseries(ifile['site'], ifile['obs'],
                                                pr, sites_for_p, color=colors[i])
plt.xlabel('Abundance (binned)', fontsize=22)
plt.xticks(np.arange(1,16), ['1', '2-3', '4-7', '8-15', '16-31', '32-63',
                             '64-127', '128-255', '256-511', '512-1023',
                             '1024-2047', '2048-4095', '4096-8191',
                             '8192-16383', '16384-32767'],
           rotation=45, fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('Deviation (% of site richness)', fontsize=22)
plt.axis([0.1, 17, -8, 3.5])
plt.legend(('BBS','CBC','FIA','Gentry','MCDB','NABC'), 'lower right')
plt.savefig('fig4.png', dpi=400, bbox_inches = 'tight', pad_inches=0.1)

#Supplementary Figure 2
sim_data_files = ['bbs_sim_r2.csv', 'cbc_sim_r2.csv', 'fia_sim_r2.csv',
                  'gentry_sim_r2.csv', 'mcdb_sim_r2.csv', 'naba_sim_r2.csv']
lowerbounds = [-0.7, -2.3, 0, 0, -0.75, -0.5]
fig = plt.figure()
for i, data_file in enumerate(sim_data_files):
    sim_data = np.genfromtxt(workdir + data_file, dtype= "f8, f8",
                             names=['simulation', 'r2'], delimiter=",")
    obs_pred_data = np.genfromtxt(input_filenames[i],
                                  usecols=(2, 3), dtype = "f8,f8",
                                  names = ['obs','pred'], delimiter = ",")
    obs_r2 = macroecotools.obs_pred_rsquare(np.log10(obs_pred_data['obs']),
                                       np.log10(obs_pred_data['pred']))
    print obs_r2
    sim_kde = stats.kde.gaussian_kde(sim_data['r2'])
    xvals = np.arange(lowerbounds[i], 1, 0.01)
    yvals = sim_kde.evaluate(xvals)
    xvals = xvals[yvals > 0.000001]
    yvals = yvals[yvals > 0.000001]
    ax = fig.add_subplot(3,2,i+1)
    longdashes = [10,5]
    plot_obj, = plt.plot(xvals, yvals, 'k--', linewidth=2, color=colors[i])
    plot_obj.set_dashes(longdashes)
    plt.plot([obs_r2, obs_r2], [0, max(yvals)], color=colors[i], linewidth=2)
    plt.axis([lowerbounds[i], 1, 0, 1.1 * max(yvals)])
    
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