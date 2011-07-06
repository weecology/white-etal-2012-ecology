from mpl_toolkits.basemap import Basemap
import macroeco_distributions as md
import mete
import mete_sads
import macroeco
import matplotlib.pyplot as plt
import numpy as np

#workdir = raw_input('Enter the directory where the data files are located:\n')
workdir = '/home/ethan/Dropbox/Research/MaxEnt/Code/data/'

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

input_filenames2 = (workdir + 'bbs_lat_long.csv',
                   workdir + 'cbc_lat_long.csv',
                   workdir + 'fia_lat_long.csv',
                   workdir + 'gentry_lat_long.csv',
                   workdir + 'mcdb_lat_long.csv',
                   workdir + 'naba_lat_long.csv')

colors = [ '#87a4ef', '#0033cc', '#97ca82', '#339900','#ff6600', '#990000']

#figure 1a
##aea383 - brown - cdaa7d - 7c5e4e - #d5a76b - #325582 - dark blue
def map_sites(input_filenames, markers = ['o'], colors = ['b', 'r', 'g', 'y', 'c'],
              markersizes=3):
    """Generate a world map with sites color-coded by database"""
    
    map = Basemap(projection='merc',llcrnrlat=-57,urcrnrlat=71,\
                llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='l')
    
    map.drawcoastlines(linewidth = 1.25)

    for i in range(0, len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "f8,f8", 
                                   names = ['lat','long'], delimiter = ",")
        lats = ifile["lat"]
        longs = ifile["long"]  
    
        x,y = map(longs,lats)
        map.plot(x,y, ls = '', marker = markers[i], markerfacecolor = colors[i], 
                 markeredgewidth = 0.25, markersize = markersizes)
    
    plt.savefig('map.png', dpi=400, facecolor='w', edgecolor='w', 
                bbox_inches = 'tight', pad_inches=0)
    
map_sites(input_filenames2, markers = ['o','o','s','s','D','v'], markersizes = 4,
                    colors = [ '#87a4ef', '#0033cc', '#97ca82', '#339900','#ff6600', '#990000'])

#figure 1b    
def map_sites_inset(input_filenames, markers = ['o'],
                    colors = ['b', 'r', 'g', 'y', 'c'], markersizes=4):
    """Generate a US map with sites color-coded by database"""
    
    map = Basemap(projection='merc',llcrnrlat=24.5,urcrnrlat=49,\
                llcrnrlon=-125,urcrnrlon=-69,lat_ts=50,resolution='l')
    
    map.drawcoastlines(linewidth=1.5)
    map.drawcountries(linewidth=1.5)
    map.drawmapboundary()     
    
    for i in range(0, len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "f8,f8", 
                                   names = ['lat','long'], delimiter = ",")
        lats = ifile["lat"]
        longs = ifile["long"]  
    
        x,y = map(longs,lats)
        map.plot(x,y, ls = '', marker = markers[i], markerfacecolor = colors[i], 
                 markeredgewidth = 0.25, markersize = markersizes)
    
    plt.savefig('map_inset.png', dpi=400, facecolor='w', edgecolor='w', 
                bbox_inches = 'tight', pad_inches=0)
    
map_sites_inset(input_filenames2, markers = ['o','o','s','s','D','v'], markersizes = 5,
                    colors = [ '#87a4ef', '#0033cc', '#97ca82',
                               '#339900','#ff6600', '#990000'])

#figure 1c
def example_plot(workdir, dataset_code, site_id, color, axis_limits):
    """Generate an example SAD plot for the map figure"""    
    obs_pred_data = np.genfromtxt(workdir + dataset_code + '_obs_pred.csv',
                          dtype = "S15,i8,f8,f8",
                          names = ['site','year','obs','pred'], delimiter = ",")    
    site = obs_pred_data["site"]
    obs = obs_pred_data["obs"]   
    pred = obs_pred_data["pred"]
    site_obs_ab = obs[site==site_id]
    site_pred_ab = pred[site==site_id]
    rank_obs, relab_obs = macroeco.get_rad_data(site_obs_ab)
    rank_pred, relab_pred = macroeco.get_rad_data(site_pred_ab)
    plt.figure(figsize=(2,2))
    plt.semilogy(rank_obs, relab_obs, 'o', markerfacecolor='none', markersize=8, 
             markeredgecolor=color, markeredgewidth=1.5)
    plt.semilogy(rank_pred, relab_pred, '-', color='black', linewidth=2)
    plt.axis(axis_limits)
    plt.xticks(fontsize = '10')
    plt.yticks(fontsize = '10')
    plt.savefig(dataset_code + '_example.png', dpi=400, facecolor='w', edgecolor='w', 
            bbox_inches = 'tight', pad_inches=0.2)

dataset_codes = ['bbs', 'cbc', 'fia', 'gentry', 'mcdb', 'naba']
site_ids = ['17220', 'L13428', '211131000022', '84', '1353', 'TX_Still_ollow']
axis_limits = [[0, 16, 10 ** -4, 1], [-5, 115, 10 ** -5, 1],
               [0, 15, 10 ** -2, 1], [-10, 225, 10 ** -3, 1], 
               [0, 11, 10 ** -2, 1], [-2, 44, 10 ** -3, 1]]
for i, dataset in enumerate(dataset_codes):
    example_plot(workdir, dataset, site_ids[i], colors[i], axis_limits[i])

#figure 2
def var_plot(input_filenames, radius=2):
    """Multiple obs-predicted plotter"""
    #TODO Cleanup transformations using dictionary based approach and error
    #     checking for cases where a provided transformation is undefined
    #TODO Generalize to different numbers of subplots
    titles = ('BBS', 'CBC','FIA','Gentry','MCDB','NABA')
    
    for i in range(0,len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,i8,f8,f8", 
                           names = ['site','year', 'obs','pred'], delimiter = ",")
        site = ((ifile["site"]))
        obs = ((ifile["obs"]))    
        pred = ((ifile["pred"])) 
        
        axis_min = 0.5 * min(obs)
        axis_max = 2 * max(obs)
        axis_scale = 1
        
        macroeco.plot_color_by_pt_dens(pred, obs, radius, loglog=axis_scale, 
                                       plot_obj=plt.subplot(3,2,i+1))        
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        if i == 2:
            plt.ylabel('Observed Abundance')
        elif i == 4:
            plt.xlabel('Predicted Abundance')
        elif i == 5:        
            plt.xlabel('Predicted Abundance')
        plt.subplots_adjust(left=0.2, bottom=0.12, right=0.8, top=0.92, 
                                wspace=0.29, hspace=0.21)
        #plt.text(0.9, axis_max - ((axis_max-axis_min)/2 - axis_max*0.1), titles[i])   
        
    plt.savefig('fig2.png', dpi=400, bbox_inches = 'tight', pad_inches=0) 
    
var_plot(input_filenames, radius = 3)

#figure 3
input_filenames = (workdir + 'bbs_dist_test.csv',
                   workdir + 'cbc_dist_test.csv',
                   workdir + 'fia_dist_test.csv',
                   workdir + 'gentry_dist_test.csv',
                   workdir + 'mcdb_dist_test.csv',
                   workdir + 'naba_dist_test.csv')

def cross_taxa_weight_plot (input_filenames):
    """Plot histogram of weights across taxa
    
    Keyword arguments:
    input_filenames -- list of file names to be processed
    
    """
    n = len(input_filenames)
    indiv_dataset_bars = []
    fig = plt.figure()
    plot_obj = fig.add_subplot(111)
    for i, input_filename in enumerate(input_filenames):
        width = round(1.0/(3 + n * 3), 2)
        left = [(width * (i + 1)), (width * (i + n + 2)), (width * (i + n + 9))]
        ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8,i8,f8,f8", 
                              names = ['site', 'year', 'S', 'N', 'p', 'weight'],
                              delimiter = ",")
        weights = ((ifile["weight"]))
        bins = [0, 0.33333, 0.66667, 1]
        cts = np.histogram(weights, bins = bins)
        height = cts[0] * 100 / sum(cts[0])
        indiv_dataset_bars.append(plot_obj.bar(left, height, width,
                                               color=colors[i]))
    plt.ylabel('Percentage of sites')
    plt.xlim((width/2), (width*(3.5 + n * 3)))
    plt.xticks((((n/2 + 1) * width), (11 * width),(18 * width)), 
               ('Log-normal', 'Indeterminate', 'Log-series'))
    plt.legend((indiv_dataset_bars[0][0], indiv_dataset_bars[1][0], 
                indiv_dataset_bars[2][0], indiv_dataset_bars[3][0], 
                indiv_dataset_bars[4][0], indiv_dataset_bars[5][0]),
                ('BBS', 'CBC', 'FIA', 'Gentry', 'MCDB', 'NABC'),
                loc = 'upper left')
    plt.show()
    
cross_taxa_weight_plot (input_filenames)

#Text results associated with Figure 3
total_logser_count_onethird = 0
total_logser_count_twothird = 0
total_count = 0
for filename in input_filenames1:
    ifile = np.genfromtxt(filename, dtype = "S15,i8,i8,i8,f8,f8", 
                       names = ['site', 'year', 'S', 'N', 'p', 'weight'],
                       delimiter = ",")
    weights = ((ifile["weight"]))
    total_logser_count_onethird += len(weights[weights >= 0.333333])
    total_logser_count_twothird += len(weights[weights >= 0.666667])
    total_count += len(weights)    
total_logser_prop_better = float(total_logser_count_twothird) / total_count
total_logser_prop_equiv = float(total_logser_count_onethird) / total_count
print total_logser_prop_better, total_logser_prop_equiv

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
plt.xlabel('Preston Bin', fontsize=18)
plt.ylabel('Deviation (% of site richness)', fontsize=18)
plt.axis([0.1, 17, -8, 3.5])
plt.legend(('BBS','CBC','FIA','Gentry','MCDB','NABC'), 'lower right')