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
def example_plot(input_filenames):
    titles = ('BBS', 'CBC','FIA','Gentry','MCDB','NABC')
    
    for i in range(0,len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,i8,f8,f8", 
                           names = ['site','year','obs','pred'], delimiter = ",")
        
        ifile2 = np.genfromtxt(input_filenames1[i], dtype = "S15,i8,i8,i8,f8,f8", 
                           names = ['site','year','S','N','p','weight'], delimiter = ",")
        
        colors1 = [ '#0033cc', '#0033cc', '#339900', '#339900', '#ff6600', '#990000']
        colors2 = [ '#87a4ef', '#87a4ef', '#97ca82', '#97ca82', '#ff8c3f', '#b25353']
        site = ifile["site"]
        obs = ifile["obs"]   
        pred = ifile["pred"] 
        S = ifile2["S"]
        site2 = ifile2["site"]
        max_site = site2[S == max(S)][0]
        min_site = site2[S == min(S)][0]
        max_obs = obs[site == max_site]
        max_pred = pred[site == max_site]
        min_obs = obs[site == min_site]
        min_pred = pred[site == min_site]
        x1 = range(1,max(S) + 1)
        x2 = range(1,min(S) + 1)
        plt.figure()#facecolor='#d7dfff')
        #plt.subplot(111, axisbg='#d7dfff')
        plt.plot(x1, np.log(max_pred), '-', color = colors1[i], linewidth = 2)
        plt.plot(x1, np.log(max_obs), 'ow', markersize = 10, 
                 markeredgecolor = colors1[i], markeredgewidth = 2)        
        plt.plot(x2, np.log(min_pred), '-', color = colors2[i], linewidth = 2)
        plt.plot(x2, np.log(min_obs), 'ow', markeredgecolor = colors2[i], 
                 markersize = 10, markeredgewidth = 2)
        plt.xlim(-1.25, max(x1) + 2)
        plt.ylim(-0.15, np.log(max(max(max_obs), max(max_pred), max(min_pred), max(min_pred)))+ 0.25)
        plt.xlabel('Rank', fontsize = '22')
        plt.xticks(fontsize = '20')
        plt.yticks(fontsize = '20')
        plt.ylabel('log(Abundance)', fontsize ='22')
        plt.savefig(titles[i] + '_example.png', dpi=400, facecolor='w', edgecolor='w', 
                bbox_inches = 'tight', pad_inches=0)

example_plot(input_filenames)

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
    plt.figure(1) 
    n = len(input_filenames)
    colors = [ '1.0', '0.8', '0.6', '0.4','0.2', '0']
    titles = ['BBS', 'CBC', 'FIA', 'Gentry', 'MCDB', 'NABC']
    for i in range(0, n):
        input_filename = input_filenames[i]
        width = round(1.0/(3 + n * 3), 2)
        left = [(width * (i + 1)), (width * (i + n + 2)), (width * (i + n + 9))]
        mete_sads.plot_weights(input_filename, data = 'percent', left = left, 
                     color = colors[i], width = width)
    
    plt.ylabel('Percentage of sites')
    plt.xlim((width/2), (width*(3.5 + n * 3)))
    plt.xticks((((n/2 + 1) * width), (11 * width),(18 * width)), 
               ('Log-normal', 'Indeterminate', 'Log-series') )
    plt.legend(('BBS', 'CBC', 'FIA', 'Gentry', 'MCDB', 'NABC'), loc = 'upper left')
    
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
    pr = (ifile2["p"])
    S = ifile2["S"]
    usites = list(set(ifile["site"]))
    
    sites = ifile["site"]
    sites2 = ifile2["site"]
    obs_ab = ifile["obs"]
    pred_ab = ifile["pred"]

    plt.subplot(3,2,i+1)
    mete_sads.plot_avg_deviation_from_logseries(ifile['site'], ifile['obs'],
                                                ifile['pred'])
    plt.subplots_adjust(left=0.2, bottom=0.12, right=0.8, 
                        top=0.92, wspace=0.29,hspace=0.21)