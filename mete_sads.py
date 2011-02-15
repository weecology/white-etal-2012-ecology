"""Project-specific Code for Testing METE's SAD Predictions

Required input = Abundances per species per site for one sampling period
    
All data queries used can be found in MaxEnt/trunk/data:
    BBS_data_query
    CBC_data_query
    Gentry_data_query
    
Requires mpl_toolkits module. Install by following these instructions:
http://matplotlib.sourceforge.net/basemap/doc/html/users/installing.html
        
"""

from mpl_toolkits.basemap import Basemap
import macroeco_distributions as md
import mete
import csv
import macroeco
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import weestats
import pickle

def run_test(input_filename, output_filename1, output_filename2, cutoff = 9):
    """Use data to compare the predicted and empirical SADs and get results in csv files
    
    Keyword arguments:
    input_filename -- path to file that has raw data in the format: 
                        'site','year','sp','ab'
    output_filename1 -- file that will store the pred and observed species-level 
                        abundances for each site in the input file
    output_filename2 -- file that will store the p-values and weights from 
                        dist_test for each site in the input file
    cutoff      --  minimum number of species required to run - 1.
    
    """
    
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,S9,i8", 
                       names = ['site','year','sp','ab'], delimiter = ",")
    
    usites = np.sort(list(set(ifile["site"])))
    
    f1 = csv.writer(open(output_filename1,'ab'))
    f2 = csv.writer(open(output_filename2,'ab'))
    
    for i in range(0, len(usites)):
        subsites = ifile["site"][ifile["site"] == usites[i]]
        S = len(subsites)
        subab = np.sort(ifile["ab"][ifile["site"] == usites[i]])
        N = sum(subab)
        if S > cutoff:
            # Generate predicted values and p (e ** -lambda_sad) based on METE:
            mete_pred = mete.get_mete_sad(int(S),int(N))
            pred = np.array(mete_pred[0])
            p = mete_pred[1]
            subab = np.sort(subab)[::-1]
            # Calculate Akaike weight of log-series:
            L_logser = md.logser_ll(subab, p)        
            mu = np.mean(np.log(subab))
            sigma = np.std(np.log(subab))
            L_pln = md.pln_ll(mu,sigma,subab)        
            k1 = 1
            k2 = 2    
            AICc_logser = weestats.AICc(k1, L_logser, S)
            AICc_pln = weestats.AICc(k2, L_pln, S)
            weight = weestats.aic_weight(AICc_logser, AICc_pln, S, cutoff = 4) 
            #save results to a csv file:
            results = ((np.column_stack((subsites, subab, pred))))
            results2 = ((np.column_stack((usites[i], S, N, p, weight))))
            f1.writerows(results)
            f2.writerows(results2)
            
def plot_pred_obs(input_filename, title = ''): 
    """use output from run_test to plot observed vs. predicted abundances"""
    
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8", 
                       names = ['site','obs','pred'], delimiter = ",")
    
    pred = ((ifile["pred"]))
    obs = ((ifile["obs"]))
    
    plt.figure()
    macroeco.plot_color_by_pt_dens(pred, obs, 5, loglog=1)
    plt.title(title)
    plt.xlabel('Predicted abundances')
    plt.ylabel('Observed abundances')
    plt.show()    
    
def plot_weights(input_filename, data = 'raw', left = [0, 0.4, 0.8], 
                 width = 0.2, color = 'b', title = ''): 
    """use output from run_test to plot frequency distribution of Akaike weights"""
    
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8,f8,f8", 
                       names = ['site','S','N','p','weight'], delimiter = ",")
    
    weights = ((ifile["weight"]))
    weights=weights[weights>=0]
    bins = [0, 0.4, 0.6, 1]
    cts = np.histogram(weights, bins = bins)
    
    if data == 'raw':
        height = cts[0]  
    else:
        height = cts[0] * 100 / sum(cts[0])
    
    plot_obj = plt.subplot(111)
    plot_obj.bar(left, height, width, color = color)
    plot_obj.set_title(title)
    plot_obj.set_ylabel('Number of sites') 
    
    return plot_obj

def cross_taxa_weight_plot (input_filenames):
    """Plot histogram of weights across taxa
    
    Keyword arguments:
    input_filenames -- list of file names to be processed
    
    """     
    plt.figure(1) 
    n = len(input_filenames)
    colors = ['b', 'r', 'k', 'g', '0.75']
    
    for i in range(0, n):
        input_filename = input_filenames[i]
        width = round(1.0/(3 + n * 3), 2)
        left = [(width * (i + 1)), (width * (i + n + 2)), (width * (i + n + 8))]
        plot_weights(input_filename, data = 'percent', left = left, 
                     color = colors[i], width = width)
    
    plt.ylabel('Percentage of sites')
    # TO DO: figure out universal means of determining xtick locations
    plt.xlim((width/2), (width*(3.5 + n * 3)))
    plt.xticks((((n/2 + 1) * width), (((3 + n * 3)/2 + 0.75) * width),
                (((n * 3) + 0.5) * width)), 
               ('Log-normal', 'Indeterminate', 'Log-series') )
    #plt.xticks(((width * (n/2)), (width * (n/2) + 3.5), (width * (n/2) + 11)),
    #           ('Log-normal', 'Indeterminate', 'Log-series'))
    # TO DO: figure out how to include a color-coded legend: 
    plt.legend(('CBC', 'BBS', 'MCDB', 'FIA', 'Gentry'), loc = 'upper left')
    plt.show()
    
def rare_sp_count (input_filename, abundance_class):
    """Count and plot number of species observed and predicted in a given abundance class
    
    Keyword arguments:
    input_filename -- name of file containing observed and predicted abundances 
    in the format ['site', 'obs', 'pred'], as output from run_test
    abundance_class -- singleton, doubleton, rare(n <=10), or dominant
    
    """
    
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8", 
                       names = ['site','obs','pred'], delimiter = ",")
    
    site = ((ifile["site"]))    
    usites = list(set(site))  
    
    pred_class = []
    obs_class = []
    for i in range (0, len(usites)):
        pred = ifile["pred"][ifile["site"] == usites[i]]
        obs = ifile["obs"][ifile["site"] == usites[i]]
        if abundance_class == 'singleton':
            subpred = len(pred[pred == 1])
            subobs = len(obs[obs == 1])
        elif abundance_class == 'doubleton':
            subpred = len(pred[pred == 2])
            subobs = len(obs[obs == 2])
        elif abundance_class == 'dominant':
            subpred = max(pred)
            subobs = max(obs)
        elif abundance_class == 'rare':
            subpred = len(pred[pred <= 10])
            subobs = len(obs[obs <= 10])
        pred_class.append(subpred)
        obs_class.append(subobs)
        
    return(pred_class, obs_class)

def multi_taxa_rare_sp_plot(input_filenames):
    """Generate 2 x 2 paneled figure of pred vs obs numbers of rare species 
    for all 5 databases
    
    """
    titles = ('BBS', 'CBC','Mammals','Trees')
    for i, filename in enumerate(input_filenames):
        results = rare_sp_count (input_filenames[i], 'rare')
        pred = results[0]
        obs = results[1]      
        plot_obj = plt.subplot(2,2,i+1)
        macroeco.plot_color_by_pt_dens(np.array(pred), np.array(obs), 2, 0, plot_obj)
        plotmax = max(max(pred),max(obs)) + 5;
        plt.plot([0, plotmax], [0, plotmax], 'k-')
        plt.xlim(0, plotmax)
        plt.ylim(0, plotmax)
        plt.title(titles[i])
        if i == 0:
            plt.ylabel('Observed number of species')            
        elif i == 2:
            plt.xlabel('Predicted number of species')
            plt.ylabel('Observed number of species')  
        elif i == 3:
            plt.xlabel('Predicted number of species')            

    plt.show()

def ab_class_test_plot(input_filename):
    """Regress number of species predicted vs. observed in a given abundance class
    
    Keyword arguments:
    input_filename -- name of file containing observed and predicted abundances 
    in the format ['site', 'obs', 'pred'], as output from run_test
    
    """
    abundance_classes = ['singleton', 'doubleton', 'rare', 'dominant']

    regr_results = []
    for i in range (0, len(abundance_classes)):
        results = rare_sp_count (input_filename, abundance_classes[i])
        pred = results[0]
        obs = results[1] 
        slope, intercept, r_value, p_value, std_err = stats.linregress(pred, obs)
        results2 = ((np.column_stack((slope, intercept, r_value, p_value, std_err))))
        regr_results.append(results2)
        # had to change macroeco.py to accept plot_obj to generate subplots   
        plot_obj = plt.subplot(2,2,i+1)
        if i < 3:
            macroeco.plot_color_by_pt_dens(np.array(pred), np.array(obs), 1, 0, plot_obj)
            plt.plot([0,max(max(pred),max(obs))+ 1], 
                     [0,max(max(pred),max(obs)) + 1], 'k-')
            plt.xlim(0, max(max(pred),max(obs)) + 1)
            plt.ylim(0, max(max(pred),max(obs)) + 1)
            plt.title(abundance_classes[i])
            plt.xlabel('Predicted number of species')
            plt.ylabel('Observed number of species')
        else:
            macroeco.plot_color_by_pt_dens(np.array(pred), np.array(obs), 1, 1, plot_obj)
            plt.loglog([1, (max(pred)) * 2], [1, (max(pred)) * 2], 'k-')
            plt.xlim(1, max(max(pred),max(obs)) * 2)
            plt.ylim(1, max(max(pred),max(obs)) * 2)
            plt.title(abundance_classes[i])
            plt.xlabel('Predicted abundance')
            plt.ylabel('Observed abundance')
        r2 = ('r2 = ' + str(round(r_value**2, 2)))
        b = ('y = ' + str(round(slope, 2)) + 'x + ' + str(round(intercept)))
        plt.annotate(b, xy=(-10, 10), xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=14)
        plt.annotate(r2, xy=(-10, 30), xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=14)

    plt.show()   
    return(regr_results)

def multi_taxa_conf_hulls(input_filenames, radius, conf_interval, logscale=0):
    colors = ['r', 'b', 'k', 'g', 'c']
    plotmax = 0
    for i, filename in enumerate(input_filenames):
        infile = np.genfromtxt(filename, dtype = "S15,i8,i8", 
                       names = ['site','obs','pred'], delimiter = ",")
        hull_points = macroeco.confidence_hull(infile['pred'], infile['obs'], radius,
                                               logscale=logscale,
                                               color=colors[i], alpha=0.75-float(i)/10)
        plotmax = max([plotmax, np.max(hull_points)])
    plt.loglog([0.5, plotmax * 2], [0.5, plotmax * 2], 'k-')
    plt.xlim(.5, plotmax * 2)
    plt.ylim(.5, plotmax * 2)
    plt.xlabel('Predicted abundance')
    plt.ylabel('Observed abundance')
    plt.show()

def evar_pred_obs(input_filenames, output_filenames):
    """Calculate Evar for observed and predicted SADs and save to file
    
    Keyword arguments:
    input_filename -- use file output from run_test (db_obs_pred.csv)
    output_filename -- name of file in which to store Evar results for all sites
    
    """
    
    for i in range(0,len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,i8,i8", 
                           names = ['site','obs','pred'], delimiter = ",")
        
        site = ((ifile["site"]))    
        usites = list(set(site))  
        
        f1 = csv.writer(open(output_filenames[i],'a'))
    
        for i in range (0, len(usites)):
            pred = ifile["pred"][ifile["site"] == usites[i]]
            obs = ifile["obs"][ifile["site"] == usites[i]]
            evar_obs = macroeco.e_var(obs)
            evar_pred = macroeco.e_var(pred)
            results = (usites[i], evar_obs, evar_pred)
            f1.writerow(results)  
            
def var_plot(input_filenames, radius=2, transform='no'):
    """Multiple obs-predicted plotter"""
    #TODO Cleanup transformations using dictionary based approach and error
    #     checking for cases where a provided transformation is undefined
    #TODO Generalize to different numbers of subplots
    titles = ('CBC', 'BBS', 'Gentry', 'MCDB', 'FIA')
    
    for i in range(0,len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,f8,f8", 
                           names = ['site','obs','pred'], delimiter = ",")
    
        obs = ((ifile["obs"]))    
        pred = ((ifile["pred"])) 
        
        if transform == 'arcsin':
            obs_trans = np.arcsin(np.sqrt(obs))
            pred_trans = np.arcsin(np.sqrt(pred))
            axis_min = 0
            axis_max = 1
            axis_scale = 0
        elif transform == 'log10':
            obs_trans = np.log10(obs)
            pred_trans = np.log10(pred)
            axis_min = 0.5 * min(obs)
            axis_max = 2 * max(obs)
            axis_scale = 1
        else:
            obs_trans = obs
            pred_trans = pred
            axis_min = min(obs) - 1
            axis_max = max(obs) + 1
            axis_scale = 0
        slope, intercept, r_value, p_value, std_err = stats.linregress(pred_trans,
                                                                       obs_trans)
            
        plt.subplot(2,2,i+1)
        macroeco.plot_color_by_pt_dens(pred, obs, radius, loglog=axis_scale, 
                                       plot_obj=plt.subplot(2,2,i+1))        
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        if i == 0:
            plt.ylabel('Observed')
        elif i == 2:
            plt.xlabel('Predicted')
            plt.ylabel('Observed')
        elif i == 3:        
            plt.xlabel('Predicted Evar')
        plt.title(titles[i])
        r2 = ('r2 = ' + str(round(r_value**2, 2)))
        b = ('y = ' + str(round(slope, 2)) + 'x + ' + str(round(intercept)))
        plt.annotate(b, xy=(-10, 10), xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=14)
        plt.annotate(r2, xy=(-10, 30), xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=14)
        
    plt.show() 
    
def SN_diff_plot(input_filenames):

    titles = ('CBC', 'BBS', 'Gentry', 'MCDB', 'FIA')
    
    for i in range(0,len(input_filenames)/2):
        ifile = np.genfromtxt(input_filenames[i], dtype = "S15,f8,f8", 
                           names = ['site','obs','pred'], delimiter = ",")
        
        ifile2 = np.genfromtxt(input_filenames[i + len(input_filenames)/2], 
                               dtype = "S15,i8,i8,f8,f8", 
                               names = ['site','S','N','p','weight'], delimiter = ",")
        
        usites = ((ifile["site"]))
        
        diffs=[]
        xvar=[]
        for a in range(0, len(usites)):
            pred = ifile["pred"][ifile["site"] == usites[a]]
            obs = ifile["obs"][ifile["site"] == usites[a]]
            S = ifile2["S"][ifile2["site"] == usites[a]]
            N = ifile2["N"][ifile2["site"] == usites[a]]
            if len(S) > 0:
                var_diff = pred - obs
                S_over_N = float(S) / float(N)
                diffs.append(var_diff)
                xvar.append(S_over_N)
        
        plt.subplot(2,2,i+1)
        plt.xscale('log')
        plt.yscale('log')
        plt.scatter(xvar, diffs)        
        plt.plot([0, 1],[0, 1], 'k-')
        plt.xlim(min(xvar),10 ** -0.1)
        plt.ylim(min(diffs),max(diffs))
        if i == 0:
            plt.ylabel('Predicted - Observed')
        elif i == 2:
            plt.xlabel('S/N')
            plt.ylabel('Predicted - Observed')
        elif i == 3:        
            plt.xlabel('S/N')
        plt.title(titles[i])
        
    plt.show() 
    
def map_sites(input_filenames, markers = ['o'], colors = ['b', 'r', 'g', 'y', 'c'],
              markersizes=3):
    """Generate a world map with sites color-coded by database"""
    
    map = Basemap(projection='merc',llcrnrlat=-50,urcrnrlat=70,\
                llcrnrlon=-175,urcrnrlon=175,lat_ts=20,resolution='c')
    
    # bluemarble(ax=None, scale=None, **kwargs) display blue marble image 
    # (from http://visibleearth.nasa.gov) as map background. Default image size is 
    # 5400x2700, which can be quite slow and use quite a bit of memory. The scale 
    # keyword can be used to downsample the image (scale=0.5 downsamples to 2700x1350).
    map.bluemarble(scale=0.75)    
    
    # 20 degree graticule.
    map.drawparallels(np.arange(-80,81,20), linewidth = 0.25)
    map.drawmeridians(np.arange(-180,180,20), linewidth = 0.25)       
    
    for i in range(0, len(input_filenames)):
        ifile = np.genfromtxt(input_filenames[i], dtype = "f8,f8", 
                                   names = ['lat','long'], delimiter = ",")
        lats = ifile["lat"]
        longs = ifile["long"]  
    
        x,y = map(longs,lats)
        map.plot(x,y, ls = '', marker = markers[i], markerfacecolor = colors[i], 
                 markeredgewidth = 0.25, markersize = markersizes)
    
    plt.savefig('map.png', dpi=400, edgecolor='k', bbox_inches = 'tight', pad_inches=0)
    
def sim_null(S0, N0, lib_lambda):
    """Abundances simulated from a discrete uniform and associated METE predictions"""
    N_sim = sorted(np.random.random_integers(1, (2 * N0 - S0) / S0, S0), reverse = True)
    NS_ratio = round(sum(N_sim) / S0, 4) 
    if NS_ratio not in lib_lambda:
        lib_lambda[NS_ratio] = mete.get_lambda_sad(S0, sum(N_sim))
    N_pred = mete.get_mete_sad(S0, sum(N_sim), lib_lambda[NS_ratio])[0]     
    return N_sim, N_pred

def compare_null(obs, pred, Niter):
    """Compare R^2 of observed and simulated abundances"""
    r2_sim = []
    i = 0
    while i < Niter:
        r2_sim.append(sim_null(len(obs), sum(obs)))
        i += 1
    r2_sim_avg = sum(r2_sim) / len(r2_sim)
    r2_obs = macroeco.obs_pred_rsquare(obs, pred)
    return r2_sim_avg, r2_obs

def create_null_dataset(input_filename, output_filename, lib_filename, Niter):
    """Create list of R^2 values for simulated observed vs. predicted 
    abundance relationships for a dataset.
    
    iter: number of simulations
    input_filename: same format as output_filename2 from run_test
    output_filename: 1 column - R^2 for simulated data, one value per iteration
    
    """
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8", 
                       names = ['site','S','N'], delimiter = ",")  
    site = sorted(list(set(ifile['site'])))
    result = open(output_filename, 'ab')
    out = csv.writer(result, dialect = 'excel')
    lib_file = open(lib_filename, 'r')
    lib_lambda = pickle.load(lib_file)
    lib_file.close()
    i = 0
    while i < Niter:
        sim_obs = list()
        sim_pred = list()        
        for isite in site:
            idata = ifile[ifile['site'] == isite]
            iS = idata['S']
            iN = idata['N']
            ires = sim_null(iS, iN, lib_lambda)
            sim_obs.extend((ires[0]))
            sim_pred.extend((ires[1]))
        r2 = macroeco.obs_pred_rsquare(np.array(sim_obs), np.array(sim_pred))
        results = ((np.column_stack((i, r2))))
        out.writerows(results)
        i += 1
    lib_output = open(lib_filename, 'w')
    pickle.dump(lib_lambda, lib_output)
    lib_output.close()
    result.close()