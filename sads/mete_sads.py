"""Project-specific Code for Testing METE's SAD Predictions

This code allows the replication of analyses performed in the paper:
White, Ethan P., Katherine M. Thibault and Xiao Xiao. In review. Characterizing
species-abundance distributions across taxa and ecosystems using a simple
maximum entropy model.

Command line usage:
python mete_sads.py /path/to/data/ args

args:
empir - generate empirical fits, obs-pred data, and compare to log-normal
sims - generate simulated datasets and compare to empirical results

The main data required is abundances of each species at each site for one 
sampling period
    
All data queries used can be found in MaxEnt/trunk/projects/mete_sads_data_export.py
           
"""

from __future__ import division

import csv
import cPickle
import sys
import multiprocessing
import itertools
import os
import matplotlib.pyplot as plt
import numpy as np
from math import log, exp
from scipy import stats

import mete
import macroecotools
import macroeco_distributions as md

def import_latlong_data(input_filename):
    data = np.genfromtxt(input_filename, dtype = "f8,f8", 
                         names = ['lat','long'], delimiter = ",")
    return data

def import_obs_pred_data(input_filename):
    data = np.genfromtxt(input_filename, dtype = "S15,i8,f8,f8",
                                  names = ['site','year','obs','pred'],
                                  delimiter = ",")
    return data

def import_raw_data(input_filename):
    """Import csv data of the form: Site, Year, Species, Abundance"""
    raw_data = np.genfromtxt(input_filename, dtype = "S15,i8,S10,i8",
                          names = ['site','year','sp','ab'], delimiter = ",")   
    return raw_data

def import_obs_SN_data(input_filename):
    """Import csv data with the richness and abundance data for each site"""
    rich_abund_data = np.genfromtxt(input_filename, usecols=(2, 3),
                                    dtype="i8,i8",
                                    names=['Svals','Nvals'],
                                    delimiter=",")
    return rich_abund_data

def run_test(raw_data, dataset_name, data_dir='./data/', cutoff = 9):
    """Use data to compare the predicted and empirical SADs and get results in csv files
    
    Keyword arguments:
    raw_data : numpy structured array with 4 columns: 'site','year','sp','ab'
    dataset_name : short code that will indicate the name of the dataset in
                    the output file names
    data_dir : directory in which to store output
    cutoff : minimum number of species required to run - 1.
    
    """
    
    usites = np.sort(list(set(raw_data["site"])))
    f1 = csv.writer(open(data_dir + dataset_name + '_obs_pred.csv','wb'))
    f2 = csv.writer(open(data_dir + dataset_name + '_dist_test.csv','wb'))
    
    for i in range(0, len(usites)):
        subsites = raw_data["site"][raw_data["site"] == usites[i]]        
        subab = raw_data["ab"][raw_data["site"] == usites[i]]
        subyr = raw_data["year"][raw_data["site"] == usites[i]]
        uyr = np.sort(list(set(subyr)))
        for a in range(0, len(uyr)):
            subsubsites = subsites[subyr == uyr[a]]
            subsubab = subab[subyr == uyr[a]]
            subsubyr = subyr[subyr == uyr[a]]
            N = sum(subsubab)
            S = len(subsubsites)
            if S > cutoff:
                # Generate predicted values and p (e ** -beta) based on METE:
                mete_pred = mete.get_mete_rad(int(S), int(N))
                pred = np.array(mete_pred[0])
                p = mete_pred[1]
                p_untruncated = exp(-mete.get_beta(S, N, version='untruncated'))
                subab3 = np.sort(subsubab)[::-1]
                # Calculate Akaike weight of log-series:
                L_logser = md.logser_ll(subab3, p)
                L_logser_untruncated = md.logser_ll(subab3, p_untruncated)
                mu, sigma = md.pln_solver(subab3)
                L_pln = md.pln_ll(mu,sigma,subab3)        
                k1 = 1
                k2 = 2    
                AICc_logser = macroecotools.AICc(k1, L_logser, S)
                AICc_logser_untruncated = macroecotools.AICc(k1, L_logser_untruncated, S)
                AICc_pln = macroecotools.AICc(k2, L_pln, S)
                weight = macroecotools.aic_weight(AICc_logser, AICc_pln, S, cutoff = 4)
                weight_untruncated = macroecotools.aic_weight(AICc_logser_untruncated,
                                                         AICc_pln, S, cutoff = 4)
                #save results to a csv file:
                results = ((np.column_stack((subsubsites, subsubyr, subab3, pred))))
                results2 = ((np.column_stack((np.array(usites[i], dtype='S20'),
                                                       uyr[a], S, N, p, weight,
                                                       p_untruncated,
                                                       weight_untruncated))))
                f1.writerows(results)
                f2.writerows(results2)
    
def hist_mete_r2(sites, obs, pred):
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in sites:
        obs_site = obs[sites==site]
        pred_site = pred[sites==site]
        r2 = macroecotools.obs_pred_rsquare(obs_site, pred_site)
        r2s.append(r2)
    hist_r2 = np.histogram(r2s, range=(0, 1))
    xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
    xvals = xvals[0:len(xvals)-1]
    yvals = hist_r2[0]
    plt.plot(xvals, yvals, 'k-', linewidth=2)
    plt.axis([0, 1, 0, 1.1 * max(yvals)])

def plot_numsp_obs_pred(sites, obs_ab, min_abundance, max_abundance):
    """Observed vs. predicted plot of the number of species in an abundance range

    Drops communities where there are 0 species that occur within the range so
    that the results can be displayed on log-scaled axes. Prints the number of
    dropped communities to the screen.
    
    """
    sites = np.array(sites)
    usites = np.unique(sites)
    obs_ab = np.array(obs_ab)
    pred = []
    obs = []
    for site in usites:
        site_abs = obs_ab[sites==site]
        site_range_abundances = site_abs[(site_abs >= min_abundance) &
                                            (site_abs <= max_abundance)]
        obs_richness = len(site_range_abundances)
        pred_richness = mete.get_mete_sad(len(site_abs), sum(site_abs),
                                          bin_edges=[min_abundance,
                                                     max_abundance + 1])
        obs.append(obs_richness)
        pred.append(pred_richness)
    pred = np.array(list(itertools.chain.from_iterable(pred)))
    obs = np.array(obs)
    obs_pred_data = np.column_stack((obs, pred))
    np.savetxt('temp_sp_obs_pred_data', obs_pred_data)
    num_dropped_communities = len(obs[obs==0])
    pred = pred[obs > 0]
    obs = obs[obs > 0]
    print("%s communities out of a total of %s communities were dropped because no species were observed in the given abundance range"
          % (num_dropped_communities, num_dropped_communities + len(obs)))
    macroecotools.plot_color_by_pt_dens(pred, obs, 3, loglog=1)
    
def sim_null_curry(tup):
    """Wrapping function to allow sim_null to work with multiprocessing"""
    return sim_null(*tup)
    
def sim_null(S0, N0, dic_beta):
    """Abundances simulated from a discrete uniform and associated METE predictions"""
    N_sim = sorted(np.random.random_integers(1, (2 * N0 - S0) / S0, S0), reverse = True)
    N_tot = sum(N_sim)
    
    #In cases where N and S are nearly equal it is possible for random draws to
    #yield all singletons which breaks the numerical solutions for Beta.
    #If this is the case make one species a doubleton.
    if N_tot == S0:
        N_sim[0] = 2
        
    if (S0, N0) not in dic_beta:
        dic_beta[(S0, N0)] = mete.get_beta(S0, sum(N_sim))
    N_pred = mete.get_mete_rad(S0, sum(N_sim), dic_beta[(S0, N0)])[0] 
    np.random.seed()
    return N_sim, N_pred

def create_null_dataset(Svals, Nvals, Niter, dataset_name, data_dir='./data/',
                        dic_filename='beta_lookup_table.pck', return_obs_pred=0):
    """Create simulated fits to uniform abundance distribution data
    
    Create list of coefficients of determination for simulated observed vs.
    predicted abundance relationships for a dataset. If these values are
    similar to those observed then the constraints alone largely determine the
    fit to the data. If they are weaker than the observed then the application
    of maximum entropy is also important.
    
    Svals : a list of values of observed species richnesses to match. Each row
            is a community (e.g., a site, year, etc.)
    Nvals : a list of values of observed community abundances to match. The
            ordering of rows should match that of Svals so that each row
            represents the S and N values for a single community.
    Niter: number of simulations
    dataset_name : short code that will indicate the name of the dataset in
                    the output file names
    data_dir : directory in which to store output    
    
    """
    resultfile = open(data_dir + dataset_name + '_sim_r2.csv', 'wb')
    out = csv.writer(resultfile, dialect = 'excel')
    dic_beta = mete.get_beta_dict(dic_filename)    
    for i in range(Niter):
        pool = multiprocessing.Pool()
        curried_args = itertools.izip(Svals, Nvals, itertools.repeat(dic_beta))
        site_sim_results = pool.map(sim_null_curry, curried_args)
        pool.close()
        
        sim_obs = []
        sim_pred = []
        for site in site_sim_results:
            sim_obs.extend((site[0]))
            sim_pred.extend((site[1]))
        r2 = macroecotools.obs_pred_rsquare(np.array(np.log10(sim_obs)),
                                       np.array(np.log10(sim_pred)))
        results = ((np.column_stack((i, r2))))
        out.writerows(results)
    resultfile.close()
    if return_obs_pred == 1:
        return sim_obs, sim_pred
    
def plot_avg_deviation_from_logseries(sites, obs_ab, p=None, sites_for_p=None,
                                      error_bars=0, color='b'):
    """Plot a figure showing deviations from the log-series as a function of ab
    
    Takes the obs-pred data for individual sites, groups them into Preston bins,
    stores the difference between observed and predicted data within each bin
    for each site, and then plots the average deviation against the center of
    the bin. Deviations are calculated as the percentage deviation within each
    bin, so if there is a difference of one species in a bin with 10 predicted
    species the deviation = 0.1.
    
    """
    usites = np.unique(sites)
    max_N = max(obs_ab)
    max_integer_logN = int(np.ceil(np.log2(max_N)) + 1)
    log_bin_edges = np.array(range(0, max_integer_logN))
    bin_edges = np.exp2(log_bin_edges)
    deviations = np.zeros((len(usites), len(bin_edges)-1))
    for i, site in enumerate(usites):
        site_abundances = obs_ab[sites == site]
        S = len(site_abundances)
        N = sum(site_abundances)
        obs_sad = macroecotools.preston_sad(site_abundances, b=bin_edges)
        if p==None:
            pred_sad = mete.get_mete_sad(S, N, bin_edges=bin_edges)
        else:
            beta = -log(p[sites_for_p==site])
            pred_sad = mete.get_mete_sad(S, N, beta=beta,
                                         bin_edges=bin_edges)
        deviation_from_predicted = (obs_sad[0] - pred_sad) / S * 100
        deviations[i,:] = deviation_from_predicted
    bin_numbers = range(1, max_integer_logN)
    mean_deviations = stats.nanmean(deviations)
    if error_bars == 1:
        std_deviations = stats.nanstd(deviations)
        plt.errorbar(bin_numbers, mean_deviations, yerr=std_deviations, fmt='b-')
    else:
        plt.plot(bin_numbers, mean_deviations, color=color, linewidth=3)
    plt.show()
    
def get_combined_obs_pred_data(datasets, data_dir='./data/'):
    """Combine obs-pred data from multiple datasets"""
    for i, dataset in enumerate(datasets):
        file_data = np.genfromtxt(data_dir + dataset + '_obs_pred.csv',
                                  dtype = "S15,i8,i8,i8",
                                  names = ['site','year','obs','pred'],
                                  delimiter = ",")
        #file_data = np.column_stack([i * np.ones(len(file_data)), file_data])
        if i == 0:
            data = file_data
        else:
            data = np.concatenate([data, file_data])
    return data

def run_empir_analysis(datasets, data_dir='./data/'):
    for dataset in datasets:
        raw_data = import_raw_data(workdir + dataset + '_spab.csv')
        run_test(raw_data, dataset, data_dir=workdir)

def run_sim_analysis(datasets, workdir, Niter):
    for dataset in datasets:
        obs_SN_data = import_obs_SN_data(workdir + dataset +
                                         '_dist_test.csv')
        create_null_dataset(obs_SN_data['Svals'], obs_SN_data['Nvals'],
                            Niter, dataset, data_dir=workdir)
    
def map_sites(datasets, data_dir='./data/', markers = ['o'],
              colors = ['b', 'r', 'g', 'y', 'c'], markersizes=3):
    """Generate a world map with sites color-coded by database"""
    
    from mpl_toolkits.basemap import Basemap
    
    map = Basemap(projection='merc',llcrnrlat=-57,urcrnrlat=71,\
                llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='l')
    map.drawcoastlines(linewidth = 1.25)

    for i, dataset in enumerate(datasets):
        latlong_data = import_latlong_data(data_dir + dataset + '_lat_long.csv')
        lats = latlong_data["lat"]
        longs = latlong_data["long"]
        x,y = map(longs,lats)
        map.plot(x,y, ls = '', marker = markers[i], markerfacecolor = colors[i], 
                 markeredgewidth = 0.25, markersize = markersizes)
    
    plt.savefig('map.png', dpi=400, facecolor='w', edgecolor='w', 
                bbox_inches = 'tight', pad_inches=0)
    
def map_sites_inset(datasets, data_dir='./data/', markers = ['o'],
                    colors=['b', 'r', 'g', 'y', 'c'], markersizes=4):
    """Generate a US map with sites color-coded by database"""
    
    from mpl_toolkits.basemap import Basemap
    
    map = Basemap(projection='merc',llcrnrlat=24.5,urcrnrlat=49,\
                llcrnrlon=-125,urcrnrlon=-69,lat_ts=50,resolution='l')
    
    map.drawcoastlines(linewidth=1.5)
    map.drawcountries(linewidth=1.5)
    map.drawmapboundary()     
    
    for i, dataset in enumerate(datasets):
        latlong_data = import_latlong_data(data_dir + dataset + '_lat_long.csv')
        lats = latlong_data["lat"]
        longs = latlong_data["long"]  
    
        x,y = map(longs,lats)
        map.plot(x,y, ls = '', marker = markers[i], markerfacecolor = colors[i], 
                 markeredgewidth = 0.25, markersize = markersizes)
    
    plt.savefig('map_inset.png', dpi=400, facecolor='w', edgecolor='w', 
                bbox_inches = 'tight', pad_inches=0)

def example_sad_plot(dataset, site_id, color, axis_limits, data_dir='./data/'):
    """Generate an example SAD plot for the map figure"""    
    obs_pred_data = import_obs_pred_data(data_dir + dataset + '_obs_pred.csv')    
    site = obs_pred_data["site"]
    obs = obs_pred_data["obs"]   
    pred = obs_pred_data["pred"]
    site_obs_ab = obs[site==site_id]
    site_pred_ab = pred[site==site_id]
    rank_obs, relab_obs = macroecotools.get_rad_data(site_obs_ab)
    rank_pred, relab_pred = macroecotools.get_rad_data(site_pred_ab)
    plt.figure(figsize=(2,2))
    plt.semilogy(rank_obs, relab_obs, 'o', markerfacecolor='none', markersize=8, 
             markeredgecolor=color, markeredgewidth=1.5)
    plt.semilogy(rank_pred, relab_pred, '-', color='black', linewidth=2)
    plt.axis(axis_limits)
    plt.xticks(fontsize = '10')
    plt.yticks(fontsize = '10')
    plt.savefig(dataset + '_example_' + str(site_id) + '.png', dpi=400,
                facecolor='w', edgecolor='w', bbox_inches = 'tight',
                pad_inches=0.2)
    
if __name__ == '__main__':
    assert len(sys.argv) >= 3, """You must provide at least two arguments:
    1. a path to the where the data is or will be stored
    2. argument for the type of analysis to be conducted"""
    workdir = sys.argv[1]
    
    #Determine which analyses to run
    analysis = sys.argv[2]
    if analysis not in ('empir', 'sim', 'figs', 'all'):
        print("The second argument should be empir, sim, figs, or all. See the docs")    
    if analysis == 'all':
        analyses = ('empir', 'sim', 'figs')
    else:
        analyses = (analysis)
        
    #Determine which datasets to use
    if os.path.exists(workdir + 'dataset_config.txt'):
        dataset_config_file = open(workdir + 'dataset_config.txt', 'r')
        datasets = []
        for line in dataset_config_file:
            datasets.append(line.strip())
    else:
        datasets = ['bbs', 'cbc', 'fia', 'gentry', 'mcdb', 'naba']  
    
    #Run selected analysss
    if 'empir' in analyses:
        run_empir_analysis(datasets, workdir)
    if 'sim' in analyses:
        assert len(sys.argv) == 4, """Running simulation analyses requires a
        third argument specifying the number of iterations"""
        iterations = int(sys.argv[3])
        run_sim_analysis(datasets, workdir, iterations)
    if 'figs' in analyses:
        colors = ['#87a4ef', '#0033cc', '#97ca82', '#339900','#ff6600', '#990000']
        
        #Figure 1 (Map & Example Fits)
        #Global Map
        map_sites(datasets, markers = ['o','o','s','s','D','v'], colors=colors,
                  markersizes=4)
        #US Map
        plt.figure()
        map_sites_inset(datasets, markers=['o','o','s','s','D','v'],
                        colors=colors, markersizes=5)
        #Example SADs
        site_ids = ['17220', 'L13428', '211131000022', '84', '1353', 'TX_Still_ollow']
        axis_limits = [[0, 16, 10 ** -4, 1], [-5, 115, 10 ** -5, 1],
                       [0, 15, 10 ** -2, 1], [-10, 225, 10 ** -3, 1], 
                       [0, 11, 10 ** -2, 1], [-2, 44, 10 ** -3, 1]]
        for i, dataset in enumerate(datasets):
            example_sad_plot(dataset, site_ids[i], colors[i], axis_limits[i])        

        #Figure 2 (Observed-predicted plots for each dataset)
        #Figure 3 ()
        #Figure 4 ()
        #Supplemental Figure X ()
    plt.show()    