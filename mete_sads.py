"""Project-specific Code for Testing METE's SAD Predictions

Command line usage:
python mete_sads.py /path/to/data/ args

args:
empir - generate empirical fits, obs-pred data, and compare to log-normal
sims - generate simulated datasets and compare to empirical results

The main data required is abundances of each species at each site for one 
sampling period
    
All data queries used can be found in MaxEnt/trunk/data:
    BBS_data_query
    CBC_data_query
    Gentry_data_query
           
"""

#TODO 1. integrate lookup dictionary more broadly so that rerunning the analysis
#        avoids the time consuming numerical solutions.
#     2. update notation terminology to match Harte's 2011 book

from __future__ import division
import macroeco_distributions as md
import mete
import csv
import macroeco
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import weestats
import cPickle
import sys
import multiprocessing
import itertools
import os
from math import log, exp
import itertools

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
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,S10,i8", 
                       names = ['site','year','sp','ab'], delimiter = ",")
    
    usites = np.sort(list(set(ifile["site"])))
    
    f1 = csv.writer(open(output_filename1,'wb'))
    f2 = csv.writer(open(output_filename2,'wb'))
    
    for i in range(0, len(usites)):
        subsites = ifile["site"][ifile["site"] == usites[i]]        
        subab = ifile["ab"][ifile["site"] == usites[i]]
        subyr = ifile["year"][ifile["site"] == usites[i]]
        uyr = np.sort(list(set(subyr)))
        for a in range(0, len(uyr)):
            subsubsites = subsites[subyr == uyr[a]]
            subsubab = subab[subyr == uyr[a]]
            subsubyr = subyr[subyr == uyr[a]]
            N = sum(subsubab)
            S = len(subsubsites)
            if S > cutoff:
                # Generate predicted values and p (e ** -lambda_sad) based on METE:
                mete_pred = mete.get_mete_rad(int(S), int(N))
                pred = np.array(mete_pred[0])
                p = mete_pred[1]
                p_untruncated = exp(-mete.get_lambda_sad(S, N, version='untruncated'))
                subab3 = np.sort(subsubab)[::-1]
                # Calculate Akaike weight of log-series:
                L_logser = md.logser_ll(subab3, p)
                L_logser_untruncated = md.logser_ll(subab3, p_untruncated)
                mu, sigma = md.pln_solver(subab3)
                L_pln = md.pln_ll(mu,sigma,subab3)        
                k1 = 1
                k2 = 2    
                AICc_logser = weestats.AICc(k1, L_logser, S)
                AICc_logser_untruncated = weestats.AICc(k1, L_logser_untruncated, S)
                AICc_pln = weestats.AICc(k2, L_pln, S)
                weight = weestats.aic_weight(AICc_logser, AICc_pln, S, cutoff = 4)
                weight_untruncated = weestats.aic_weight(AICc_logser_untruncated,
                                                         AICc_pln, S, cutoff = 4)
                #save results to a csv file:
                results = ((np.column_stack((subsubsites, subsubyr, subab3, pred))))
                results2 = ((np.column_stack((np.array(usites[i], dtype='S20'),
                                                       uyr[a], S, N, p, weight,
                                                       p_untruncated,
                                                       weight_untruncated))))
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
    
    ifile = np.genfromtxt(input_filename, dtype = "S15,i8,i8,i8,f8,f8", 
                       names = ['site', 'year', 'S', 'N', 'p', 'weight'],
                       delimiter = ",")
    weights = ((ifile["weight"]))
    weights = weights[weights >= 0]
    bins = [0, 0.33333, 0.66667, 1]
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
    
def hist_mete_r2(sites, obs, pred):
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in sites:
        obs_site = obs[sites==site]
        pred_site = pred[sites==site]
        r2 = macroeco.obs_pred_rsquare(obs_site, pred_site)
        r2s.append(r2)

    #density_estimate = stats.kde.gaussian_kde(r2s)
    #xvals = np.arange(0, 1, 0.01)
    #yvals = density_estimate.evaluate(xvals)
    hist_r2 = np.histogram(r2s, range=(0, 1))
    xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
    xvals = xvals[0:len(xvals)-1]
    yvals = hist_r2[0]
    plt.plot(xvals, yvals, 'k-', linewidth=2)
    plt.axis([0, 1, 0, 1.1 * max(yvals)])
    
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

def plot_numsp_obs_pred (sites, obs_ab, min_abundance, max_abundance):
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
    macroeco.plot_color_by_pt_dens(pred, obs, 3, loglog=1)
    
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
        
        f1 = csv.writer(open(output_filenames[i],'ab'))
    
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
    titles = ('BBS', 'CBC','FIA','Gentry','MCDB','NABC')
    
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
            axis_min = 0 #min(obs) - 1
            axis_max = 1 #max(obs) + 1
            axis_scale = 0
        #slope, intercept, r_value, p_value, std_err = stats.linregress(pred_trans,
        #                                                               obs_trans)
        #r_squared = macroeco.obs_pred_rsquare(obs_trans, pred_trans)
            
        plt.subplot(3,2,i+1)
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
            plt.xlabel('Predicted Evar')
        plt.title(titles[i])
        #r2 = ('r2 = ' + str(round(r_squared, 2)))
        #b = ('y = ' + str(round(slope, 2)) + 'x + ' + str(round(intercept)))
        #plt.annotate(b, xy=(-10, 10), xycoords='axes points',
                #horizontalalignment='right', verticalalignment='bottom',
                #fontsize=14)
        #plt.annotate(r2, xy=(-10, 30), xycoords='axes points',
                #horizontalalignment='right', verticalalignment='bottom',
                #fontsize=14)
        
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
    
def sim_null_curry(tup):
    """Wrapping function to allow sim_null to work with multiprocessing"""
    return sim_null(*tup)
    
def sim_null(S0, N0, dic_lambda):
    """Abundances simulated from a discrete uniform and associated METE predictions"""
    N_sim = sorted(np.random.random_integers(1, (2 * N0 - S0) / S0, S0), reverse = True)
    N_tot = sum(N_sim)
    NS_ratio = N_tot / S0
    
    #In cases where N and S are nearly equal it is possible for random draws to
    #yield all singletons which breaks the numerical solutions. If this is the
    #case make one species a doubleton
    if N_tot == S0:
        N_sim[0] = 2
        
    if NS_ratio not in dic_lambda:
        dic_lambda[NS_ratio] = mete.get_lambda_sad(S0, sum(N_sim))
    N_pred = mete.get_mete_rad(S0, sum(N_sim), dic_lambda[NS_ratio])[0] 
    np.random.seed()
    return N_sim, N_pred

def create_null_dataset(input_filename, output_filename, Niter,
                        dic_filename='lambda_library.pck', return_obs_pred=0):
    """Create list of R^2 values for simulated observed vs. predicted 
    abundance relationships for a dataset.
    
    Niter: number of simulations
    input_filename: same format as output_filename2 from run_test
    output_filename: 1 column - R^2 for simulated data, one value per iteration
    
    """
    #TODO: We are currently not updating the lookup table. This is not
    #      trivial to do in parallel. Here is a starting point:
    #http://stackoverflow.com/questions/2545961/how-to-synchronize-a-python-dict-with-multiprocessing
    dist_test_results = np.genfromtxt(input_filename, usecols=(2, 3),
                                      dtype="i8,i8",
                                      names=['Svals','Nvals'],
                                      delimiter=",")
    resultfile = open(output_filename, 'wb')
    out = csv.writer(resultfile, dialect = 'excel')
    dic_lambda = mete.get_lambda_dict(dic_filename)    
    for i in range(Niter):
        pool = multiprocessing.Pool()
        curried_args = itertools.izip(dist_test_results['Svals'],
                                      dist_test_results['Nvals'],
                                      itertools.repeat(dic_lambda))
        site_sim_results = pool.map(sim_null_curry, curried_args)
        pool.close()
        
        sim_obs = []
        sim_pred = []
        for site in site_sim_results:
            sim_obs.extend((site[0]))
            sim_pred.extend((site[1]))
        r2 = macroeco.obs_pred_rsquare(np.array(np.log10(sim_obs)),
                                       np.array(np.log10(sim_pred)))
        results = ((np.column_stack((i, r2))))
        out.writerows(results)
    resultfile.close()
    if return_obs_pred == 1:
        return sim_obs, sim_pred
    
def plot_sads(sites, obs_ab, pred_ab, num_sites = 25):
    """Plot Preston SADs for both observed and predicted SADs for a subset of sites (num_sites)"""
    
    usites = sorted(list(set(sites)))
    step = len(usites)/(num_sites)
    a = 1
    for i in range (0, len(usites), step):
        if a <= num_sites:
            iobs = obs_ab[sites == usites[i]]
            ipred = pred_ab[sites == usites[i]]
            sad_obs = macroeco.preston_sad(iobs, normalized = 'yes')
            sad_pred = macroeco.preston_sad(ipred, normalized = 'yes')
            
            plot_obj = plt.subplot(5,5,a)
            height_obs = sad_obs[0]
            height_pred = sad_pred[0]
            left_obs = range(1,len(height_obs)+ 1)  
            left_pred = range(1,len(height_pred)+ 1)
            plt.bar(left_obs, height_obs, width = 0.4, color = 'b')
            plt.bar((np.array(left_pred) + 0.4), height_pred, width = 0.4, color = 'r')
            
            plot_obj.set_title(usites[i])
            plot_obj.set_ylabel('Number of species') 
            plot_obj.set_xlim(0.5, max(max(left_obs),max(left_pred)) + 1)
            plot_obj.set_xticks(np.array(left_obs) + 0.4)
            plot_obj.set_xticklabels(np.array(sad_obs[1][0:len(height_obs)], dtype = int))
            plot_obj.set_ylim(0, (max(max(height_obs),max(height_pred)) + 1))
            
            a += 1
    plt.show()
    
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
        obs_sad = macroeco.preston_sad(site_abundances, b=bin_edges)
        if p==None:
            pred_sad = mete.get_mete_sad(S, N, bin_edges=bin_edges)
        else:
            lambda_sad = -log(p[sites_for_p==site])
            pred_sad = mete.get_mete_sad(S, N, lambda_sad=lambda_sad,
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
        
def plot_sad_fit(sites, obs_ab, pred_ab, sites2, pr, dist = 'pln', 
                 normalized = 'yes', num_sites = 25):
    """Plot pmfs of fits of observed abundances"""
    
    usites = sorted(list(set(sites)))
    step = len(usites)/(num_sites)
    a = 1
    for i in range (1, len(usites), step):
        if a <= num_sites:
            iobs = obs_ab[sites == usites[i]]
            ipred = pred_ab[sites == usites[i]]
            sad_obs = macroeco.preston_sad(iobs, normalized = 'yes')
            sad_pred = macroeco.preston_sad(ipred, normalized = normalized)
            
            ipr = pr[sites2 == usites[i]]
            
            plot_obj = plt.subplot(5,5,a)
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, 
                                wspace=0.55, hspace=0.33)
            height_obs = sad_obs[0]
            height_pred = sad_pred[0]
            left_obs = range(1,len(height_obs)+ 1)  
            left_pred = range(1,len(height_pred)+ 1)
            plt.bar(left_obs, height_obs, width = 0.3, color = 'b')
            plt.bar((np.array(left_pred) + 0.3), height_pred, width = 0.3, color = 'r')
            
            if dist == 'logser':                
                L = stats.logser.pmf(range(1, max(iobs)+1),ipr)
            if dist == 'pln':
                mu = np.mean(np.log(iobs))
                sigma = np.std(np.log(iobs))
                L = md.pln_lik(mu,sigma,np.sort(iobs),approx_cut = 10)  
                L = md.pln_lik(mu,sigma,range(1, max(iobs)+1),approx_cut = 10)
            L_scaled = L * len(iobs)
            q = np.exp2(list(range(0, 25)))    
            b = q [(q <= max(iobs)*2)]
            height_dist=[]
            for c in range(0,len(b)-1):
                count = float(sum(L_scaled[b[c]-1:b[c+1]-1]))
                binwidth = b[c+1]-b[c]
                height_dist.append(count/binwidth)                
            left_dist = range(1,len(height_dist)+ 1)            
            plt.bar((np.array(left_dist) + 0.6), height_dist, width = 0.3, color = 'k')
                        
            plot_obj.set_title(usites[i]) 
            plot_obj.set_xlim(0.5, max(max(left_obs),max(left_pred), max(left_dist)) + 1)
            plot_obj.set_xticks(np.array(left_obs) + 0.4)
            plot_obj.set_xticklabels(np.array(sad_obs[1][0:len(height_obs)], dtype = int))
            plot_obj.set_ylim(0, (max(max(height_obs),max(height_pred), max(height_dist)) + 0.5))
            if a == 1 or a == 6 or a == 11 or a == 16 or a == 21:
                plot_obj.set_ylabel('Number of species')
            if max(max(height_obs), max(height_pred), max(height_dist)) > 59:
                tk = 20
            elif max(max(height_obs), max(height_pred), max(height_dist)) > 19:
                tk = 5
            elif max(max(height_obs), max(height_pred), max(height_dist)) > 10:
                tk = 2
            else:
                tk = 1
            plot_obj.set_yticks(range(tk,max(max(height_obs),max(height_pred)) + 1,tk))     
            print a
            a += 1
    plt.show()
    
def get_combined_obs_pred_data(inputfilenames):
    """Combine all obs-pred data from a list of run_test files"""
    for i, filename in enumerate(inputfilenames):
        file_data = np.genfromtxt(filename, dtype = "S15,i8,i8,i8",
                                  names = ['site','year','obs','pred'],
                                  delimiter = ",")
        #file_data = np.column_stack([i * np.ones(len(file_data)), file_data])
        if i == 0:
            data = file_data
        else:
            data = np.concatenate([data, file_data])
    return data

if __name__ == '__main__':
    assert len(sys.argv) >= 3, """You must provide at least two arguments,
    a path to the where the data is or will be stored, and an argument for the
    type of analysis to be conducted"""
    workdir = sys.argv[1]
    if os.path.exists(workdir + 'dataset_config.txt'):
        dataset_config_file = open(workdir + 'dataset_config.txt', 'r')
        datasets = []
        for line in dataset_config_file:
            datasets.append(line.strip())
    else:
        datasets = ['bbs', 'cbc', 'fia', 'gentry', 'mcdb', 'naba']
    if sys.argv[2] == 'empir':
        for dataset in datasets:
            run_test(workdir + dataset + '_spab.csv',
                     workdir + dataset + '_obs_pred.csv',
                     workdir + dataset + '_dist_test.csv')
    elif sys.argv[2] == 'sim':
        if len(sys.argv) == 4:
            Niter = int(sys.argv[3])
        else:
            Niter = 10
        for dataset in datasets:
            create_null_dataset(workdir + dataset + '_dist_test.csv',
                                workdir + dataset + '_sim_r2.csv', Niter)
    else:
        print "The second argument should be either empir or sim. See the docs"