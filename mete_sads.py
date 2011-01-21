"""Project-specific Code for Testing METE's SAD Predictions

Required input = Abundances per species per site for one sampling period
    
All data queries used can be found in MaxEnt/trunk/data:
    BBS_data_query
    CBC_data_query
    Gentry_data_query
        
"""
    
import macroeco_distributions as md
import mete
import csv
from macroeco import plot_bivar_color_by_pt_density_relation as densityplt
import matplotlib.pyplot as plt
import numpy as np
import weestats

def run_test(input_filename, output_filename1, output_filename2, cutoff = 9):
    """Use data to compare the predicted and empirical SADs and get results in csv files
    
    Keyword arguments:
    input_filename -- path to file that has raw data in the format: 
                        'site','year','sp','ab'
    output_filename1 -- file that will store the pred and observed species-level 
                        abundances for each site in the input file
    output_filename2 -- file that will store the p-values and weights from 
                        dist_test for each site in the input file
    cutoff      --  minimum number of species required to run.
    
    """
    
    ifile = np.genfromtxt(input_filename, dtype = "S9,i8,S9,i8", 
                       names = ['site','year','sp','ab'], delimiter = ",")
    
    usites = np.sort(list(set(ifile["site"])))
    
    f1 = csv.writer(open(output_filename1,'a'))
    f2 = csv.writer(open(output_filename2,'a'))
    
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
    
    ifile = np.genfromtxt(input_filename, dtype = "S9,i8,i8", 
                       names = ['site','obs','pred'], delimiter = ",")
    
    pred = ((ifile["pred"]))
    obs = ((ifile["obs"]))
    
    plt.figure()
    densityplt(pred, obs, 1, loglog=1)
    plt.title(title)
    plt.xlabel('Predicted abundances')
    plt.ylabel('Observed abundances')
    plt.show()    
    
def plot_weights(input_filename, data = 'raw', left = [0, 0.4, 0.8], 
                 width = 0.2, color = 'b', title = ''): 
    """use output from run_test to plot frequency distribution of Akaike weights"""
    
    ifile = np.genfromtxt(input_filename, dtype = "S9,i8,i8,f8,f8", 
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
    colors = ['k', 'r', 'g', 'b']
    
    for i in range(0, n):
        input_filename = input_filenames[i]
        width = round(1.0/(3 + n * 3), 2)
        left = [(width * (i + 1)), (width * (i + n + 2)), (width * (i + n + 6))]
        plot_weights(input_filename, data = 'percent', left = left, 
                     color = colors[i], width = width)
    
    plt.ylabel('Percentage of sites')
    # TO DO: figure out universal means of determining xtick locations
    plt.xticks(((((3 + n * 3)/4.8) * width), (((3 + n * 3)/1.85) * width), 
               (((3 + n * 3)/1.14) * width)), 
               ('Log-normal', 'Indeterminate', 'Log-series') )
    # TO DO: figure out how to include a color-coded legend: 
    plt.legend(('CBC', 'BBS', 'Gentry'), loc = 'upper left')
    plt.show()