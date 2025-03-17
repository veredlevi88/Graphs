
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
from scipy.stats import ks_2samp

sns.set(rc={"figure.dpi":300, 'savefig.dpi':300, "figure.figsize":(4,4)})
sns.set_context('notebook')
sns.set_style("ticks")

BINS = 1000
z_score_max = 2


"""--------------------area--------------------"""

area_11_minus = pd.read_csv('area_11_minus.txt')
area_11_plus = pd.read_csv('area_11_plus.txt')
area_14_minus = pd.read_csv('area_14_minus.txt')
area_14_plus = pd.read_csv('area_14_plus.txt')

area_11_minus = area_11_minus[(np.abs(stats.zscore(area_11_minus)) < z_score_max).all(axis=1)]
area_11_plus = area_11_plus[(np.abs(stats.zscore(area_11_plus)) < z_score_max).all(axis=1)]
area_14_minus = area_14_minus[(np.abs(stats.zscore(area_14_minus)) < z_score_max).all(axis=1)]
area_14_plus = area_14_plus[(np.abs(stats.zscore(area_14_plus)) < z_score_max).all(axis=1)]


"""--------------------circ--------------------"""

circ_11_minus = pd.read_csv('circ_11_minus.txt')
circ_11_plus = pd.read_csv('circ_11_plus.txt')
circ_14_minus = pd.read_csv('circ_14_minus.txt')
circ_14_plus = pd.read_csv('circ_14_plus.txt')

circ_11_minus = circ_11_minus[(np.abs(stats.zscore(circ_11_minus)) < z_score_max).all(axis=1)]
circ_11_plus = circ_11_plus[(np.abs(stats.zscore(circ_11_plus)) < z_score_max).all(axis=1)]
circ_14_minus = circ_14_minus[(np.abs(stats.zscore(circ_14_minus)) < z_score_max).all(axis=1)]
circ_14_plus = circ_14_plus[(np.abs(stats.zscore(circ_14_plus)) < z_score_max).all(axis=1)]


"""--------------------coloc--------------------"""

coloc_11_14_minus = pd.read_csv('coloc_11_14_minus.txt')
coloc_11_14_plus = pd.read_csv('coloc_11_14_plus.txt')

coloc_11_14_minus = coloc_11_14_minus[(np.abs(stats.zscore(coloc_11_14_minus)) < z_score_max).all(axis=1)]
coloc_11_14_plus = coloc_11_14_plus[(np.abs(stats.zscore(coloc_11_14_plus)) < z_score_max).all(axis=1)]



def ks_test (wt, ko, x_label): 
    
    D_auto, p_value = ks_2samp(wt, ko)
    
    all_values = np.concatenate((wt, ko))
    min_all_values = min(all_values)
    max_all_values = max(all_values)
    
    count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    plt.step(bins_count[1:], cdf, label="WT")
    
    count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
    pdf2 = count2 / sum(count2)
    cdf2 = np.cumsum(pdf2)
    plt.step(bins_count2[1:], cdf2, label="KO") # bins_count2 is the same as bins_count
        
    D = abs(cdf-cdf2)
    maximum = max(D)
    
    x_maximum = bins_count[np.argmax(D)]
    
    y_min = cdf[np.argmax(D)]
    y_max = cdf2[np.argmax(D)]
    
    plt.plot((x_maximum, x_maximum), (y_min, y_max), 'k' , linestyle="--")
    plt.xlabel(x_label)
    plt.ylabel("CDF")

    plt.legend(frameon=False)

    plt.show()
    
    return D_auto, maximum , p_value


area_11_plus = area_11_plus.to_numpy().flatten()
area_11_minus = area_11_minus.to_numpy().flatten()

area_14_plus = area_14_plus.to_numpy().flatten()
area_14_minus = area_14_minus.to_numpy().flatten()

circ_11_plus = circ_11_plus.to_numpy().flatten()
circ_11_minus = circ_11_minus.to_numpy().flatten()

circ_14_plus = circ_14_plus.to_numpy().flatten()
circ_14_minus = circ_14_minus.to_numpy().flatten()

coloc_11_14_plus = coloc_11_14_plus.to_numpy().flatten()
coloc_11_14_minus = coloc_11_14_minus.to_numpy().flatten()

D_area_11, D_area_11_manual_calc, p_value_area_11 = ks_test (area_11_plus, area_11_minus, "Relative territory Chr. #11")
D_area_14, D_area_14_manual_calc, p_value_area_14 = ks_test (area_14_plus, area_14_minus, "Relative territory Chr. #14")
D_circ_11, D_circ_11_manual_calc, p_value_circ_11 = ks_test (circ_11_plus, circ_11_minus, "Circularity of Chr. #11")
D_circ_14, D_circ_14_manual_calc, p_value_circ_14 = ks_test (circ_14_plus, circ_14_minus, "Circularity of Chr. #14")
D_coloc_11_14, D_coloc_11_14_manual_calc, p_value_coloc_11_14 = ks_test (coloc_11_14_plus, coloc_11_14_minus, "Relative colocalization")



def ks_test_subplot (all_data, all_xlabels): 
    
    fig, axs = plt.subplots(nrows = 1, ncols = 5, sharey=True, figsize=(9, 2))
    fig.tight_layout()    
    sns.set_theme(style="ticks",font="arial") #,fontsize=28)
   
    for i in range(5):
          
        wt = all_data[i][0]
        ko = all_data[i][1]
        x_label = all_xlabels[i]
        
        D_auto, p_value = ks_2samp(wt, ko)
        
        all_values = np.concatenate((wt, ko))
        min_all_values = min(all_values)
        max_all_values = max(all_values)
        
        count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        axs[i].plot(bins_count[1:], cdf, label="WT")
        
        count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
        pdf2 = count2 / sum(count2)
        cdf2 = np.cumsum(pdf2)
        axs[i].plot(bins_count2[1:], cdf2, label="KO") # bins_count2 is the same as bins_count  
        
        D = abs(cdf-cdf2)
        maximum = max(D)
        
        x_maximum = bins_count[np.argmax(D)]
        
        y_min = cdf[np.argmax(D)]
        y_max = cdf2[np.argmax(D)]
            
        axs[i].plot((x_maximum, x_maximum), (y_min, y_max), 'k' , linestyle="--")     
        axs[i].set_xlabel(x_label)
   
        if i == 0:
            axs[i].set_ylabel('CDF')
    
        plt.legend(frameon=False)
    

sns.set()
sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_context('notebook')
sns.set_style("ticks")


all_data = [(area_11_plus, area_11_minus),
            (area_14_plus, area_14_minus),
            (circ_11_plus, circ_11_minus),
            (circ_14_plus, circ_14_minus),
            (coloc_11_14_plus, coloc_11_14_minus)] 
             
    
all_xlabels = ["Relative territory\nChr. #11", 
               "Relative territory\nChr. #14", 
               "Circularity of\nChr. #11", 
               "Circularity of\nChr. #14", 
               "Relative\ncolocalization"]


ks_test_subplot (all_data, all_xlabels)
