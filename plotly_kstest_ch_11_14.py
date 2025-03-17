
import numpy as np
import pandas as pd

import plotly
import plotly.graph_objs
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from scipy import stats
from scipy.stats import ks_2samp

import matplotlib.pyplot as plt

BINS = 1000
z_score_max = 2
# px.colors.DEFAULT_PLOTLY_COLORS

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



def ks_test (wt, ko, x_label, xaxis_dtick): 
    
    D_auto, p_value = ks_2samp(wt, ko)
    
    all_values = np.concatenate((wt, ko))
    min_all_values = min(all_values)
    max_all_values = max(all_values)
    
    count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    
    count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
    pdf2 = count2 / sum(count2)
    cdf2 = np.cumsum(pdf2)
        
    D = abs(cdf-cdf2)
    maximum = max(D)
    
    x_maximum = bins_count[np.argmax(D)]
    
    y_min = cdf[np.argmax(D)]
    y_max = cdf2[np.argmax(D)]
    

    d = {'WT': cdf, 'KO': cdf2}
    df = pd.DataFrame(data=d)
    
    fig = px.line(df, x=bins_count[1:], y=['WT', 'KO'], template="simple_white", width=800, height=800, line_shape="vh")

    fig.add_trace(
        
        go.Scatter(
        x=[x_maximum, x_maximum],
        y=[y_min, y_max], 
        mode='lines', 
        line=dict(color='Black', dash='dash'), 
        showlegend=False))
    
    fig.update_layout(
        
        autosize= False, width=1000, height=1000, margin = dict(l=100, r=100, b=100, t=100, pad=4),    

        xaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = xaxis_dtick),         
        yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.2),  
        
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.2), 
        legend_title_text=None, 
        
        #title="Plot Title",
        xaxis_title=x_label,
        yaxis_title="CDF",

        font=dict(family="arial", size=36, color="Black"))
    
    plotly.offline.plot(fig, filename= "html_kstest_" +x_label + ".html")
    plt.figure() 
    
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

D_area_11, D_area_11_manual_calc, p_value_area_11 = ks_test (area_11_plus, area_11_minus, "Relative territory Chr. #11", 0.1)
D_area_14, D_area_14_manual_calc, p_value_area_14 = ks_test (area_14_plus, area_14_minus, "Relative territory Chr. #14", 0.1)
D_circ_11, D_circ_11_manual_calc, p_value_circ_11 = ks_test (circ_11_plus, circ_11_minus, "Circularity of Chr. #11", 0.5)
D_circ_14, D_circ_14_manual_calc, p_value_circ_14 = ks_test (circ_14_plus, circ_14_minus, "Circularity of Chr. #14", 0.5)
D_coloc_11_14, D_coloc_11_14_manual_calc, p_value_coloc_11_14 = ks_test (coloc_11_14_plus, coloc_11_14_minus, "Relative colocalization", 0.05)






def ks_test_subplot (all_data, all_xlabels): 
    
    fig = make_subplots(rows=1, cols=5, shared_yaxes=True) #,horizontal_spacing=0.05, vertical_spacing=0.1

    for i in range(5):

        wt = all_data[i][0]
        ko = all_data[i][1]

        D_auto, p_value = ks_2samp(wt, ko)
        
        all_values = np.concatenate((wt, ko))
        min_all_values = min(all_values)
        max_all_values = max(all_values)
        
        count, bins_count = np.histogram(wt, range = (min_all_values,max_all_values), bins=BINS)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        
        count2, bins_count2 = np.histogram(ko, range = (min_all_values,max_all_values), bins=BINS)
        pdf2 = count2 / sum(count2)
        cdf2 = np.cumsum(pdf2)

        D = abs(cdf-cdf2)
        maximum = max(D)
        
        x_maximum = bins_count[np.argmax(D)]
        
        y_min = cdf[np.argmax(D)]
        y_max = cdf2[np.argmax(D)]
  
        fig.add_trace(
            
            go.Scatter(
            x= bins_count[1:],
            y= cdf, 
            name= 'WT',
            line= dict(color='rgb(31, 119, 180)'),
            showlegend= i==0), 
            row=1, col=i+1)            
        
        fig.add_trace(
            
            go.Scatter(            
            x= bins_count[1:],
            y= cdf2, 
            name= 'KO',
            line= dict(color='rgb(255, 127, 14)'),
            showlegend= i==0), 
            row=1, col=i+1)           

        
        fig.add_trace(
            
            go.Scatter(
            x= [x_maximum, x_maximum],
            y= [y_min, y_max], 
            mode= 'lines', 
            line= dict(color='Black', dash='dash'), 
            showlegend= False), 
            row=1, col=i+1)            
        
        if i ==0:
            fig['layout']['xaxis']['title'] = all_xlabels[i]
        else:
            fig['layout']['xaxis'+str(i+1)]['title'] = all_xlabels[i]
        
    fig.update_layout(
        
        autosize=False, width=2000, height=565, margin = dict(l=100, r=100, b=100, t=100, pad=4),    
        
        xaxis= dict(mirror=True, ticks='outside', showline=True),        
        yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.5),         

        legend_title_text=None,
        #legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.2),
         
        #title="Plot Title",
        #xaxis_title="x_title",
        yaxis_title="CDF",
        
        font=dict(family="arial", size=36, color="Black"),
        
        template="simple_white")
    
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

    plotly.offline.plot(fig, filename= "html_kstest_subplot.html")


all_data = [(area_11_plus, area_11_minus),
            (area_14_plus, area_14_minus),
            (circ_11_plus, circ_11_minus),
            (circ_14_plus, circ_14_minus),
            (coloc_11_14_plus, coloc_11_14_minus)] 
             
    
all_xlabels = ["Relative territory<br>Chr. #11", 
               "Relative territory<br>Chr. #14", 
               "Circularity of<br>Chr. #11", 
               "Circularity of<br>Chr. #14", 
               "Relative<br>colocalization"]


ks_test_subplot (all_data, all_xlabels)

