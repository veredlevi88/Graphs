
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plotly
import plotly.graph_objs
import plotly.express as px
import plotly.graph_objects as go

from scipy import stats
from scipy.stats import ttest_ind
from plotly.subplots import make_subplots

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

area_all_Relative_territory = len(area_11_minus) + len(area_11_plus) + len(area_14_minus) + len(area_14_plus)

df_area = pd.DataFrame(np.zeros((area_all_Relative_territory,3)), columns = ["Chromosome number", 'Cell type', 'Relative territory'])

df_area.loc[:(len(area_11_minus) + len(area_11_plus) -1),"Chromosome number"] = 'Chr. #11'
df_area.loc[len(area_11_minus) + len(area_11_plus):,"Chromosome number"] = 'Chr. #14'

df_area.loc[:len(area_11_plus) -1 ,'Cell type'] = 'WT'
df_area.loc[len(area_11_plus):len(area_11_minus) + len(area_11_plus) -1,'Cell type'] = 'KO'
df_area.loc[len(area_11_minus) + len(area_11_plus):len(area_11_minus) + len(area_11_plus) + len(area_14_plus) -1 ,'Cell type'] = 'WT'
df_area.loc[len(area_11_minus) + len(area_11_plus) + len(area_14_plus):,'Cell type'] = 'KO'

df_area['Relative territory'] = area_11_plus.append(area_11_minus, ignore_index = True).append(area_14_plus, ignore_index = True).append(area_14_minus, ignore_index = True)



fig1 = px.violin(df_area, y="Relative territory", x="Chromosome number", color="Cell type", box=True, hover_data=df_area.columns, points=False, template="simple_white")        

fig1.update_layout(
    
    autosize=False, width=1000, height=1000, margin = dict(l=100, r=100, b=100, t=100, pad=4),    

    xaxis= dict(mirror=True, ticks='outside', showline=True),
    yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.1),

    legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95), 
    legend_title_text=None,        

    #title="Plot Title",
    #yaxis_title="Y Axis Title",
    xaxis_title=None,

    font=dict(family="arial", size=36, color="Black"))

plotly.offline.plot(fig1, filename= 'html_violins_area' + ".html")
plt.figure()

# statistics
# area_order_grouped = df_area.groupby(by=["Chromosome number", 'Cell type'])

# area_order_grouped.mean()['Relative territory']
# area_order_grouped.std()['Relative territory']

# area_order_grouped.median()['Relative territory']

# area_order_grouped.quantile(0.25)['Relative territory']
# area_order_grouped.quantile(0.75)['Relative territory']

# area_order_grouped.min()['Relative territory']
# area_order_grouped.max()['Relative territory']

# ttest_ind(area_11_minus, area_11_plus)
# ttest_ind(area_14_minus, area_14_plus)


"""--------------------circ--------------------"""

circ_11_minus = pd.read_csv('circ_11_minus.txt')
circ_11_plus = pd.read_csv('circ_11_plus.txt')
circ_14_minus = pd.read_csv('circ_14_minus.txt')
circ_14_plus = pd.read_csv('circ_14_plus.txt')

circ_11_minus = circ_11_minus[(np.abs(stats.zscore(circ_11_minus)) < z_score_max).all(axis=1)]
circ_11_plus = circ_11_plus[(np.abs(stats.zscore(circ_11_plus)) < z_score_max).all(axis=1)]
circ_14_minus = circ_14_minus[(np.abs(stats.zscore(circ_14_minus)) < z_score_max).all(axis=1)]
circ_14_plus = circ_14_plus[(np.abs(stats.zscore(circ_14_plus)) < z_score_max).all(axis=1)]

circ_all_Circularity = len(circ_11_minus) + len(circ_11_plus) + len(circ_14_minus) + len(circ_14_plus)

df_circ = pd.DataFrame(np.zeros((circ_all_Circularity,3)), columns = ["Chromosome number", 'Cell type', 'Circularity'])

df_circ.loc[:(len(circ_11_minus) + len(circ_11_plus) -1),"Chromosome number"] = 'Chr. #11'
df_circ.loc[len(circ_11_minus) + len(circ_11_plus):,"Chromosome number"] = 'Chr. #14'

df_circ.loc[:len(circ_11_plus) -1 ,'Cell type'] = 'WT'
df_circ.loc[len(circ_11_plus):len(circ_11_minus) + len(circ_11_plus) -1,'Cell type'] = 'KO'
df_circ.loc[len(circ_11_minus) + len(circ_11_plus):len(circ_11_minus) + len(circ_11_plus) + len(circ_14_plus) -1 ,'Cell type'] = 'WT'
df_circ.loc[len(circ_11_minus) + len(circ_11_plus) + len(circ_14_plus):,'Cell type'] = 'KO'

df_circ['Circularity'] = circ_11_plus.append(circ_11_minus, ignore_index = True).append(circ_14_plus, ignore_index = True).append(circ_14_minus, ignore_index = True)



fig2 = px.violin(df_circ, y="Circularity", x="Chromosome number", color="Cell type", box=True, hover_data=df_circ.columns, points=False, template="simple_white")        

fig2.update_layout(
    
    autosize=False, width=1000, height=1000, margin = dict(l=100, r=100, b=100, t=100, pad=4),    

    xaxis= dict(mirror=True, ticks='outside', showline=True),
    yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 1),

    legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95), 
    legend_title_text=None,        

    #title="Plot Title",
    #yaxis_title="Y Axis Title",
    xaxis_title=None,

    font=dict(family="arial", size=36, color="Black"))

plotly.offline.plot(fig2, filename= 'html_violins_circ' + ".html")
plt.figure()

# statistics
# circ_order_grouped = df_circ.groupby(by=["Chromosome number", 'Cell type'])

# circ_order_grouped.mean()['Circularity']
# circ_order_grouped.std()['Circularity']

# circ_order_grouped.median()['Circularity']

# circ_order_grouped.quantile(0.25)['Circularity']
# circ_order_grouped.quantile(0.75)['Circularity']

# circ_order_grouped.min()['Circularity']
# circ_order_grouped.max()['Circularity']

# ttest_ind(circ_11_minus, circ_11_plus)
# ttest_ind(circ_14_minus, circ_14_plus)


"""--------------------coloc--------------------"""

coloc_11_14_minus = pd.read_csv('coloc_11_14_minus.txt')
coloc_11_14_plus = pd.read_csv('coloc_11_14_plus.txt')

coloc_11_14_minus = coloc_11_14_minus[(np.abs(stats.zscore(coloc_11_14_minus)) < z_score_max).all(axis=1)]
coloc_11_14_plus = coloc_11_14_plus[(np.abs(stats.zscore(coloc_11_14_plus)) < z_score_max).all(axis=1)]

coloc_all_Relative_colocalization = len(coloc_11_14_minus) + len(coloc_11_14_plus) 
df_coloc = pd.DataFrame(np.zeros((coloc_all_Relative_colocalization,2)), columns = ['Cell type', 'Relative colocalization'])

df_coloc.loc[:(len(coloc_11_14_plus) -1),'Cell type'] = 'WT'
df_coloc.loc[len(coloc_11_14_plus):,'Cell type'] = 'KO'

df_coloc['Relative colocalization'] = coloc_11_14_plus.append(coloc_11_14_minus, ignore_index = True)



fig3 = px.violin(df_coloc, y="Relative colocalization", x="Cell type", color="Cell type", box=True, hover_data=df_coloc.columns, points=False, template="simple_white")

fig3.update_layout(
    
    autosize=False, width=1000, height=1000, margin = dict(l=100, r=100, b=100, t=100, pad=4),    

    xaxis= dict(mirror=True, ticks='outside', showline=True),
    yaxis= dict(mirror=True, ticks='outside', showline=True, tickmode = 'linear', tick0 = 0, dtick = 0.05),
    
    showlegend= False,
    #legend= dict(yanchor="top", y=0.95, xanchor="right", x=0.95), 
    #legend_title_text=None,        

    #title="Plot Title",
    #yaxis_title="Y Axis Title",
    #xaxis_title=None,

    font=dict(family="arial", size=36, color="Black"))



plotly.offline.plot(fig3, filename= 'html_violins_coloc' + ".html")
plt.figure()

# statistics
# coloc_order_grouped = df_coloc.groupby(by=['Cell type'])

# coloc_order_grouped.mean()['Relative colocalization']
# coloc_order_grouped.std()['Relative colocalization']

# coloc_order_grouped.median()['Relative colocalization']

# coloc_order_grouped.quantile(0.25)['Relative colocalization']
# coloc_order_grouped.quantile(0.75)['Relative colocalization']

# coloc_order_grouped.min()['Relative colocalization']
# coloc_order_grouped.max()['Relative colocalization']

# ttest_ind(coloc_11_14_minus, coloc_11_14_plus)
