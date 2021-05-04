--- 
layout: posts 
classes: wide 
sidebar:
  nav: "content" 
---
```
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', help = 'csv_dataframe')
parser.add_argument('sample_column', help = 'what is the name of the sample_column')
args = parser.parse_args()

#Reading in file
df = pd.read_csv(args.infile)

#Setting font and color pallette
#current_palette = sns.color_palette("Greys", 200)  #grey palette
#current_palette = sns.color_palette("ch:start=.2,rot=-.3", 200) # dark blues
current_palette = sns.color_palette("Blues", 200)  #grey palette
#current_palette = sns.diverging_palette(220, 20, as_cmap=True) #diverging bluered
#sns.set(font_scale = 0.3)

#set sample name and color
sample = df.pop(args.sample_column)


#Designing the plot
g = sns.clustermap(df, cmap = current_palette, robust = True, yticklabels = sample, xticklabels = 1)
plt.setp(g.ax_heatmap.xaxis.get_ticklabels(), fontsize = 15)
plt.setp(g.ax_heatmap.yaxis.get_ticklabels(), fontsize = 8)
plt.gcf().subplots_adjust(bottom=0.35)
fig = plt.gcf()
fig.set_size_inches((20, 16))


#Showing the plot
plt.savefig('depletolites_hc_lightblues.png', dpi = 600)
plt.show()```