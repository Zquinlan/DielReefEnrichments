import plotly.express as px
import pandas as pd


df = pd.read_csv('~/Documents/GitHub/DORCIERR/data/analysis/depletionSunburst.csv')

organisms = ['Dictyota', 'Turf', 'CCA', 'Pocillopora verrucosa', 'Porites lobata']
diel = ['Day', 'Night']

# for time in diel:
for organism in organisms:
    orgFilt = df[df['Organism'] == organism]
    orgDayNight = orgFilt.groupby(by="DayNight").sum()[["xic"]].reset_index()
    # dayOrgFilt = orgFilt[orgFilt['DayNight'] == time]

    fig = px.sunburst(orgFilt, 
        path = ['DayNight', 'superclass', 'class', 'subclass'], 
        values='xic', 
        color = 'DayNight',
        color_discrete_map = {'Day': 'FF9300', 'Night': 'grey'})
    
    name = str(organism + '_sunburst.png') 
    fig.write_image(name, width=1000, height=1000)        