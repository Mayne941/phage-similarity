---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.2
  kernelspec:
    display_name: Python 3.8.12 ('base')
    language: python
    name: python3
---



```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
import plotly.express as px
import networkx as nx 

le = LabelEncoder()
df = pd.read_csv("../MASH_dist_01Mar2022.tsv", delimiter="\t", engine="python")
df_no_dups = df[~(df["binA"] == df["binB"])]
df_labs = df_no_dups.copy()
```

```python
df_labs
```

```python
df_adj = pd.DataFrame(df_labs[["binA","binB","distance"]].sample(1000))
#df_adj = df_adj.set_index("binA")
```

```python
genome_names = df_adj["binA"].tolist()
```

```python
G = nx.Graph()
for i in genome_names:
    G.add_node(i)
```

```python
binB = df_adj["binB"].tolist()
dist = df_adj["distance"].tolist()
edges = []
for i in range(len(genome_names)):
    edges.append([(genome_names[i], binB[i], {"distance": dist[i]})])
```

```python
for i in edges:
    G.add_edges_from(i)
```

```python
import plotly.graph_objects as go
class DemoPlot3d:
    def __init__(self, G):#, colour_lut):
        self.k = 0.6 # Node separation
        self.G = G
        self.pos = nx.spring_layout(self.G, k = self.k, dim=3)
        #self.colour_lut = colour_lut

    def make_edge(self, x, y, z, text, width):
        return  go.Scatter3d(x         = x,
                        y         = y,
                        z         = z,
                        line      = dict(width = width,
                                    color = '#00bc9c'), # Osprey light gn
                        hoverinfo = 'text',
                        text      = ([text]),
                        mode      = 'lines')

    def make_edges(self):
        '''For each edge, make an edge_trace, append to list'''
        edge_trace = []
        for edge in self.G.edges():
            #if G.edges()[edge]['weight'] > 0:
            char_1 = edge[0]
            char_2 = edge[1]
            x0, y0, z0 = self.pos[char_1]
            x1, y1, z1 = self.pos[char_2]
            text = char_1 + '--' + char_2 + ': ' #+ str(G.edges()[edge]['weight'])
            trace  = self.make_edge([x0, x1, None], [y0, y1, None], [z0, z1, None], text, width = 0.3)  #*G.edges()[edge]['weight']**1.75)
            edge_trace.append(trace)
        return edge_trace

    def make_nodes(self):
        '''Make a node trace'''
        node_trace = go.Scatter3d(x         = [],
                                y         = [],
                                z         = [],
                                text      = [],
                                textposition = "middle center",
                                textfont_size = 10,
                                mode      = 'markers+text',
                                hoverinfo = 'text',
                                marker    = dict(color = [],
                                                size  = [],
                                                line  = None))
        # For each node in midsummer, get the position and size and add to the node_trace
        for node in self.G.nodes():
            x, y, z = self.pos[node]
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_trace['z'] += tuple([z])
            node_trace['marker']['color'] += tuple(['cornflowerblue'])
            #node_trace['marker']['size'] += tuple([5*self.G.nodes()[node]['degree']]) # was [5*G.nodes()[node]['size']]
            node_trace['text'] += tuple(['<b>' + node + '</b>'])
        return node_trace

    def make_plot(self, edge_trace, node_trace):
        axis=dict(showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title='',
            )
        layout = go.Layout(
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
                ),
            title = "lots of genomes",
            paper_bgcolor="#333333", # osprey grey bg
            plot_bgcolor="#333333", 
            font_color='#f7f7f7',
            #width=1200,
            #height=650,
            xaxis =  {'showgrid': False, 'zeroline': False}, # no gridlines
            yaxis = {'showgrid': False, 'zeroline': False}, 
            hovermode='closest',
        )
        fig = go.Figure(layout = layout)# Add all edge traces
        for trace in edge_trace:
            fig.add_trace(trace)# Add node trace
        trace2 = go.Scatter3d(x=node_trace["x"],
               y=node_trace["y"],
               z=node_trace["z"],
               mode='markers',
               name='actors',
               marker=dict(symbol='circle',
                             size=node_trace['marker']['size'],
                             #color=self.colour_lut,
                             colorscale='YlGnBu',
                             line=dict(color='rgb(50,50,50)', width=0.5)
                             ),
               text=node_trace["text"],
               hoverinfo='text'
               )
        fig.add_trace(trace2)
        fig.update_layout(showlegend = False)
        fig.update_xaxes(showticklabels = False)
        fig.update_yaxes(showticklabels = False)
        return fig

    def make_all(self):
        edge_trace = self.make_edges()
        node_trace = self.make_nodes()
        fig = self.make_plot(edge_trace, node_trace)
        fig.show()
        return edge_trace, node_trace

demo = DemoPlot3d(G)#, colour_by_type) # Initialise plot method, give it graph object and colour lookup table
edge_trace, node_trace = demo.make_all() # Call plotting function
```

```python
import plotly.graph_objects as go 
import networkx as nx

class DemoPlot2d_2:
    def __init__(self, G):
        self.k = 0.6 # Node separation
        self.G = G
        self.pos = nx.kamada_kawai_layout(self.G) #, k = self.k)

    def make_edge(self, x, y, text, width):
        return  go.Scatter(x         = x,
                        y         = y,
                        line      = dict(width = width,
                                    color = '#00bc9c'), # Osprey light gn
                        hoverinfo = 'text',
                        text      = ([text]),
                        mode      = 'lines')

    def make_edges(self):
        '''For each edge, make an edge_trace, append to list'''
        edge_trace = []
        for edge in self.G.edges():
            #if G.edges()[edge]['weight'] > 0:
            char_1 = edge[0]
            char_2 = edge[1]
            x0, y0 = self.pos[char_1]
            x1, y1 = self.pos[char_2]
            text = char_1 + '--' + char_2 + ': ' #+ str(G.edges()[edge]['weight'])
            trace  = self.make_edge([x0, x1, None], [y0, y1, None], text, width = 0.3)  #*G.edges()[edge]['weight']**1.75)
            edge_trace.append(trace)
        return edge_trace

    def make_nodes(self):
        '''Make a node trace'''
        node_trace = go.Scatter(x         = [],
                                y         = [],
                                text      = [],
                                textposition = "middle center",
                                textfont_size = 10,
                                mode      = 'markers',#+text',
                                hoverinfo = 'text',
                                marker    = dict(color = [],
                                                size  = [],
                                                line  = None))
        # For each node in midsummer, get the position and size and add to the node_trace
        for node in self.G.nodes():
            x, y = self.pos[node]
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_trace['marker']['color'] += tuple(['cornflowerblue'])
            #node_trace['marker']['size'] += tuple([5*self.G.nodes()[node]['degree']]) # was [5*G.nodes()[node]['size']]
            node_trace['text'] += tuple(['<b>' + node + '</b>'])
        return node_trace

    def make_plot(self, edge_trace, node_trace):
        axis=dict(showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title='',
            )
        layout = go.Layout(
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
                ),
            title = "2D Kamada-Kawai layout (path length cost function)",
            paper_bgcolor="#333333", # osprey grey bg
            plot_bgcolor="#333333", 
            font_color='#f7f7f7',
            width=1200,
            height=650,
            xaxis =  {'showgrid': False, 'zeroline': False}, # no gridlines
            yaxis = {'showgrid': False, 'zeroline': False}, 
            hovermode='closest',
        )
        fig = go.Figure(layout = layout)# Add all edge traces
        for trace in edge_trace:
            fig.add_trace(trace)# Add node trace
        fig.add_trace(node_trace)# Remove legend
        fig.update_layout(showlegend = False)
        fig.update_xaxes(showticklabels = False)
        fig.update_yaxes(showticklabels = False)
        fig.show()

    def make_all(self):
        edge_trace = self.make_edges()
        node_trace = self.make_nodes()
        self.make_plot(edge_trace, node_trace)

plot2 = DemoPlot2d_2(G)
plot2.make_all()
```
