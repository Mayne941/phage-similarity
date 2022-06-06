import plotly
import plotly.graph_objects as go


class ClusterPlot:
    def __init__(self, G, par, pos, sample_size, fpath):
        self.k = 0.6  # Node separation. Not used currently
        self.G = G
        self.pos = pos
        self.partition = par
        self.sample_size = sample_size
        self.fpath = fpath

    def make_edge(self, x, y, text, width):
        return go.Scatter(x=x,
                          y=y,
                          line=dict(width=width,
                                    color='#D3D3D3'),
                          hoverinfo='text',
                          text=([text]),
                          mode='lines')

    def make_edges(self):
        '''For each edge, make an edge_trace, append to list'''
        edge_trace = []
        for edge in self.G.edges():
            char_1 = edge[0]
            char_2 = edge[1]
            x0, y0 = self.pos[char_1]
            x1, y1 = self.pos[char_2]
            text = char_1 + '--' + char_2 + ': '
            trace = self.make_edge(
                [x0, x1, None], [y0, y1, None], text, width=1.0)  
            edge_trace.append(trace)
        return edge_trace

    def make_nodes(self):
        '''Make a node trace'''
        node_trace = go.Scatter(x=[],
                                y=[],
                                text=[],
                                textposition="middle center",
                                textfont_size=10,
                                mode='markers',  
                                hoverinfo='text',
                                marker=dict(color=[],
                                             size=5,
                                             line=None),
                                marker_colorscale=plotly.colors.sequential.Viridis)
        # For each node, get the position and size and add to the node_trace
        for node in self.G.nodes():
            x, y = self.pos[node]
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_trace['marker']['color'] += tuple([self.partition[node]])
            node_trace['text'] += tuple(['<b>' + node + '</b>\n' +
                                        "cluster: " + str(self.partition[node])])
        return node_trace

    def make_plot(self, edge_trace, node_trace):
        axis = dict(showbackground=False,
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
            title=f"2D projection, n={self.sample_size}",
            paper_bgcolor="#333333",  
            plot_bgcolor="#333333",
            font_color='#f7f7f7',
            width=1200,
            height=650,
            xaxis={'showgrid': False, 'zeroline': False},  # no gridlines
            yaxis={'showgrid': False, 'zeroline': False},
            hovermode='closest',
        )
        fig = go.Figure(layout=layout)  # Add all edge traces
        for trace in edge_trace:
            fig.add_trace(trace)  # Add node trace
        fig.add_trace(node_trace)  # Remove legend
        fig.update_layout(showlegend=False)
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        fig.write_image(f"{self.fpath}fig1.png")
        fig.show()

    def make_all(self):
        edge_trace = self.make_edges()
        node_trace = self.make_nodes()
        self.make_plot(edge_trace, node_trace)
