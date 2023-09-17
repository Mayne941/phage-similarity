import pandas as pd
import igraph as ig
import leidenalg as la
import community as louvain
import numpy as np

def do_louvain_partitions(gr, save_partition_data, fpath) -> pd.DataFrame():
    '''Use Louvain method to get partitions, save as JSON'''
    df = pd.DataFrame([louvain.best_partition(gr,weight='distance')]).T
    if save_partition_data == True:
        df.to_json(f"{fpath}full_louvain_partition.json")
    return df

def do_leiden_partitions(gr, save_partition_data, fpath) -> pd.DataFrame():
    '''Use Leiden method to get partitions, save as JSON'''
    #Convert networkx to igraph
    g = ig.Graph.from_networkx(gr)
    nodes_list = np.asarray(g.vs['_nx_name'])
    partitions = dict()
    for partition, nodes in enumerate(la.find_partition(g, partition_type=la.ModularityVertexPartition, seed=0)):
        mapping = dict(zip(nodes_list[nodes], [partition]*len(nodes)))
        partitions.update(mapping)
    df = pd.DataFrame(partitions.items())
    if save_partition_data == True:
        df.to_csv(f"{fpath}full_leiden_partition.json")
    return df
