import pandas as pd
import os
import community as louvain
import numpy as np
from sklearn.preprocessing import LabelEncoder
import plotly.express as px
import networkx as nx
import time
from datetime import datetime
from .cluster_graph import ClusterPlot
from .dendrogram import create_dendrogram

class GetPartitions:
    def __init__(self, params) -> None:
        '''Set up exp folder, decode runtime params'''
        self.le = LabelEncoder()
        self.sample_size = params["sample_size"]
        self.distance_value = params["distance_value"]    
        self.input_file_path = params["input_file_path"]
        self.taxonomy_file_path = params["taxonomy_file_path"]
        self.save_partition_data = params["save_partition_data"]
        self.plot_graph = params["plot_graph"]
        self.dendro_threshold = params["dendro_threshold"]
        self.trunc = params["dendro_truncation"]
        self.dendro_width = params["dendro_width"]
        self.hash_threshold = params["hash_threshold"]
        self.p_threshold = params["p_threshold"]
        if not os.path.exists("./output/"):
            os.mkdir("./output/")
        self.fpath = f"./output/exp_{datetime.now().strftime('%d-%b-%Y--%H-%M-%S')}/"
        if not os.path.exists(self.fpath):
            os.mkdir(self.fpath)

    def load_and_clean(self) -> pd.DataFrame:
        '''Read in data, filter duplicates and by p value threshold'''
        df = pd.read_csv(self.input_file_path, delimiter="\t", engine="python", header=0, names=["binA", "binB", "distance", "p-value", "tX"]) 
        df_no_dups = df[~(df["binA"] == df["binB"])]
        '''Hash filter'''
        df_no_dups["t"] = df_no_dups["tX"].str.split("/").str[0].astype(int)
        df_no_dups["X"] = df_no_dups["tX"].str.split("/").str[1].astype(int)
        df_no_dups["hashes"] = df_no_dups["t"]/df_no_dups["X"]
        print(f"Filtering by hashes has removed {df_no_dups.shape[0] - df_filtered.shape[0]} samples.")
        '''Sample size filter'''
        if self.sample_size <= 0:
            self.sample_size = df_filtered.shape[0]
        return df_filtered.head(self.sample_size)

    def make_graph(self, df) -> nx.Graph():
        '''Make NetworkX graph object'''
        genome_names = df["binA"].tolist()  
        binB = df["binB"].tolist()
        dist = df["distance"].tolist()
        gr = nx.Graph()
        [gr.add_node(i) for i in df["binA"].unique().tolist()]
        edges = [[(genome_names[i], binB[i], {"distance": dist[i], "hashes": hashes[i], "p-value": pvalue[i]})] for i in range(
            len(genome_names))]
        [gr.add_edges_from(i) for i in edges]
        # Filter a list of edges from the graph that are greater than the distance threshold
        long_edges = list(filter(lambda e: e[2] > self.distance_value, (e for e in gr.edges.data('distance'))))
        le_ids = list(e[:2] for e in long_edges)
        gr.remove_edges_from(le_ids)
        # Filter a list of edges from the graph that are less than the hashes threshold
        hash_edges = list(filter(lambda e: e[2] < self.hash_threshold, (e for e in gr.edges.data('hashes'))))
        he_ids = list(e[:2] for e in hash_edges)
        # Filter a list of edges from the graph that are greater than the p-value threshold
        pval_edges = list(filter(lambda e: e[2] > self.p_threshold, (e for e in gr.edges.data('p-value'))))
        pval_ids = list(e[:2] for e in pval_edges)
        gr.remove_edges_from(pval_ids)
        print(nx.info(gr))
        ''' Write out a graphml file for visualisation in cytoscape'''
        nx.write_graphml(gr, f"{self.fpath}networkx.xml")
        return gr

    def do_partitions(self, gr) -> pd.DataFrame():
        '''Use Louvain method to get partitions, save as JSON'''
        df = pd.DataFrame([louvain.best_partition(gr,weight='distance')]).T
        if self.save_partition_data == True:
            df.to_json(f"{self.fpath}full_louvain_partition.json")
        return df

    def output_conditioning(self, df) -> pd.DataFrame():
        '''Produce structured data, save'''
        unique, counts = np.unique(df, return_counts=True)
        print(f"Saving output. There are {len(unique)} communities detected")
        genomes = [df[0].iloc[np.where(np.isin(df.values, int(i)))[0]].index.to_list() for i in unique]
        output_table = pd.DataFrame()
        output_table["community_id"] = unique
        output_table["n_genomes"] = counts
        output_table["threshold"] = self.distance_value
        output_table.to_csv(f"{self.fpath}community_counts_table.csv")
        output_table["accession_ids"] = genomes
        output_table["accession_ids"] = output_table["accession_ids"].apply(lambda x: str(x).replace("[","").replace("]","").replace("'",""))
        output_table.to_csv(f"{self.fpath}community_members_table.csv")
        return output_table #len(unique)

    def output_mapping(self, output_table) -> None:
        '''Map to current BVS taxonomy, save'''
        # load the data
        communities_df = pd.DataFrame(output_table)
        taxonomy_df = pd.read_csv(self.taxonomy_file_path, header=0)
        # split the accession_ids field into lists
        communities_df['accession_ids'] = communities_df['accession_ids'].str.split(', ')
        # explode the lists
        communities_df = communities_df.explode('accession_ids', ignore_index = True)
        # merge and export
        mapped_results_df = pd.merge(communities_df, taxonomy_df[['GENBANK accession', 'Genus', 'Type']], left_on='accession_ids', right_on="GENBANK accession", how='left')
        mapped_results_df.to_csv(f"{self.fpath}mapped_taxa.csv", encoding='utf-8', index=False)
        # Count the number of unique communities by genus name
        count_unique = mapped_results_df.groupby('Genus')['community_id'].nunique().reset_index(name='N_communities')
        # Obtain a list of the community id(s) for each genus and export
        genus_communities_df = mapped_results_df.groupby('Genus')['community_id'].unique().apply(list).reset_index(name='community_ids')
        genus_communities_df = genus_communities_df.merge(count_unique[['Genus', 'N_communities']], left_on='Genus', right_on='Genus', how='left') 
        genus_communities_df.to_csv(f"{self.fpath}genus_community_ids.csv", encoding='utf-8', index=False)
    
    def calculate_layout(self, gr) -> dict:
        '''Get NetworkX graph layout'''
        print("Calculating graph layout. This may take a while...")
        return nx.spring_layout(gr)

    def main(self) -> str:
        '''Entrypoint'''
        st = time.time()
        print(f"Starting job at {st}\nDoing data extraction, loading & transformation")

        '''Clean & structure data'''
        df_clean = self.load_and_clean()
        graph_object = self.make_graph(df_clean[["binA", "binB", "distance"]])
        partition_df = self.do_partitions(graph_object)
        p = self.output_conditioning(partition_df)
        d = self.output_mapping(p)

        if self.plot_graph == True:
            '''Do dendrogram'''
            th = time.time()
            if self.trunc == "lastp":
                dendro_labels = partition_df[0]
                p_value = p
            else:
                dendro_labels = partition_df.index
                p_value = 0
            create_dendrogram(partition_df, self.fpath, self.sample_size, self.dendro_threshold, self.trunc, self.dendro_width, p=p_value, labels=dendro_labels)
            print(f"Dendrogram time, {self.sample_size} samples = {time.time()-th}")

            # return # Only uncomment me if you really want to plot a massive network graph

            '''Do network graph'''
            pos = self.calculate_layout(graph_object)
            plot = ClusterPlot(
                graph_object, partition_df[0], pos, self.sample_size, self.fpath)
            plot.make_all()

        return f"Finished in {round(time.time() - st, 2)} Sec. Results saved to {self.fpath}"


if __name__ == "__main__":
    '''Dev use only'''
    ...
