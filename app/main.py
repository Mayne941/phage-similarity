import pandas as pd
import os
import numpy as np
import networkx as nx
import time
from datetime import datetime

from app.cluster_graph import ClusterPlot
from app.dendrogram import create_dendrogram
from app.utils.call_graph_objs import do_leiden_partitions, do_louvain_partitions

class GetPartitions:
    def __init__(self, params, is_louvain = False) -> None:
        '''Set up exp folder, decode runtime params'''
        self.is_louvain = is_louvain
        self.params = params
        if not os.path.exists("./output/"):
            os.mkdir("./output/")
        self.fpath = f"./output/exp_{datetime.now().strftime('%d-%b-%Y--%H-%M-%S')}/"
        if not os.path.exists(self.fpath):
            os.mkdir(self.fpath)

    def load_and_clean(self) -> pd.DataFrame:
        '''Read in data, filter duplicates and by p value threshold'''
        try:
            df = pd.read_csv(self.params["input_file_path"], delimiter="\t", engine="python", header=0, names=["binA", "binB", "distance", "p-value", "tX"]) 
        except:
            raise FileNotFoundError(f'Cant find your input file: {self.params["input_file_path"]}\nCheck that file exists and is in correct format')
        df_no_dups = df[~(df["binA"] == df["binB"])]
        '''Hash filter'''
        df_no_dups["t"] = df_no_dups["tX"].str.split("/").str[0].astype(int)
        df_no_dups["X"] = df_no_dups["tX"].str.split("/").str[1].astype(int)
        df_no_dups["hashes"] = df_no_dups["t"]/df_no_dups["X"]
        '''Sample size filter'''
        if self.params["sample_size"] <= 0:
            self.params["sample_size"] = df_no_dups.shape[0]
        return df_no_dups.head(self.params["sample_size"])

    def make_nx_graph(self, df) -> nx.Graph():
        '''Make NetworkX graph object'''
        genome_names = df["binA"].tolist()  
        binB = df["binB"].tolist()
        dist = df["distance"].tolist()
        hashes = df["hashes"].tolist()
        pvalue = df["p-value"].tolist()
        gr = nx.Graph()
        [gr.add_node(i) for i in df["binA"].unique().tolist()]
        edges = [[(genome_names[i], binB[i], {"distance": dist[i], "hashes": hashes[i], "p-value": pvalue[i]})] for i in range(
            len(genome_names))]
        [gr.add_edges_from(i) for i in edges]

        '''Filter a list of edges from the graph that are greater than the distance threshold'''
        long_edges = list(filter(lambda e: e[2] > self.params["distance_value"], (e for e in gr.edges.data('distance'))))
        le_ids = list(e[:2] for e in long_edges)
        gr.remove_edges_from(le_ids)

        '''Filter a list of edges from the graph that are greater than the p-value threshold'''
        pval_edges = list(filter(lambda e: e[2] > self.params["p_threshold"], (e for e in gr.edges.data('p-value'))))
        pval_ids = list(e[:2] for e in pval_edges)
        gr.remove_edges_from(pval_ids)
        
        '''Filter a list of edges from the graph that are less than the hashes threshold'''
        hash_edges = list(filter(lambda e: e[2] < self.params["hash_threshold"], (e for e in gr.edges.data('hashes'))))
        he_ids = list(e[:2] for e in hash_edges)
        gr.remove_edges_from(he_ids)
                
        print(nx.info(gr))

        ''' Write out a graphml file for visualisation in cytoscape'''
        nx.write_graphml(gr, f"{self.fpath}networkx.xml")
        return gr, do_louvain_partitions(gr, self.params["save_partition_data"], self.fpath)

    def make_i_graph(self, df) -> nx.Graph():
        '''Make NetworkX graph object'''
        genome_names = df["binA"].tolist() 
        binB = df["binB"].tolist()
        dist = df["distance"].tolist()
        hashes = df["hashes"].tolist()
        pvalue = df["p-value"].tolist()
        df.to_csv(f"{self.fpath}passed_to_networkx.tsv", sep='\t', encoding='utf-8') 
        gr = nx.Graph()
        [gr.add_node(i) for i in df["binA"].unique().tolist()]
        edges = [[(genome_names[i], binB[i], {"distance": dist[i], "hashes": hashes[i], "p-value": pvalue[i]})] for i in range(
            len(genome_names))]
        [gr.add_edges_from(i) for i in edges]

        '''Filter a list of edges from the graph that are greater than the distance threshold'''
        long_edges = list(filter(lambda e: e[2] > self.params["distance_value"], (e for e in gr.edges.data('distance'))))
        le_ids = list(e[:2] for e in long_edges)
        gr.remove_edges_from(le_ids)

        '''Filter a list of edges from the graph that are greater than the p-value threshold'''
        pval_edges = list(filter(lambda e: e[2] > self.params["p_threshold"], (e for e in gr.edges.data('p-value'))))
        pval_ids = list(e[:2] for e in pval_edges)
        gr.remove_edges_from(pval_ids)

        '''Filter a list of edges from the graph that are less than the hashes threshold'''
        hash_edges = list(filter(lambda e: e[2] < self.params["hash_threshold"], (e for e in gr.edges.data('hashes'))))
        he_ids = list(e[:2] for e in hash_edges)
        gr.remove_edges_from(he_ids)
        
        print(nx.info(gr))
        ''' Write out a graphml file for visualisation in cytoscape'''
        nx.write_graphml(gr, f"{self.fpath}networkx.xml")
        return gr, do_leiden_partitions(gr, self.params["save_partition_data"], self.fpath)

    def output_conditioning(self, df) -> pd.DataFrame():
        '''Produce structured data, save'''
        if self.is_louvain:
            genome_match_df = df[0]
        else:
            df.set_index([0], inplace=True)
            genome_match_df = df.copy()
        unique, counts = np.unique(genome_match_df, return_counts=True)
        genomes = [genome_match_df.iloc[np.where(np.isin(genome_match_df.values, int(i)))[0]].index.to_list() for i in unique]
        print(f"Saving output. There are {len(unique)} communities detected")

        output_table = pd.DataFrame()
        output_table["community_id"] = unique
        output_table["n_genomes"] = counts
        output_table["distance_threshold"] = self.params["distance_value"]
        output_table["pvalue_threshold"] = self.params["p_threshold"]
        output_table["hashes_threshold"] = self.params["hash_threshold"]
        output_table.to_csv(f"{self.fpath}community_counts_table.csv")
        output_table["accession_ids"] = genomes
        output_table["accession_ids"] = output_table["accession_ids"].apply(lambda x: ", ".join(x))
        output_table.to_csv(f"{self.fpath}community_members_table.csv")
        return output_table, len(unique)

    def output_mapping(self, output_table) -> None:
        '''Map to current BVS taxonomy, save'''
        '''load the data & custom VMR with strain list'''
        communities_df = pd.DataFrame(output_table)
        try:
            taxonomy_df = pd.read_csv(self.params["taxonomy_file_path"], header=0)
        except:
            raise FileNotFoundError(f"Can't open your taxonomy file: {self.params['taxonomy_file_path']}\nCheck that it exists and is in csv format")
        
        '''split the accession_ids field into lists, explode the lists, merge and export'''
        communities_df['accession_ids'] = communities_df['accession_ids'].str.split(', ')
        communities_df = communities_df.explode('accession_ids', ignore_index = True)
        mapped_results_df = pd.merge(communities_df, taxonomy_df[['GENBANK accession', 'Genus', 'Type']], left_on='accession_ids', right_on="GENBANK accession", how='left')
        mapped_results_df.to_csv(f"{self.fpath}mapped_taxa.csv", encoding='utf-8', index=False)

        '''Count the number of unique communities by genus name'''
        count_unique = mapped_results_df.groupby('Genus')['community_id'].nunique().reset_index(name='N_communities')

        '''Obtain a list of the community id(s) for each genus and export'''
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
        if self.is_louvain:
            graph_object, partition_df = self.make_nx_graph(df_clean[["binA", "binB", "distance", "hashes", "p-value"]])
        else: 
            graph_object, partition_df = self.make_i_graph(df_clean[["binA", "binB", "distance", "hashes", "p-value"]])

        output_table, p = self.output_conditioning(partition_df)
        self.output_mapping(output_table)

        if self.params["plot_graph"] == True:
            '''Do dendrogram'''
            th = time.time()
            if self.params["dendro_truncation"] == "lastp":
                dendro_labels = partition_df[0]
                p_value = p
            else:
                dendro_labels = partition_df.index
                p_value = 0
            create_dendrogram(partition_df, self.fpath, self.params["sample_size"], self.params["dendro_threshold"], self.params["dendro_truncation"], self.params["dendro_width"], p=p_value, labels=dendro_labels)
            print(f"Dendrogram time, {self.params['sample_size']} samples = {time.time()-th}")

            # return # Only uncomment me if you really want to plot a massive network graph

            '''Do network graph'''
            pos = self.calculate_layout(graph_object)
            plot = ClusterPlot(
                graph_object, partition_df[0], pos, self.params["sample_size"], self.fpath)
            plot.make_all()

        msg = f"Finished in {round(time.time() - st, 2)} Sec. Results saved to {self.fpath}"
        print(msg)
        return msg


if __name__ == "__main__":
    '''Dev use only'''
    ...
