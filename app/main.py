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
        self.p_value = params["p_value"]
        self.input_file_path = params["input_file_path"]
        self.save_partition_data = params["save_partition_data"]
        self.plot_graph = params["plot_graph"]
        self.dendro_threshold = params["dendro_threshold"]
        self.trunc = params["dendro_truncation"]
        self.dendro_width = params["dendro_width"]
        if not os.path.exists("./output/"):
            os.mkdir("./output/")
        self.fpath = f"./output/exp_{datetime.now().strftime('%d-%b-%Y--%H-%M-%S')}/"
        if not os.path.exists(self.fpath):
            os.mkdir(self.fpath)

    def load_and_clean(self) -> pd.DataFrame:
        '''Read in data, filter duplicates and by p value threshold'''
        df = pd.read_csv(self.input_file_path, delimiter="\t", engine="python")
        df_no_dups = df[~(df["binA"] == df["binB"])]
        df_no_dups = df_no_dups[~(df_no_dups["tX1000"] == "1000/1000")]
        df_p_filter = df_no_dups[df_no_dups["distance"] <= self.p_value]
        if self.sample_size <= 0:
            self.sample_size = df_p_filter.shape[0]
        return df_p_filter.head(self.sample_size)

    def make_graph(self, df) -> nx.Graph():
        '''Make NetworkX graph object'''
        genome_names = df["binA"].tolist()
        binB = df["binB"].tolist()
        dist = df["distance"].tolist()
        gr = nx.Graph()

        [gr.add_node(i) for i in df["binA"].unique().tolist()]
        edges = [[(genome_names[i], binB[i], {"distance": dist[i]})] for i in range(
            len(genome_names))]
        [gr.add_edges_from(i) for i in edges]
        return gr

    def do_partitions(self, gr) -> pd.DataFrame():
        '''Use Louvain method to get partitions, save as JSON'''
        partition = louvain.best_partition(gr)
        df = pd.DataFrame([partition]).T
        if self.save_partition_data == True:
            df.to_json(f"{self.fpath}full_louvain_partition.json")
        return df

    def output_conditioning(self, df) -> None:
        '''Produce structured data, save'''
        unique, counts = np.unique(df, return_counts=True)
        print(f"Saving output. There are {len(unique)} communities detected")
        genomes = [df[0].iloc[np.where(np.isin(df.values, int(i)))[0]].index.to_list() for i in unique]
        output_table = pd.DataFrame()
        output_table["community_id"] = unique
        output_table["n_genomes"] = counts
        output_table["threshold"] = self.p_value
        output_table.to_csv(f"{self.fpath}community_counts_table.csv")
        output_table["accession_ids"] = genomes
        output_table["accession_ids"] = output_table["accession_ids"].apply(lambda x: str(x).replace("[","").replace("]","").replace("'",""))
        output_table.to_csv(f"{self.fpath}community_members_table.csv")
        return len(unique)

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

        if self.plot_graph == True:
            '''Do dendrogram'''
            th = time.time()
            if self.trunc == "lastp":
                create_dendrogram(partition_df, self.fpath, self.sample_size, self.dendro_threshold, self.trunc, self.dendro_width, p=p, labels=partition_df[0])
            else:
                create_dendrogram(partition_df, self.fpath, self.sample_size, self.dendro_threshold, self.trunc, self.dendro_width, p=0, labels=partition_df.index)
            print(f"Dendrogram time, {self.sample_size} samples = {time.time()-th}")

            return # Only uncomment me if you really want to plot a massive network graph

            '''Do network graph'''
            pos = self.calculate_layout(graph_object)
            plot = ClusterPlot(
                graph_object, partition_df[0], pos, self.sample_size, self.fpath)
            plot.make_all()

        return f"Finished in {round(time.time() - st, 2)} Sec. Results saved to {self.fpath}"


if __name__ == "__main__":
    '''Dev use only'''
    params = {"sample_size": 0,
              "p_value": 0.1,
              "input_file_path": "./MASH_dist_1Jun2022_0.3.k31.tsv", #./MASH_dist_01Mar2022.tsv",
              "save_partition_data": True,
              "plot_graph": True,
              "dendro_threshold": 1.5,
              "dendro_truncation": "lastp", # "none" is subs kw for None
              "dendro_width": 25000} 
    gp = GetPartitions(params)
    status = gp.main()
    print(status)
