# Grouping MASH results into communities
## Rich Mayne & Dann Turner 2022

### Description
MASHUP creates a networkx object from the MASH distances, where the nodes are genomes and the edges are the calculated MASH distance. The matching hashes and p-values are also associated with each edge. Edges that do not meet the user specified thresholds are dropped from the network prior to paritioning into communities using the Louvain algorithm. Each genome is then mapped against the existing ICTV BVS taxonomy. Six output files are created:

| Output File                                     | Description                                                  |
| ----------------------------------------------- | ------------------------------------------------------------ |
|community_counts_table.csv                       |CSV file detailing community id, number of genomes and user supplied thresholds|
|community_members_table.csv                      |CSV file detailing community id, number of genomes, user supplied thresholds and a comma-separate list of accession numbers for members of each community |
|mapped_taxa.csv                                  |Exploded list where each row represents a single community member (accession number) with mapped taxonomy data (genus, species or strain) |
|genus_community_ids.csv                          |The number of identities of communities mapping to each genus in the taxonomy data | 
|networkx.xml                                     |The MASH distance network after dropping edges based on user thresholds. This can be imported directly into network visualisation software such as Cytoscape |
|full_louvain_parition.json                       |The output of the Louvain community paritioning dataframe in json format | 

### To install:
1. cd {repo directory}
1. python3 -m pip install -r requirements.txt

### To run:
1. MASHup takes the results of an all-vs-all MASH comparison of the INPHARED dataset.
Download the lastest [INPHARED genomes excluding RefSeq file](https://github.com/RyanCook94/inphared), create a MASH sketch file and distance comparisons. This step does take a significant amount of time due to the number of comparisons performed.
```bash
mash sketch -i -k 15 -s 25000 5Jan2023_genomes_excluding_refseq.fa
mash dist -i -d 0.3 5Jan2023_genomes_excluding_refseq.msh 5Jan2023_genomes_excluding_refseq.mash > 5Jan2023.d0.3.k15.s25000.tsv
```  
1. Copy the MASH distance file to the root of the cloned repository.
1. python3 -m uvicorn app.api:app --reload
1. Navigate to http://127.0.0.1:8000/docs in browser
1. Expand the "/do_partitions/" window, then hit "try it out".
1. Specify your run parameters, including defining the path to your MASH output file. Then, hit "execute".

### Run parameters:
1. Distance value. The MASH distance threshold, used to restrict the sensitivity of genome similarity, i.e. any edges with similarity scores > threshold will be dropped.
1. Hash threshold. Any edges where the fraction of matching hashes is < than the threshold will be dropped.
1. p_threshold. Any edges where the p-value is > than the threshold will be dropped. 
1. Input file path. Point to MASH output TSV file. N.b. uses Linux dot notation: recommended to place your data file in the top level directory of this repository.
1. Taxonomy file path. Point to the supplied taxonomy csv file. 
1. Save partition data. Set to True to save an unmodified dataset for cluster groupings, e.g. for making Dendrograms. Saves as JSON.
1. Plot graph. Set to True to project clusters into 2 dimensions as a nx graph (spring layout by default) and plot. Will open an interactive figure in a browser and save a still PNG image. N.b. For large numbers of samples (>10k) this will take a very long time.
1. Sample size. Used to restrict input size: set to 0 for no restriction. N.b. Running full pipeline with graph generation, on large (>10k) samples will take a very long time.
1. Dendrogram threshold. Let *t* be the color_threshold. Colors all the descendent links below a cluster node *k* the same color if *k* is the first node below the cut threshold *t*. All links connecting nodes with distances greater than or equal to the threshold are colored with de default matplotlib color 'C0'. Ignored if plot_graph = False.
1. Dendrogram truncaton. The dendrogram can be hard to read when the original observation matrix from which the linkage is derived is large. Truncation is used to condense the dendrogram. Set as either "none" (NOT None) for no truncation, or "lastp". For the latter, the last p non-singleton clusters formed in the linkage are the only non-leaf nodes in the linkage; they correspond to rows Z[n-p-2:end] in Z. All other non-singleton clusters are contracted into leaf nodes. p is automatically set to the number of communities detected.
1. Dendrogram width. Manually expand or contract width of dendrogram to aid readability, in pixels.
