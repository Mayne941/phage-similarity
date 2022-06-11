# Grouping MASH results into communities
## Rich Mayne 2022

To install:

1. cd {repo directory}
1. python3 -m pip install -r requirements.txt

To run:
1. python3 -m uvicorn app.api:app --reload
1. Navigate to http://127.0.0.1:8000/docs in browser
1. Expand the "/do_partitions/" window, then hit "try it out".
1. Specify your run parameters, including defining the path to your MASH output file. Then, hit "execute".

Run parameters:
1. P value. Threshold or probability value, used to restrict the sensitivity of genome similarity, i.e. any rows with similarity scores > threshold will be discounted.
1. Input file path. Point to MASH output TSV file. N.b. uses Linux dot notation: recommended to place your data file in the top level directory of this repo.
1. Save partition data. Set to True to save an unmodified dataset for cluster groupings, e.g. for making Dendrograms. Saves as JSON.
1. Plot graph. Set to True to project clusters into 2 dimensions as a nx graph (spring layout by default) and plot. Will open an interactive figure in a browser and save a still PNG image. N.b. For large numbers of samples (>10k) this will take a very long time.
1. Sample size. Used to restrict input size: set to 0 for no restriction. N.b. Running full pipeline with graph generation, on large (>10k) samples will take a very long time.
1. Dendrogram threshold. Let *t* be the color_threshold. Colors all the descendent links below a cluster node *k* the same color if *k* is the first node below the cut threshold *t*. All links connecting nodes with distances greater than or equal to the threshold are colored with de default matplotlib color 'C0'. Ignored if plot_graph = False.
1. Dendrogram truncaton. The dendrogram can be hard to read when the original observation matrix from which the linkage is derived is large. Truncation is used to condense the dendrogram. Set as either "none" (NOT None) for no truncation, or "lastp". For the latter, the last p non-singleton clusters formed in the linkage are the only non-leaf nodes in the linkage; they correspond to rows Z[n-p-2:end] in Z. All other non-singleton clusters are contracted into leaf nodes. p is automatically set to the number of communities detected.
1. Dendrogram width. Manually expand or contract width of dendrogram to aid readability, in pixels.