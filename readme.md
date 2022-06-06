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