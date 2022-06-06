# Grouping MASH results into communities
## Rich Mayne 2022

To install:

1. cd to PHAGE-CODING
1. python3 -m pip install -r requirements.txt

To run:
1. python3 -m uvicorn app.api:app --reload
1. Navigate to http://127.0.0.1:8000/docs in browser
1. Expand the "/do_partitions/" window, then hit "try it out".
1. Specify your run parameters, including defining the path to your MASH output file. Then, hit "execute".
