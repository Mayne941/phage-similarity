from fastapi import FastAPI
from .main import GetPartitions
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel

app = FastAPI()


class PayloadData(BaseModel):
    distance_value: float = 0.05
    input_file_path: str = "./MASH_dist_01Mar2022.tsv"
    save_partition_data: bool = True
    plot_graph: bool = False
    sample_size: int = 0
    dendro_threshold: int = 1.5
    dendro_truncation: str = "none"
    dendro_width = int = 2500


@app.get("/")
async def health_check():
    return {"message": "Healthy"}


@app.post("/do_partitions/")
async def do_partitions(payload: PayloadData):
    params = jsonable_encoder(payload)
    gp = GetPartitions(params)
    status = gp.main()
    return status
