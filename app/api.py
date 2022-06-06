from fastapi import FastAPI
from .main import GetPartitions
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel

app = FastAPI()


class PayloadData(BaseModel):
    sample_size: int = 10000
    p_value: float = 0.05
    input_file_path: str = "./MASH_dist_01Mar2022.tsv"
    save_partition_data: bool = True
    plot_graph: bool = True


@app.get("/")
async def health_check():
    return {"message": "Healthy"}


@app.post("/do_partitions/")
async def do_partitions(payload: PayloadData):
    params = jsonable_encoder(payload)
    gp = GetPartitions(params)
    status = gp.main()
    return status
