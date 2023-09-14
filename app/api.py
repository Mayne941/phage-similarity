from fastapi import FastAPI
from .main import GetPartitions
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel

app = FastAPI()


class PayloadData(BaseModel):
    distance_value: float = 0.15
    input_file_path: str = "./2Jul2023.d0.3.k15.s25000.tsv"
    taxonomy_file_path: str = "./taxonomy.csv"
    save_partition_data: bool = True
    plot_graph: bool = False
    sample_size: int = 0
    dendro_threshold: int = 1.5
    dendro_truncation: str = "none"
    dendro_width: int = 2500
    hash_threshold: float = 0.30
    p_threshold: float = 1e-10


@app.get("/")
async def health_check():
    return {"message": "Healthy"}


@app.post("/do_partitions/")
async def do_partitions(payload: PayloadData):
    params = jsonable_encoder(payload)
    gp = GetPartitions(params)
    status = gp.main()
    return status
