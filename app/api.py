from fastapi import FastAPI
from app.main import GetPartitions
from fastapi.encoders import jsonable_encoder

from app.utils.api_classes import LouvainData

app = FastAPI()




@app.get("/")
async def health_check():
    return {"message": "Healthy"}


@app.post("/louvain_partitions/")
async def do_partitions(payload: LouvainData):
    params = jsonable_encoder(payload)
    gp = GetPartitions(params, is_louvain = True)
    status = gp.main()
    return status

@app.post("/leiden_partitions/")
async def do_partitions(payload: LouvainData):
    params = jsonable_encoder(payload)
    gp = GetPartitions(params)
    status = gp.main()
    return status
