from pydantic import BaseModel

class LouvainData(BaseModel):
    distance_value: float = 0.15
    input_file_path: str = "./2Jul2023.d0.3.k15.s25000.tsv"
    taxonomy_file_path: str = "./taxonomy.csv"
    save_partition_data: bool = True
    plot_graph: bool = False
    sample_size: int = 0
    dendro_threshold: int = 1.5
    dendro_truncation: str = "none"
    dendro_width: int = 2500
    hash_threshold: float = 0.15
    p_threshold: float = 1e-10
