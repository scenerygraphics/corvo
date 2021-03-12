import pandas as pd
import numpy as np
import scanpy as sc


class PandasReindex:
    def __init__(self):
        adata = sc.read_h5ad("file_conversion/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad")

        pass


if __name__ == "__main__":
    PandasReindex()
