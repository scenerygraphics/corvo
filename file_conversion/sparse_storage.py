import h5py
import numpy as np
import pandas as pd
import scanpy as sc


class SparseStorage:
    def __init__(self):
        # f = h5py.File("file_conversion/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad")

        # f["X"].visititems(print)
        adata = sc.read_h5ad("tabula-muris-senis-facs-processed-official-annotations-Kidney.h5ad")
        print(adata)
        print(adata.var["highly_variable"][1])
        # print(len(adata["X"]))

        # print(len(adata.obs["n_genes"]))
        # print(len(adata.var["n_cells"]))
        # adata.write_zarr("aorta.zarr", chunks=(0, len(adata.obs["n_genes"])))


if __name__ == "__main__":
    SparseStorage()
