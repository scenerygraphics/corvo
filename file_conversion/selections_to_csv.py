import numpy as np
import pandas as pd
import scanpy as sc


class SelectionToCSV:
    def __init__(self):
        adata = sc.read_h5ad("temporary.h5ad")

        adata.write_csvs("test_export.csv", True, ",")  # write all annotations to csv
        umap_df = pd.DataFrame(adata.obsm["X_umap"])
        umap_df.to_csv("test_export/obsm.csv")  # overwrite obsm with only 3d umap


if __name__ == "__main__":
    SelectionToCSV()
