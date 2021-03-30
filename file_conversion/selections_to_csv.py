import numpy as np
import pandas as pd
import scanpy as sc
import csv


class SelectionsToCSV:
    def __init__(self, file_path, tissue):
        export_folder = "file_conversion/" + tissue + "_annotations"

        adata = sc.read_h5ad(file_path)

        adata.write_csvs(dirname=export_folder, skip_data=True, sep=",")  # write all annotations to csv
        umap_df = pd.DataFrame(adata.obsm["X_umap"])
        umap_df.to_csv(export_folder + "/obsm.csv")  # overwrite obsm with only 3d umap


if __name__ == "__main__":
    SelectionsToCSV("temporary.h5ad", "liver")
