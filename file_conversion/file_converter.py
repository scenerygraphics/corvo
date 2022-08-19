#!/usr/bin/env python3

import scanpy as sc
import sys

# https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb
# https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728
class FileConverter:
    """
"""
    def __init__(self):
        original_file = sys.argv[1]
        results_file = sys.argv[1].rstrip(".h5ad") + "_vr_processed.h5ad"

        # file to add 3d umap to
        adata = sc.read_h5ad(original_file)

        adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        # making X sparse in csc format
        adata.X = adata.X.tocsc()

        adata.write(results_file)


if __name__ == "__main__":
    FileConverter()

