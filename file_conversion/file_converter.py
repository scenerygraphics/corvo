import scanpy as sc
from file_conversion.test.file_converter_test_2d import FileConverterTest2D
from selections_to_csv import SelectionsToCSV


class FileConverter:
    """
"""
    def __init__(self):
        original_file = "file_conversion/tabula-muris-senis-facs-processed-official-annotations-Liver.h5ad"
        results_file = "file_conversion/liver_vr_processed.h5ad"

        # file to add 3d umap to
        adata = sc.read_h5ad(original_file)

        adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()
        adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]

        sc.tl.umap(adata, n_components=3)

        # adding csc layer
        adata.layers["X_csc"] = adata.X.tocsc()
        adata.write(results_file)

        # FileConverterTest2D(original_file, results_file)

        # annotations to csv
        SelectionsToCSV(results_file, "liver")


if __name__ == "__main__":
    FileConverter()

