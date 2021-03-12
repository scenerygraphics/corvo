import scanpy as sc
from file_converter_test_2d import FileConverterTest2D
import h5py


class FileConverter:
    """
"""
    def __init__(self):
        original_file = "file_conversion/tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad"
        results_file = "file_conversion/temporary.h5ad"

        # file to add 3d umap to
        adata = sc.read_h5ad(original_file, True)

        adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()
        adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        sc.tl.umap(adata, 3)

        FileConverterTest2D(original_file, results_file)


if __name__ == "__main__":
    FileConverter()

