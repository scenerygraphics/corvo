import scanpy as sc
from file_conversion.test.file_converter_test_2d import FileConverterTest2D
from file_conversion.selections_to_csv import SelectionsToCSV


class FileConverter:
    """
"""
    def __init__(self):
        original_file = "original_file_name"
        results_file = "vr_processed.h5ad"

        # file to add 3d umap to
        adata = sc.read_h5ad(original_file)

        adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        # adding csc layer
        adata.layers["X_csc"] = adata.X.tocsc()

        # saving categoricals as ordered lists (effectively maps)
        for i in adata.var_keys():
            try:
                sub_cat = []
                for category in adata.var[i].cat.categories:
                    sub_cat.append(category)
                adata.uns[i + "_categorical"] = sub_cat
            except AttributeError:
                pass

        for i in adata.obs_keys():
            try:
                sub_cat = []
                for category in adata.obs[i].cat.categories:
                    sub_cat.append(category)
                adata.uns[i + "_categorical"] = sub_cat
            except AttributeError:
                pass

        for i in adata.uns_keys():
            try:
                sub_cat = []
                for category in adata.uns[i].cat.categories:
                    sub_cat.append(category)
                adata.uns[i + "_categorical"] = sub_cat
            except AttributeError:
                pass

        adata.write(results_file)


if __name__ == "__main__":
    FileConverter()

