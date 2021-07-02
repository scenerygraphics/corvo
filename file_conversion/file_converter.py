import scanpy as sc


class FileConverter:
    """
"""
    def __init__(self):
        original_file = "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad"
        results_file = "marrow_vr_processed.h5ad"

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

