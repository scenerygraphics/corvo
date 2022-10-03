import pathlib
import scanpy as sc
import numpy as np
# get request default (give me index as html) Get only works for requests. Post request most commonly used for submitting form data. Could have been a get request as no additional data is sent. Could have generated a cookie based on get request to make difficult, but they don't.


class CompleteProcess:
    """
"""
    def __init__(self, dataset):
        # self.parent = parent
        self.dataset: str = dataset

        adata = sc.read_h5ad("../resources/datasets/" + dataset)
        # adata.var_names = adata.var["feature_name"]

        # file to add 3d umap to
        out_file: str = dataset.rstrip(".h5ad") + "_processed.h5ad"

        # will fail for obs with only one category or when a category has only one data point
        n_genes = 10

        for observation in adata.obs:
            try:

                sc.tl.rank_genes_groups(adata, groupby=observation, n_genes=n_genes, method='t-test', use_raw=True)

                print("success: " + observation)
                names_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                pvals_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                logfoldchanges_list = [0] * adata.obs[observation].cat.categories.size * n_genes

                for cat in [["names", names_list], ["pvals", pvals_list], ["logfoldchanges", logfoldchanges_list]]:

                    rank_counter = 0
                    for rank in adata.uns["rank_genes_groups"][cat[0]]:  # rank is list of 1st, 2nd, 3d most expr etc
                        name_counter = 0
                        if cat[0] == "pvals":
                            for index in rank:
                                cat[1][(name_counter * n_genes) + rank_counter] = round(index, 4)
                                name_counter += 1

                        else:
                            for index in rank:
                                cat[1][(name_counter * n_genes) + rank_counter] = index
                                name_counter += 1
                                # print(index)

                        rank_counter += 1

                    adata.uns[(observation + "_" + cat[0])] = cat[1]

            # pass for any obs that do not cluster
            except (AttributeError, ValueError, ZeroDivisionError):
                print("fail: " + observation)
                pass

        # adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        # making X sparse in csc format
        adata.X = adata.X.tocsc()

        adata.write(pathlib.Path("../resources/processed_datasets/" + out_file))


if __name__ == "__main__":
    CompleteProcess("marrow.h5ad")
