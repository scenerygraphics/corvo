import pathlib
import scanpy as sc


class PreProcess:
    """
"""
    def __init__(self, parent, dataset):
        self.parent = parent
        self.dataset: str = dataset

        self.process_task()

    def process_task(self):
        adata = sc.read_h5ad("../resources/datasets/" + self.dataset)
        # adata.var_names = adata.var["feature_name"]
        if self.parent.shutdown:
            return
            # file to add 3d umap to
        out_file: str = self.dataset.rstrip(".h5ad") + "_processed.h5ad"

        # will fail for obs with only one category or when a category has only one data point
        n_genes = 10

        for observation in adata.obs:
            if self.parent.shutdown:
                break
            try:

                sc.tl.rank_genes_groups(adata, groupby=observation, n_genes=n_genes, method='t-test', use_raw=True)
                if self.parent.shutdown:
                    break
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
                    if self.parent.shutdown:
                        break
                    adata.uns[(observation + "_" + cat[0])] = cat[1]

            # pass for any obs that do not cluster
            except (AttributeError, ValueError, ZeroDivisionError):
                print("fail: " + observation)
                pass

        if self.parent.shutdown:
            return
        # adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        try:
            adata.uns['neighbors']
        except Exception:
            sc.pp.neighbors(adata)

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        if self.parent.shutdown:
            return

        # making X sparse in csc format
        adata.X = adata.X.tocsc()

        if self.parent.shutdown:
            return

        adata.write(pathlib.Path("../resources/processed_datasets/" + out_file))
