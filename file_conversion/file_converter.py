#!/usr/bin/env python3

import scanpy as sc
import sys


class FileConverter:
    """
"""

    def __init__(self):
        original_file = sys.argv[1]
        results_file = sys.argv[1].rstrip(".h5ad") + "_vr_processed.h5ad"

        # file to add 3d umap to
        adata = sc.read_h5ad(original_file)

        n_genes = 10
        for observation in adata.obs:
            try:
                sc.tl.rank_genes_groups(adata, groupby=observation, n_genes=n_genes, method='wilcoxon', use_raw=True)

                names_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                pvals_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                logfoldchanges_list = [0] * adata.obs[observation].cat.categories.size * n_genes

                for cat in [["names", names_list], ["pvals", pvals_list], ["logfoldchanges", logfoldchanges_list]]:

                    rank_counter = 0
                    for rank in adata.uns["rank_genes_groups"][cat[0]]:  # rank is list of 1st, 2nd, 3d most expr etc
                        name_counter = 0

                        for index in rank:
                            cat[1][(name_counter * n_genes) + rank_counter] = index
                            name_counter += 1

                        rank_counter += 1

                    adata.uns[(observation + "_" + cat[0])] = cat[1]

            # pass for any obs that do not cluster
            except (AttributeError, ValueError, ZeroDivisionError):
                pass

        adata.uns["rank_genes_groups"] = []

        # print(adata.uns)
        # adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        # making X sparse in csc format
        adata.X = adata.X.tocsc()

        # saving categoricals as ordered lists (effectively maps)
        for i in adata.obs_keys():
            try:
                sub_cat = []
                for category in adata.obs[i].cat.categories:
                    sub_cat.append(category)
                adata.uns[i + "_categorical"] = sub_cat
            except AttributeError:
                pass

        adata.write(results_file)


if __name__ == "__main__":
    FileConverter()
