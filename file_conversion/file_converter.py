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

        for observation in adata.obs:
            try:
                sc.tl.rank_genes_groups(adata, groupby=observation, n_genes=10, method='wilcoxon', use_raw=True)

                names_list = []
                for rank in adata.uns["rank_genes_groups"]["names"]:
                    for name in rank:
                         names_list.append(name)

                adata.uns[(observation + "_names")] = names_list

                pvals_list = []
                for rank in adata.uns["rank_genes_groups"]["pvals"]:
                    for pval in rank:
                        pvals_list.append(pval)

                adata.uns[(observation + "_pvals")] = pvals_list

                logfoldchanges_list = []
                for rank in adata.uns["rank_genes_groups"]["logfoldchanges"]:
                    for logfoldchange in rank:
                        logfoldchanges_list.append(logfoldchange)

                adata.uns[(observation + "_logfoldchanges")] = logfoldchanges_list

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

