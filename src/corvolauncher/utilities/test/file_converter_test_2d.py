import scanpy as sc


class FileConverterTest2D:
    def __init__(self, original, result):
        self.adata = sc.read_h5ad(original, True)
        self.adata_out = sc.read_h5ad(result, True)

        sum_list = self.adata.obsm["X_umap"][:, 0] - self.adata_out.obsm["X_umap_2d"][:, 0]
        diff_sum = 0
        for i in sum_list:
            diff_sum += i

        assert 0.5 > diff_sum > -0.5, "difference too large"

    def show_plots_ontology(self):
        sc.pl.umap(self.adata, color=['cell_ontology_class'])
        sc.pl.umap(self.adata_out, color=['cell_ontology_class'])


if __name__ == "__main__":
    FileConverterTest2D()
