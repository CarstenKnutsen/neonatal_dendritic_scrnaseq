"""Goal:Initial analysis of dendritic cell scRNAseq, creating h5ad, filtering, and embedding
Date:241220
Author: Carsten Knutsen
conda_env:dendritic
"""

import pandas as pd
import os
import scanpy as sc
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from scipy.stats import median_abs_deviation
import scanpy.external as sce
import itertools
import string
from gtfparse import read_gtf
import anndata
from collections import defaultdict

gtf_fn = "data/genes.gtf"
input_data = "data/single_cell_files/cellranger_output"
output_data = "data/single_cell_files/scanpy_files"
os.makedirs(output_data, exist_ok=True)
adata_name = "dendritic"
figures = "data/figures"
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures
gene_dict = {
    "immune": [
        "Cd68",
        "Gal",
        "Itgax",
        "Car4",
        "C1qa",
        "C1qb",
        "Plac8",
        "Batf3",
        "Itgae",
        "Cd209a",
        "Mreg",
        # "Mcpt8",
        "Retnlg",
        "Ms4a1",
        "Gzma",
        "Cd3e",
        "Areg",
        "Mki67",
        "Col1a1",
        "Epcam",
        "Ptprc",
        "Cdh5",
        "Cd14",
        "Fcgr3",
        "Il3ra",
    ],
}
leiden_ct_dict = {
    "0": "DC",
    "1": "DC",
    "2": "DC",
    "3": "DC",
    "4": "DC",
    "5": "DC",
    "6": "DC",
    "7": "DC",
    "8": "DC",
    "9": "DC",
    "10": "DC",
    "11": "T cell",
    "12": "Alveolar macrophage",
    "13": "Monocyte",
    "14": "B cell",
    "15": "DC/T doublet",
    "16": "Neutrophil",
    "17": "Endothelial",
}
if __name__ == "__main__":
    runs = os.listdir(input_data)
    adatas = []
    gtf = read_gtf(gtf_fn)
    for run in runs:
        print(run)
        if run == "Undetermined":
            continue
        folder = f"{input_data}/{run}/outs"
        adata = sc.read_10x_h5(f"{folder}/filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()
        adata.obs_names = run + "_" + adata.obs_names
        adata.obs["Treatment"] = run[0]
        adata.obs["Library"] = run
        adata.obs["Sex"] = run[-1]
        adatas.append(adata.copy())
    adata = anndata.concat(adatas)
    adata.obs["Treatment"].replace({"H": "Hyperoxia", "N": "Normoxia"}, inplace=True)
    for column in ["gene_id", "gene_type", "seqname", "transcript_name", "protein_id"]:
        temp_dict = pd.Series(gtf[column], index=gtf["gene_name"]).to_dict()
        temp_dict.update(pd.Series(gtf[column], index=gtf["gene_id"]).to_dict())
        temp_dict = defaultdict(lambda: None, temp_dict)
        adata.var[column] = [temp_dict[x] for x in adata.var.index]
    adata.var["mt"] = [True if x == "chrM" else False for x in adata.var["seqname"]]
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains(("^Hb[^(p)]"))
    print(adata)
    adata.write(
        f"{output_data}/{adata_name}_cellranger_all_cells.gz.h5ad", compression="gzip"
    )
    figures_qc = f"{figures}/qc"
    os.makedirs(figures_qc, exist_ok=True)
    sc.settings.figdir = figures_qc

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        expr_type="umis",
        percent_top=[20],
        log1p=True,
        inplace=True,
    )

    sc.pl.scatter(
        adata,
        x="log1p_total_umis",
        y="log1p_n_genes_by_umis",
        show=False,
        save="_genes_by_counts_pretrim_log",
    )
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f"_pretrim")
    sc.pl.violin(
        adata,
        ["log1p_total_umis"],
        groupby="Library",
        rotation=90,
        show=False,
        save="_umis_pretrim",
    )
    sc.pl.violin(
        adata,
        ["log1p_n_genes_by_umis"],
        groupby="Library",
        rotation=90,
        show=False,
        save="_genes_pretrim",
    )
    sc.pl.violin(
        adata,
        ["pct_umis_mt"],
        groupby="Library",
        rotation=90,
        show=False,
        save="_mt_pretrim",
    )
    # We are running scrublet here to detect any doublets on the raw count matrix
    sc.pp.scrublet(adata, batch_key="Library")
    print(adata)
    sc.pl.scrublet_score_distribution(adata, show=False, save="scrublet_scores")

    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_umis", 5)
        | is_outlier(adata, "log1p_n_genes_by_umis", 5)
        | is_outlier(adata, "pct_umis_in_top_20_genes", 5)
    )  # https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#
    adata.obs["mt_outlier"] = adata.obs["pct_umis_mt"] > 30
    print(adata.obs.mt_outlier.value_counts())
    print(adata.obs.outlier.value_counts())
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pl.scatter(
        adata,
        x="log1p_total_umis",
        y="log1p_n_genes_by_umis",
        show=False,
        save="_genes_by_umis_posttrim_log",
    )
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f"_posttrim")
    sc.pl.violin(
        adata,
        ["log1p_total_umis"],
        rotation=90,
        show=False,
        groupby="Library",
        save="_umis_posttrim",
    )
    sc.pl.violin(
        adata,
        ["log1p_n_genes_by_umis"],
        rotation=90,
        show=False,
        groupby="Library",
        save="_genes_posttrim",
    )
    sc.pl.violin(
        adata,
        ["pct_umis_mt"],
        groupby="Library",
        rotation=90,
        show=False,
        save="_mt_posttrim",
    )
    print(f"Median UMIS: {adata.obs['total_umis'].median()}")
    print(f"Median Genes: {adata.obs['n_genes_by_umis'].median()}")
    figures_embed = f"{figures}/initial_embedding"
    os.makedirs(figures_embed, exist_ok=True)
    sc.settings.figdir = figures_embed
    adata.layers["raw"] = adata.X.copy()
    sc.pp.normalize_total(adata, key_added=None, target_sum=1e4)
    adata.layers["cp10k"] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, batch_key="Library")
    sc.pp.pca(adata, use_highly_variable=True)
    # sce.pp.harmony_integrate(adata, 'Library',adjusted_basis='X_pca')
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, key_added="leiden", resolution=0.5)
    adata.obs["celltype_rough"] = [leiden_ct_dict[x] for x in adata.obs["leiden"]]
    sc.pl.dotplot(
        adata,
        gene_dict["immune"],
        groupby="leiden",
        dendrogram=False,
        show=False,
        save="useful_genes.png",
    )
    sc.tl.umap(adata, min_dist=0.3)
    for color in [
        "log1p_total_umis",
        "log1p_n_genes_by_umis",
        "Library",
        "Treatment",
        "Sex",
        "leiden",
        "doublet_score",
        "predicted_doublet",
        "celltype_rough",
    ]:
        sc.pl.umap(adata, color=color, show=False, save=f"_{color}.png")
        sc.pl.pca(adata, color=color, show=False, save=f"_{color}.png")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", pts=True)
    with pd.ExcelWriter(
        f"{figures_embed}/leiden_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["leiden"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
            df.set_index("names")
            df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby="leiden",
        n_genes=int(150 / len(adata.obs["leiden"].unique())),
        show=False,
        save=f"leiden_markers.png",
    )
    sc.pl.pca_overview(adata, color="leiden", show=False, save=True)
    sc.pl.pca_variance_ratio(adata, show=False, save=True)
    sc.pl.pca_loadings(
        adata,
        components=",".join([str(x) for x in range(1, 10)]),
        show=False,
        save=True,
    )
    sc.pl.umap(adata, color=gene_dict["immune"], show=False, save=f"_genes.png")
    with pd.ExcelWriter(
        f"{figures_embed}/metadata_counts.xlsx", engine="xlsxwriter"
    ) as writer:
        obs_list = ["Library", "Treatment", "leiden"]
        num_obs = len(obs_list) + 1
        for ind in range(0, num_obs):
            for subset in itertools.combinations(obs_list, ind):
                if len(subset) != 0:
                    subset = list(subset)
                    if len(subset) == 1:
                        key = subset[0]
                        adata.obs[key].value_counts().to_excel(writer, sheet_name=key)
                    else:
                        key = "_".join(subset)
                        adata.obs.groupby(subset[:-1])[subset[-1]].value_counts(
                            normalize=True
                        ).to_excel(writer, sheet_name=key[:31])
    adata.write(
        f"{output_data}/{adata_name}_filtered_embed.gz.h5ad", compression="gzip"
    )
    adata = sc.read(f"{output_data}/{adata_name}_filtered_embed.gz.h5ad")
