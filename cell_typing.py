"""Goal:Adding cell types and reemvedding DCs
Date:250114
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
from functions import plot_obs_abundance,compare_obs_values_within_groups_to_excel

#we set hardcoded paths here
output_data = "data/single_cell_files/scanpy_files"
os.makedirs(output_data, exist_ok=True)
adata_name = "dendritic"
figures = "data/figures/cell_typing"
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures
#usefula gene lists below
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
#leiden dictionary to assign cell types
leiden_ct_dict = {
    "0": "cDC1",
    "1": "Interstitial Mac",
    "2": "Reg. DC",
    "3": "cDC2",
    "4": "Pro. Reg. DC",
    "5": "Plasmacytoid DC",
    "6": "Pro. cDC1-A",
    "7": "Pro. cDC2",
    "8": "mig DC",
    "9": "doublet DC",
    "10": "Pro. cDC1-B",
    "11": "T cell",
    "12": "Alveolar Mac",
    "13": "Monocyte",
    "14": "B cell",
    "15": "doublet DC_T",
    "16": "Neutrophil",
    "17": "Endothelial",
}
leiden_no_cc_ct_dict = {
    "0": "cDC1",
    "1": "Reg. DC",
    "2": "Interstitial Mac",
    "3": "cDC2",
    "4": "Plasmacytoid DC",
    "5": "cDC1",
    "6": "mig DC",
    "7": "Alveolar Mac",
    "8": "Monocyte",
    "9": "T Cell",
    "10": "doublet_DC_Tcell",
    "11": "B cell",
    "12": "Neutrophil",
    "13": "doublet_DC",
    "14": "Endothelial",
}
if __name__ == "__main__":
    figures_all = f"{figures}/all_genes"
    os.makedirs(figures_all, exist_ok=True)
    sc.settings.figdir = figures_all
    adata = sc.read(f"{output_data}/{adata_name}_filtered_embed.gz.h5ad")
    adata.obs['Cell Subtype'] = [leiden_ct_dict[x] for x in adata.obs["leiden"]]
    adata = adata[(~adata.obs['Cell Subtype'].str.startswith('doublet'))&(adata.obs['Cell Subtype']!='Endothelial')]
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
        "Cell Subtype",
    ]:
        sc.pl.umap(adata, color=color, show=False, save=f"_{color}.png")
    sc.pl.umap(adata,color=gene_dict['immune'],show=False,save='_genes.png')
    sc.tl.rank_genes_groups(adata, "Cell Subtype", method="wilcoxon", pts=True)
    with pd.ExcelWriter(
            f"{figures_all}/celltype_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["Cell Subtype"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
            df.set_index("names")
            df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
    adata.write(f"{output_data}/{adata_name}_celltyped.gz.h5ad")


    figures_dc = f"{figures_all}/dc_alone"
    os.makedirs(figures_dc, exist_ok=True)
    sc.settings.figdir = figures_dc
    adata_dc = adata[adata.obs['Cell Subtype'].str.contains('DC')]
    compare_obs_values_within_groups_to_excel(adata_dc, f"Cell Subtype",
                                              output_prefix=f"{figures_dc}/dc_comparisons")
    plot_obs_abundance(adata_dc,f"Cell Subtype",hue="Library",ordered=True,as_percentage=True,save=f'{figures_dc}/celltype_abundance.png',)
    sc.pp.highly_variable_genes(adata_dc, batch_key="Library")
    sc.pp.pca(adata_dc, mask_var="highly_variable")
    sc.pp.neighbors(adata_dc, use_rep="X_pca")
    sc.tl.umap(adata_dc, min_dist=0.3)
    for color in [
        "log1p_total_umis",
        "log1p_n_genes_by_umis",
        "Library",
        "Treatment",
        "Sex",
        "leiden",
        "doublet_score",
        "predicted_doublet",
        "Cell Subtype",
    ]:
        sc.pl.umap(adata_dc, color=color, show=False, save=f"_{color}.png")
    sc.pl.umap(adata_dc,color=gene_dict['immune'],show=False,save='_genes.png')
    sc.tl.rank_genes_groups(adata_dc, "Cell Subtype", method="wilcoxon", pts=True)
    with pd.ExcelWriter(
            f"{figures_dc}/dc_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata_dc.obs["Cell Subtype"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata_dc, key="rank_genes_groups", group=ct)
            df.set_index("names")
            df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
    adata_dc.write(f"{output_data}/{adata_name}_celltyped_just_dcs.gz.h5ad")


    ## Cell typing no cc
    figures_cc = f"{figures}/no_cc"
    os.makedirs(figures_cc, exist_ok=True)
    sc.settings.figdir = figures_cc


    adata = sc.read(f"{output_data}/{adata_name}_filtered_embed.gz.h5ad")
    adata.obs['Cell Subtype'] = [leiden_ct_dict[x] for x in adata.obs["leiden"]]
    cell_cycle_genes = [x.strip() for x in open('data/outside_data/regev_lab_cell_cycle_genes.txt')]
    cell_cycle_genes = [x.lower().capitalize() for x in cell_cycle_genes]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    sc.tl.score_genes(adata, ['Mki67', 'Top2a', 'Birc5', 'Hmgb2', 'Cenpf'], score_name='proliferation_score')
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    adata.layers['no_cc'] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, batch_key="Library")
    sc.pp.pca(adata, mask_var="highly_variable")
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, key_added="leiden_no_cc", resolution=0.5)
    adata.obs['Cell Subtype_no_cc'] = [leiden_no_cc_ct_dict[x] for x in adata.obs["leiden_no_cc"]]
    adata = adata[(~adata.obs['Cell Subtype_no_cc'].str.startswith('doublet'))&(adata.obs['Cell Subtype_no_cc']!='Endothelial')]

    sc.tl.umap(adata, min_dist=0.5)
    for color in [    "Cell Subtype",
                      'leiden',
                      'leiden_no_cc',
                      "Cell Subtype_no_cc",]:
        sc.tl.paga(adata,color)
        sc.pl.paga(adata, save=color,show=False)
    adata.X = adata.layers['log1p'].copy()
    for color in [
        "log1p_total_umis",
        "log1p_n_genes_by_umis",
        "Library",
        "Treatment",
        "Sex",
        "leiden",
        "leiden_no_cc",
        "doublet_score",
        "predicted_doublet",
        "Cell Subtype",
        "Cell Subtype_no_cc",
        'proliferation_score',

    ]:
        sc.pl.umap(adata, color=color, show=False, save=f"_{color}.png")
    sc.pl.umap(adata,color=gene_dict['immune'],show=False,save='_genes.png')
    sc.tl.rank_genes_groups(adata,"leiden_no_cc",method="wilcoxon", pts=True)
    sc.tl.dendrogram(adata,"leiden_no_cc")
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, save='no_cc_leiden', show=False)
    sc.tl.rank_genes_groups(adata, "Cell Subtype_no_cc", method="wilcoxon", pts=True)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,save='no_cc_cts',show=False)
    with pd.ExcelWriter(
            f"{figures_cc}/celltype_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["Cell Subtype_no_cc"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
            df.set_index("names")
            df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
adata.write(f"{output_data}/{adata_name}_celltyped_no_cc.gz.h5ad")

figures_dc = f"{figures_cc}/dc_alone"
os.makedirs(figures_dc, exist_ok=True)
sc.settings.figdir = figures_dc
adata_dc = adata[adata.obs['Cell Subtype_no_cc'].str.contains('DC')]
compare_obs_values_within_groups_to_excel(adata_dc, f"Cell Subtype_no_cc",
                                          output_prefix=f"{figures_dc}/dc_comparisons")
plot_obs_abundance(adata_dc,f"Cell Subtype_no_cc",hue="Library",ordered=True,as_percentage=True,save=f'{figures_dc}/celltype_abundance.png',)
adata_dc.X = adata_dc.layers['no_cc'].copy()
sc.pp.highly_variable_genes(adata_dc, batch_key="Library")
sc.pp.pca(adata_dc, mask_var="highly_variable")
sc.pp.neighbors(adata_dc, use_rep="X_pca")
sc.tl.umap(adata_dc, min_dist=0.3)
for color in ["Cell Subtype",
              'leiden',
              'leiden_no_cc',
              "Cell Subtype_no_cc", ]:
    sc.tl.paga(adata_dc, color)
    sc.pl.paga(adata_dc, save=color, show=False)
adata_dc.X = adata_dc.layers['log1p'].copy()

for color in [
    "log1p_total_umis",
    "log1p_n_genes_by_umis",
    "Library",
    "Treatment",
    "Sex",
    "leiden",
    "leiden_no_cc",
    "doublet_score",
    "predicted_doublet",
    "Cell Subtype",
    "Cell Subtype_no_cc",
    'proliferation_score',

]:
    sc.pl.umap(adata_dc, color=color, show=False, save=f"_{color}.png")
sc.pl.umap(adata_dc,color=gene_dict['immune'],show=False,save='_genes.png')

# sc.tl.rank_genes_groups(adata_dc, "leiden_no_cc", method="wilcoxon", pts=True)
# sc.pl.rank_genes_groups_dotplot(adata_dc, n_genes=10, save='no_cc_leiden',show=False)

sc.tl.rank_genes_groups(adata_dc, "Cell Subtype_no_cc", method="wilcoxon", pts=True)
sc.tl.dendrogram(adata_dc,"Cell Subtype_no_cc")
sc.pl.rank_genes_groups_dotplot(adata_dc, n_genes=10, save='no_cc_cts',show=False)

with pd.ExcelWriter(
        f"{figures_dc}/dc_markers.xlsx", engine="xlsxwriter"
) as writer:
    for ct in adata_dc.obs["Cell Subtype_no_cc"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_dc, key="rank_genes_groups", group=ct)
        df.set_index("names")
        df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
        df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
adata_dc.write(f"{output_data}/{adata_name}_celltyped_just_dcs_no_cc.gz.h5ad")



