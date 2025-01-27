'''Goal:Run trajectory inference on DCs including RNA velocity and pseudotime
Date:250127
Author: Carsten Knutsen
conda_env:trajectory_inference
'''


import os
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import palantir
import cellrank as cr
import matplotlib.pyplot as plt
data = 'data/single_cell_files/scanpy_files'
figures = 'data/figures/trajectory_inference/dendritic'
adata_name = "dendritic"

os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, format="png")
sc.settings.figdir = figures
scv.settings.figdir = figures

if __name__ == '__main__':
    # adata_dc = sc.read(f'{data}/{adata_name}_celltyped_just_dcs_no_cc.gz.h5ad')
    # adata_velocity = sc.read(f'{data}/{adata_name}_all_cells_velocyto.gz.h5ad')
    # adata_dc_velocity = adata_velocity[[x for x in adata_dc.obs_names if x in adata_velocity.obs_names], :]
    # adata_dc_velocity = scv.utils.merge(adata_dc, adata_dc_velocity)
    # adata_dc_velocity.X = adata_dc_velocity.layers['raw'].copy()
    # del adata_dc_velocity.uns['log1p']
    # print(adata_dc_velocity)
    # scv.pp.filter_and_normalize(adata_dc_velocity, enforce=True)
    # scv.pp.moments(adata_dc_velocity, n_pcs=30, n_neighbors=30)
    # scv.tl.recover_dynamics(adata_dc_velocity)
    adata_dc_velocity=sc.read(f'{data}/{adata_name}_dc_velocity.gz.h5ad')

    for mode in ['stochastic','steady_state','dynamical']:
        scv.tl.velocity(adata_dc_velocity, mode=mode)
        scv.tl.velocity_graph(adata_dc_velocity)
        scv.pl.velocity_embedding_stream(adata_dc_velocity, basis='umap', density=4, color='Cell Subtype_no_cc', title='',show=False,save=f'{mode}_velocity.png')
    # adata_dc_velocity.X = adata_dc_velocity.layers['log1p'].copy()
    # adata_dc_velocity.write(f'{data}/{adata_name}_dc_velocity.gz.h5ad',compression='gzip')

    # dm_res = palantir.utils.run_diffusion_maps(adata_dc_velocity,
    #                                            n_components=5)  # 5 is the minimum components needed to auto detect root and terminal cells
    # fig = palantir.plot.plot_diffusion_components(adata_dc_velocity)[0]
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_diffusion_components.png')
    # plt.close()
    # ms_data = palantir.utils.determine_multiscale_space(adata_dc_velocity)
    #
    # imputed_X = palantir.utils.run_magic_imputation(adata_dc_velocity)
    #
    # genes = ["Kit", "Kitl", "Tbx2", "Mmp16", 'Gja5']
    #
    # sc.pl.embedding(
    #     adata_dc_velocity,
    #     basis="X_umap",
    #     layer="MAGIC_imputed_data",
    #     color=["Kit", "Kitl", "Tbx2", "Mmp16", 'Gja5'],
    #     frameon=False,
    #     show=False,
    #     save='_palantir_magic_genes.png'
    # )
    #
    # root_cell = palantir.utils.early_cell(adata_dc_velocity, 'Cap1', 'Cell Subtype')
    # terminal_states = palantir.utils.find_terminal_states(adata_dc_velocity, ['Cap2', 'Venous EC', 'Arterial EC'],
    #                                                       'Cell Subtype')
    #
    # fig = palantir.plot.highlight_cells_on_umap(adata_dc_velocity, terminal_states)[0]
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_terminal_cells.png')
    # plt.close()
    #
    # pr_res = palantir.core.run_palantir(
    #     adata_dc_velocity, root_cell, num_waypoints=500, terminal_states=terminal_states
    # )
    #
    # fig = palantir.plot.plot_palantir_results(adata_dc_velocity, s=3)
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_results.png')
    # plt.close()
    #
    # cells = pd.Series(
    #     ["CD31_HyOx_P3_AAACCATTCTAGTGCA-1",
    #      "Bst1_HyOx_P3_TGTGGTTGTTGACCCA-1",
    #      "CD31_HyOx_P3_AAACGGGTCCCCTTGA-1", ],
    #     index=["CD31_HyOx_P3_AAACCATTCTAGTGCA-1",
    #            "Bst1_HyOx_P3_TGTGGTTGTTGACCCA-1",
    #            "CD31_HyOx_P3_AAACGGGTCCCCTTGA-1",
    #
    #            ])
    # fig = palantir.plot.plot_terminal_state_probs(adata_dc_velocity, cells)
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_random_cell_terminal_state_probabilities.png')
    # plt.close()
    # fig = palantir.plot.highlight_cells_on_umap(adata_dc_velocity, cells)[0]
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_random_cell_umap.png')
    # plt.close()
    #
    # masks = palantir.presults.select_branch_cells(adata_dc_velocity, q=.01, eps=.01)
    #
    # fig = palantir.plot.plot_branch_selection(adata_dc_velocity)
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_branch_selection.png')
    # plt.close()
    #
    # for path in ['Cap2', 'Venous EC', 'Arterial EC']:
    #     fig,ax= plt.subplots(1,1)
    #     palantir.plot.plot_trajectory(adata_dc_velocity, path,ax=ax)
    #     fig.tight_layout()
    #     fig.savefig(f'{figures}/palantir_branch_selection_{path}.png')
    #     plt.close()
    # gene_trends = palantir.presults.compute_gene_trends(
    #     adata_dc_velocity,
    #     expression_key="MAGIC_imputed_data",
    # )
    # fig = palantir.plot.plot_gene_trends(adata_dc_velocity, genes + ['Mgp', 'Scn7a', 'Eln', 'Peg3'])
    # fig.tight_layout()
    # fig.savefig(f'{figures}/palantir_gene_trends.png')
    # plt.close()
    #
    #
    # sc.tl.diffmap(adata_dc_velocity)
    #
    # iroot = adata_dc_velocity.obs.index.get_loc(root_cell)
    # scv.pl.scatter(
    #     adata_dc_velocity,
    #     basis="diffmap",
    #     c=["Cell Subtype", iroot],
    #     legend_loc="right",
    #     components=["2, 3"],
    #     show=False,
    #     save='diffmap_Cell Subtype_root_cell.png'
    # )
    #
    # adata_dc_velocity.uns["iroot"] = iroot
    #
    # sc.tl.dpt(adata_dc_velocity)
    # sc.pl.embedding(
    #     adata_dc_velocity,
    #     basis="umap",
    #     color=["dpt_pseudotime", "palantir_pseudotime"],
    #     color_map="viridis",
    #     show=False,
    #     save='_pseudotime.png'
    # )
    #
    # pk = cr.kernels.PseudotimeKernel(adata_dc_velocity, time_key="palantir_pseudotime")
    # pk.compute_transition_matrix()
    # fig = pk.plot_projection(basis="umap", color='Cell Subtype', recompute=True,save=f'{figures}/palantir_pseudotime_stream.png')
    # adata_dc_velocity.write(f'{data}/{adata_name}_vasc_endo_velocity.gz.h5ad',compression='gzip')
    # # fig.savefig(f'{figures}/palantir_pseudotime_stream.png')
    # # plt.close()


