#************************************************************************************************************************#
############(ERROR) Finding modules of co-regulated genes ############
ciliated_cds_pr_test_res <- graph_test(cds_subTra2, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subTra2[pr_deg_ids,], resolution=1e-2)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subTra2)),
                                cell_group=partitions(cds)[colnames(cds_subTra2)])
agg_mat <- aggregate_gene_expression(cds_subTra2, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)