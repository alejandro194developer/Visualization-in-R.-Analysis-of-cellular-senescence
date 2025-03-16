create_master_table = function(em_table, annotation_table){
  master = merge(x = em_table, y= annotation_table, by.x = 0, by.y = 0) # join em with annotations 
  row.names(master) = master[,"Row.names"] # genes ids as rows id
  master = master[,-1] # delete row.names column
  master = master[,-14] # delete column biotype (no relevant information)
  master = na.omit(master) # delete na rows
  master$mean_expression = rowMeans(master[,c(1,2,3,4,5,6,7,8,9)]) # add mean expression column
  master$mean_expression_smtd_vs_p = rowMeans(master[,c(1,2,3,7,8,9)])
  master$mean_expression_s_vs_p = rowMeans(master[,c(1,2,3,4,5,6)])
  master$mean_expression_s_vs_smtd = rowMeans(master[,c(4,5,6,7,8,9)])
  return(master)
}

facet_density_plot = function(em, title, x_label, y_label, log_upset){
  master_long <- melt(em)
  master_long$value <- log10(master_long$value + log_upset)
  
  density_plot = ggplot(master_long, aes(x = value)) +
    geom_density(fill = "blue", alpha = 0.5) +
    facet_wrap(~ variable, scales = "free") +  # Separate plots for each column
    labs(title = title,
         x = x_label, y = y_label)
  save_plot("plots_3045540P/",title, density_plot)
  return(density_plot)
}

add_significant_flag_to_de_table = function(de_table, p_treshold, log2_treshold){
  de_table$significance = as.factor(de_table$p.adj < p_treshold & abs(de_table$log2fold) > log2_treshold)
  return(de_table)
}

format_de_tables = function(de_table,master_table, columns_to_select){
  df1 <- master_table %>% rownames_to_column(var = "ID")
  df2 <- de_table %>% rownames_to_column(var = "ID")
  resultado <- left_join(df1, df2, by = "ID")
  resultado <- resultado %>% column_to_rownames(var = "ID")
  resultado <- resultado[,columns_to_select]
  colnames(resultado) <- c("gene_name", "mean_expression", "log2fold", "p", "p.adj","significance")
  sorted_order = order(resultado[,5], decreasing=FALSE)
  resultado = resultado[sorted_order,]
  return(resultado)
}

plot_ma = function(de_table, title){
  # Create an MA plot
  plot_MA = ggplot(de_table, aes(x = log10(mean_expression), y = log2fold, color = significance)) +
    geom_point(alpha = 0.6) +  # Add scatter points with transparency
    scale_color_manual(values = c("grey", "red", "green", "purple")) +  # Customize category colors
    labs(title = title,
         x = "Mean Expression (Log10)",
         y = "Log2 Fold Change") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  save_plot("plots_3045540P/",title, plot_MA)
  return(plot_MA)
}

get_up_regulated_genes_table_from_de = function(de_table, p_threshold, fold_threshold){
  table_sig_up = subset(de_table, de_table$p.adj < p_threshold & de_table$log2fold > (fold_threshold))
  return(table_sig_up)
}

get_down_regulated_genes_table_from_de = function(de_table, p_threshold, fold_threshold){
  table_sig_down = subset(de_table, de_table$p.adj < p_threshold & de_table$log2fold < (-1*fold_threshold))
  return(table_sig_down)
}

plot_pca = function(table_scaled, class_name_vector, title, x_pca, y_pca){
  em_master_matrix = as.matrix(sapply(table_scaled, as.numeric))
  pca = prcomp(t(em_master_matrix))
  pca_coordinates = data.frame(pca$x)
  class_column = class_name_vector
  pca_coordinates$class = class_column
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  pca_plot = ggplot(pca_coordinates, aes(x = pca_coordinates[,x_pca], y = pca_coordinates[,y_pca], color = class)) +
    geom_point(size = 3) +
    labs(title = title,
         x = x_axis_label,
         y = y_axis_label)
  save_plot("plots_3045540P/",title, pca_plot)
  return(pca_plot)
}

plot_volcano = function(de_table, up_regulated_genes_table, down_regulated_genes_table, p_threshold, fold_threshold, title, x_label, y_label){
  table = de_table
  table$mlog10p = -log10(de_table$p.adj)
  up_regulated_genes_table$mlog10p = -log10(up_regulated_genes_table$p.adj)
  down_regulated_genes_table$mlog10p = -log10(down_regulated_genes_table$p.adj)
  
  table_sig_up_top5 = up_regulated_genes_table[1:5,]
  table_sig_down_top5 = down_regulated_genes_table[1:5,]
  
  ggp_volcano_plot = ggplot(table, aes(x = log2fold, y = mlog10p)) + 
    geom_point(colour = "black") +
    geom_point(colour = "red",size= 1, data = up_regulated_genes_table) + 
    geom_point(colour = "blue",size= 1, data = down_regulated_genes_table) + 
    geom_vline(xintercept= -1 * fold_threshold, linetype="dashed", color = "grey", size=0.5) + 
    geom_vline(xintercept= fold_threshold, linetype="dashed", color = "grey", size=0.5) + 
    geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color = "black", size=0.5) + 
    geom_text_repel(data = table_sig_up_top5, aes(label=table_sig_up_top5$gene_name)) +
    geom_text_repel(data = table_sig_down_top5 ,aes(label=table_sig_down_top5$gene_name)) + 
    labs(title = title, x = x_label, y = y_label)
  save_plot("plots_3045540P/",title, ggp_volcano_plot)
  return(ggp_volcano_plot)
}

scale_dataframe = function(em){
  em_scaled = data.frame(t(scale(data.frame(t(em)))))
  return(em_scaled)
}

get_significant_gene_names = function(de_table){
  significant_gene_names = subset(de_table, significance == TRUE)[,1]
  return(significant_gene_names)
}

filter_significant_genes =function(em_table, de_table){
  significant_gene_names = get_significant_gene_names(de_table)
  significant_genes <- subset(em_table, rownames(em_table) %in% significant_gene_names)
  return(significant_genes)
}

plot_heatmap = function(em_table, title){
  # makes a matrix
  hm.matrix = as.matrix(em_table)
  # gets the distances
  y.dist = Dist(hm.matrix, method="spearman")
  # clusters
  y.cluster = hclust(y.dist, method="average")
  # this pulls out the dendrogram
  y.dd = as.dendrogram(y.cluster)
  # this untangles the denrogram
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  # this gets the untangled gene order from the denrogram
  y.order = order.dendrogram(y.dd.reorder)
  # this reorders the original matrix in the new order
  hm.matrix_clustered = hm.matrix[y.order,]
  # makes the colour palette
  colours = c("blue","pink","red")
  palette = colorRampPalette(colours)(100)
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  heat_map = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + scale_fill_gradientn(colours = palette) + 
    labs(title=title, x="", y="") + 
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'))
  save_plot("plots_3045540P/",title, heat_map)
  return(heat_map)
}

pathway_plot = function(de_table_genes_names, title, x_label, y_label, low_color, high_color){
  # converts from ensembl IDs to Entrez - needed for cluster profiler
  sig_genes_entrez = bitr(de_table_genes_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db,readable = T,ont =
                           "CC",pvalueCutoff = 0.05,qvalueCutoff = 0.10)
  
  gene_sets = ora_results$geneID 
  description = ora_results$Description
  p.adj = ora_results$p.adjust
  
  ora_results_table <- data.frame(geneID = gene_sets, description = description, p.adj = p.adj)
  enriched_gene_set = as.character(ora_results_table [1,1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  
  df_transformed <- ora_results_table %>% mutate(gene_count = str_count(geneID, "/") + 1)
  df_transformed <- df_transformed %>% arrange(desc(gene_count))
  

  pathway_plot = ggplot(df_transformed[1:10,], aes(x = reorder(description, gene_count), y = gene_count, fill = gene_count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = title,
         x = x_label,
         y = y_label) +
    scale_fill_gradient(low = low_color, high = high_color) 
  save_plot("plots_3045540P/",title, pathway_plot)
  return(pathway_plot)
}

merge_3_de_tables = function(de_table_1, de_table_2, de_table_3, significant_genes_names, merge_column, column_names_vector){
  merge_df = merge(x=de_table_1, y=de_table_2, by.x = merge_column, by.y = merge_column)
  de_all_merge = merge(x=merge_df, y = de_table_3, by.x = merge_column, by.y = merge_column)
  colnames(de_all_merge) <- column_names_vector
  row.names(de_all_merge) = de_all_merge[,merge_column]
  de_all_merge <- subset(de_all_merge, rownames(de_all_merge) %in% significant_genes_names)
  return(de_all_merge)
}

plot_box_gitter_plot = function(em_table, title, x_label, y_label){
  df_long_signature_1 <- melt(em_table, variable.name = "Group", value.name = "Expression")
  df_long_signature_1$Group <- gsub("_\\d+", "", df_long_signature_1$Group) 
  box_gitter_plot = ggplot(df_long_signature_1, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.6) +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6) +
    labs(title = title, x = x_label, y = y_label)
  save_plot("plots_3045540P/",title, box_gitter_plot)
  return(box_gitter_plot)
  
}

save_plot = function(path, plot_title, plot){
  path <- paste0(path, plot_title,".png") 
  ggsave(path, plot = plot, width = 8, height = 6)
}