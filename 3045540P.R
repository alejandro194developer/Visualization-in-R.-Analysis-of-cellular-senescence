library(ggplot2)
library(reshape2)
library(dplyr)
library(tibble)
library(ggrepel)
library(amap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(stringr)

source("functions.r")

#get all tables
em = read.table("data\\EM.csv", header=TRUE, row.names = 1, sep="\t")
de_smtd_vs_p = read.table("data\\DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names = 1, sep="\t")
de_s_vs_smtd = read.table("data\\DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names = 1, sep="\t")
de_s_vs_p = read.table("data\\DE_Senes_vs_Prolif.csv", header=TRUE, row.names = 1, sep="\t")
human_background = read.table("data\\Human_Background_GRCh38.p13.csv", header=TRUE, row.names = 1, sep="\t")

em = na.omit(em)
de_smtd_vs_p = na.omit(de_smtd_vs_p)
de_s_vs_smtd = na.omit(de_s_vs_smtd)
de_s_vs_p = na.omit(de_s_vs_p)
human_background = na.omit(human_background)

facet_density_plot(em, "Density plots", "Density", "Expresion (Log10)", 0.1)


master = create_master_table(em, human_background)

# add significance column to de tables
de_s_vs_p = add_significant_flag_to_de_table(de_s_vs_p, 0.05, 1)
de_smtd_vs_p = add_significant_flag_to_de_table(de_smtd_vs_p, 0.05, 1)
de_s_vs_smtd = add_significant_flag_to_de_table(de_s_vs_smtd, 0.05, 1)

# format de table
de_s_vs_p = format_de_tables(de_s_vs_p,master,c(10,16,18,19,20,21))
de_smtd_vs_p = format_de_tables(de_smtd_vs_p,master,c(10,15,18,19,20,21))
de_s_vs_smtd = format_de_tables(de_s_vs_smtd,master,c(10,17,18,19,20,21))

plot_ma(de_s_vs_p, title = "MA plot (S vs P)")
plot_ma(de_smtd_vs_p, title = "MA plot (Smtd vs P)")
plot_ma(de_s_vs_smtd, title = "MA plot (S vs Smtd)")

em_scaled = scale_dataframe(em)

plot_pca(em_scaled,c("prolif","prolif","prolif","senses","senses","senses","senses_mtd","senses_mtd","senses_mtd"), "PCA plot", 1, 2)

#get up regulated genes
de_s_vs_p_up_regulated = get_up_regulated_genes_table_from_de(de_s_vs_p, 0.05, 1)
de_smtd_vs_p_up_regulated = get_up_regulated_genes_table_from_de(de_smtd_vs_p, 0.05, 1)
de_s_vs_smtd_up_regulated = get_up_regulated_genes_table_from_de(de_s_vs_smtd, 0.05, 1)

#get down regulated genes
de_s_vs_p_down_regulated = get_down_regulated_genes_table_from_de(de_s_vs_p, 0.05, 1)
de_smtd_vs_p_down_regulated = get_down_regulated_genes_table_from_de(de_smtd_vs_p, 0.05, 1)
de_s_vs_smtd_down_regulated = get_down_regulated_genes_table_from_de(de_s_vs_smtd, 0.05, 1)

plot_volcano(de_s_vs_p, de_s_vs_p_up_regulated, de_s_vs_p_down_regulated, 0.05, 1, "volacno plot (S vs P)", "log2fold", "p value adjust (log10)")
plot_volcano(de_smtd_vs_p, de_smtd_vs_p_up_regulated, de_smtd_vs_p_down_regulated, 0.05, 1, "volacno plot (Smtd vs P)", "log2fold", "p value adjust (log10)")
plot_volcano(de_s_vs_smtd, de_s_vs_smtd_up_regulated, de_s_vs_smtd_down_regulated, 0.05, 1, "volacno plot (S vs Smtd)", "log2fold", "p value adjust (log10)")

pathway_plot(get_significant_gene_names(de_s_vs_p_up_regulated), "Top 10 pathways (S vs P up regulated)", "Pathway", 'Gene count',"pink","red")
pathway_plot(get_significant_gene_names(de_smtd_vs_p_up_regulated), "Top 10 pathways (Smtd vs P up regulated)", "Pathway", 'Gene count',"pink","red")
pathway_plot(get_significant_gene_names(de_s_vs_smtd_up_regulated), "Top 10 pathways (S vs Smtd up regulated)", "Pathway", 'Gene count',"pink","red")

pathway_plot(get_significant_gene_names(de_s_vs_p_down_regulated), "Top 10 pathways (S vs P down regulated)", "Pathway", 'Gene count',"cyan","blue")
pathway_plot(get_significant_gene_names(de_smtd_vs_p_down_regulated), "Top 10 pathways (Smtd vs P down regulated)", "Pathway", 'Gene count',"cyan","blue")
pathway_plot(get_significant_gene_names(de_s_vs_smtd_down_regulated), "Top 10 pathways (S vs Smtd down regulated)", "Pathway", 'Gene count',"cyan","blue")

#em_symbols pre processing
em_master_symbols = master[,c(1,2,3,4,5,6,7,8,9,10)]
row.names(em_master_symbols) = em_master_symbols[,"SYMBOL"] # genes symbols as rows id
em_master_symbols = em_master_symbols[,-10]
em_master_symbols_scaled = scale_dataframe(em_master_symbols)

em_significant_gene_names = unique(c(get_significant_gene_names(de_s_vs_p),get_significant_gene_names(de_s_vs_smtd),get_significant_gene_names(de_smtd_vs_p)))
em_significant_genes <- subset(em_master_symbols, rownames(em_master_symbols) %in% em_significant_gene_names)
em_significant_genes_scaled <- subset(em_master_symbols_scaled, rownames(em_master_symbols_scaled) %in% em_significant_gene_names)

plot_heatmap(em_significant_genes_scaled, "Heat Map (all significant genes)")

all_de_tables = merge_3_de_tables(de_s_vs_p, de_s_vs_smtd, de_smtd_vs_p, em_significant_gene_names, 1, c("gene_name", "s_vs_p.mean_expression", "s_vs_p.log2fold", "s_vs_p.p", "s_vs_p.p.adj","s_vs_p.significance", "s_vs_smtd.mean_expression", "s_vs_smtd.log2fold", "s_vs_smtd.p", "s_vs_smtd.p.adj","s_vs_smtd.significance", "smtd_vs_p.mean_expression", "smtd_vs_p.log2fold", "smtd_vs_p.p", "smtd_vs_p.p.adj","smtd_vs_p.significance"))


signature1 = row.names(subset(all_de_tables, (all_de_tables$s_vs_p.significance == TRUE & all_de_tables$s_vs_p.log2fold > 0) & (all_de_tables$s_vs_smtd.significance == TRUE & all_de_tables$s_vs_smtd.log2fold < 0)))
em_signature1_scaled <- subset(em_master_symbols_scaled, rownames(em_master_symbols_scaled) %in% signature1)
plot_heatmap(em_signature1_scaled, "Heat Map (Signature 1)")
pathway_plot(signature1, "Top 10 pathways (Signature 1)", "Pathway", 'Gene count',"pink","red")
plot_box_gitter_plot(em_signature1_scaled,"Gene expression (Signature 1)", "Group", "Expression level (scaled)")

signature2 = row.names(subset(all_de_tables, (all_de_tables$s_vs_p.significance == TRUE & all_de_tables$s_vs_p.log2fold < 0) & (all_de_tables$s_vs_smtd.significance == TRUE & all_de_tables$s_vs_smtd.log2fold < 0)))
em_signature2_scaled <- subset(em_master_symbols_scaled, rownames(em_master_symbols_scaled) %in% signature2)
plot_heatmap(em_signature2_scaled, "Heat Map (Signature 2)")
pathway_plot(signature2, "Top 10 pathways (Signature 2) ", "Pathway", 'Gene count',"pink","red")
plot_box_gitter_plot(em_signature2_scaled,"Gene expression (Signature 2)", "Group", "Expression level (scaled)")


signature3 = row.names(subset(all_de_tables, (all_de_tables$s_vs_p.significance == TRUE & all_de_tables$s_vs_p.log2fold <  0) & (all_de_tables$smtd_vs_p.significance == TRUE & all_de_tables$smtd_vs_p.log2fold > 0)))
em_signature3_scaled <- subset(em_master_symbols_scaled, rownames(em_master_symbols_scaled) %in% signature3)
plot_heatmap(em_signature3_scaled, "Heat Map (Signature 3)")
pathway_plot(signature3, "Top 10 pathways (Signature 3) ", "Pathway", 'Gene count',"pink","red")
plot_box_gitter_plot(em_signature3_scaled,"Gene expression (Signature 3)", "Group", "Expression level (scaled)")

