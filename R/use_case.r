#use case

#read data
dat=cyto::read_expression_data(folder_path="data/")

#run differential expression
dif_exp_res=cyto::run_differncial_expression(dat,group_col_name = 'type')

#show top diff exp genes
head(dif_exp_res$dif_exp_results)

#plot norm factors vs group
dif_exp_res$norm_factor_diagnostic_plot

#plot proportion of expression from the top 10 genes
dif_exp_res$top_10_exp_plot

#plot mds
dif_exp_res$mds_plot

#plot volcano
dif_exp_res$volcano_plot

#run go term enrichment
library(enrichR)
dbs <- c("GO_Biological_Process_2021")
dif_exp_mat=dif_exp_res$dif_exp_results
up_in_group=dif_exp_mat$SYMBOL[dif_exp_mat$FDR<0.05 & dif_exp_mat$logFC>log2(1.5)]
enriched_up <- enrichR::enrichr(up_in_group, dbs)$GO_Biological_Process_2021
#remove the GO number from term name for easy plotting
enriched_up$short_name=strsplit(enriched_up$Term,'(G',fixed = T) %>% sapply('[',1)
#order the levels according to p.value to get the plot with ordered p.value
enriched_up$short_name=factor(enriched_up$short_name,levels = rev(enriched_up$short_name))
enriched_up_p=ggplot2::ggplot(enriched_up[1:10,],ggplot2::aes(x=short_name,y=-log10(P.value)))+ggplot2::geom_bar(stat="identity")+ggplot2::coord_flip()+ggplot2::theme_bw(base_size = 18)+ggplot2::labs(x='Term',y='-log10(P)',title = 'Enriched in upregulated genes')
enriched_up_p

down_in_group=dif_exp_mat$SYMBOL[dif_exp_mat$FDR<0.05 & dif_exp_mat$logFC<log2(1/1.5)]
enriched_down <- enrichR::enrichr(down_in_group, dbs)$GO_Biological_Process_2021
enriched_down$short_name=strsplit(enriched_down$Term,'(G',fixed = T) %>% sapply('[',1)
enriched_down$short_name=factor(enriched_down$short_name,levels = rev(enriched_down$short_name))
enriched_down_p=ggplot2::ggplot(enriched_down[1:10,],ggplot2::aes(x=short_name,y=-log10(P.value)))+ggplot2::geom_bar(stat="identity")+ggplot2::coord_flip()+ggplot2::theme_bw(base_size = 18)+ggplot2::labs(x='Term',y='-log10(P)',title = 'Enriched in downregulated genes')
enriched_down_p
