


#' Identify file path for count, and annotations
#'
#' The function identifies the path for the counts or annotation file in
#' the list of files in the provided folder path, if a matching file is not preset
#' or multiple matching files are found the function raises an error
#'
#' @param files_in_folder_path a list of file names
#' @param file_pattern the pattern of the file to search
#' @param file_name the file general name
#' @param folder_path the path
#' @return file_path
#'
#'
#'
test_if_file_can_be_read=function(files_in_folder_path,file_pattern,file_name,folder_path) {
  file_path=files_in_folder_path[grep(file_pattern,files_in_folder_path)]
  if(length(file_path)==0) {
    stop(paste('Could not identify a',file_name,'in',folder_path,'please provide full path to the files or make sure the count file has counts in the file name'))
  }
  if(length(file_path)>1) {
    stop(paste('More then 1',file_name, 'is present in',folder_path,'please provide a full path to the files'))
  }
  file_path=paste(folder_path,file_path,sep = '/')
  return(file_path)
}


#' Read gene expression data
#'
#' This function reads the a count matrix a gene annotation file
#'  and a sample annotation file and returns an expression set.
#'  The function accepts a path to the folder that contains all the files, or a
#'  separate path to each file.
#'  If folder_path is provided the function identifies the files by their names:
#'   counts file should contain count in the name,
#'   gene_annotation_path should contain gene and anno separated by any chr
#'   sample_annotation_path should contain sample and anno separated by any chr
#'
#'
#
#'
#' @param folder_path a path to the folder that contains the files
#' @param count_matrix_path a path to the count matrix
#' @param gene_annotation_path a path to the gene annotation
#' @param sample_annotation_path a path to the sample annotation
#'
#' @return an ExpressionSet
#' @export
#'
#'
read_expression_data=function(folder_path=NULL,count_matrix_path=NULL,gene_annotation_path=NULL,sample_annotation_path=NULL) {
  #if folder path is provided read all files from the folder
  if(!is.null(folder_path)) {
    files_in_folder_path=list.files(path =folder_path )
    count_matrix_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'count',file_name = 'counts_matrix',folder_path)
    gene_ano_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'gene.anno',file_name = 'gene annotation',folder_path)
    sample_ano_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'sample.anno',file_name = 'sample annotation',folder_path)
  } else if (is.null(count_matrix_path) | is.null(gene_annotation_path) | is.null(sample_annotation_path)) {
    stop('Please provide folder_path or all three seperathe file pathes')
  }
  #read using rio to cover bouth txt, and csv files
  count_matrix=rio::import(count_matrix_path)
  gene_ano=rio::import(gene_ano_path)
  sample_ano=rio::import(sample_ano_path)

  #test if first col of count matrix contains gene name
  is_first_col_gene_name=all(sapply(count_matrix[,1],typeof)=='character')
  if(is_first_col_gene_name) {
    gene_name=count_matrix[,1]
    count_matrix = count_matrix[,-1]
    rownames(count_matrix)=gene_name
    count_matrix=as.matrix(count_matrix)
  }
  #test that count matrix has gene name or that gene matrix is the same length as count annotation
  if(! (exists("gene_name") & length(gene_name)==nrow(count_matrix))) {
    stop('Count matrix does not contian gene name, and gene annotation row number differ then count matrix')
  }
  #if nrow gene annotation differ then nrow count matrix fix it
  if(nrow(gene_ano)!=nrow(count_matrix)) {
    #identifiy the gene name type in the count matrix
    col_ind=apply(gene_ano, 2, function(x) (sum(x %in% gene_name))) %>% which.max
    gene_name=data.frame(gene_name)
    colnames(gene_name)=colnames(gene_ano)[col_ind]
    gene_name %<>% dplyr::left_join(gene_ano)
    gene_ano=gene_name
  }
  #if number of samples in annotation differ then number in count matrix riase error
  if(ncol(count_matrix)!=nrow(sample_ano)) {
    stop('N row in sample annotaion differ then ncol of count matrix')
  }

  #create an ExpressionSet
  rownames(gene_ano)=gene_ano$ENSEMBL
  rownames(sample_ano)=sample_ano$sample_id
  expression_set=Biobase::ExpressionSet(assayData = count_matrix,phenoData =Biobase::AnnotatedDataFrame(sample_ano),featureData =Biobase::AnnotatedDataFrame(gene_ano))

  return(expression_set)
}

#' Calculate cpm
#'
#'The function accepts an ExpressionSet array or a count matrix
#'and returns a counts per million matrix
#'
#'
#' @param data
#'
#' @return a cmp matrix
#' @export
#'
#'
get_cpm=function(data) {
  #if data is an ExpressionSet use biobase to get counts
  if(class(data)[1]=='ExpressionSet') {
    count_matrix=Biobase::exprs(exp)
  } else if (class(data)[1]!='matrix') {
    stop("Data should be a count matrix or an ExpressionSet")
    count_matrix=data
  } else {
    count_matrix=data
  }
  # use edgeR to get cmp
  cpm_mat=apply(count_matrix,2,function(x) x/sum(x)*10^6)
  return(cpm_mat)
}

#' Filter out lowly expressed genes
#'
#' The function accepts an ExpressionSet array or a count matrix and keeps
#' genes that express more then cpm_saf (default 1) in at least sampels_per (default 0.9).
#' The function default is to return a filtered cpm_matrix, if return_cpm is F the function reruns the  original data filtered
#'
#'
#' @param data and ExpressionSet array or a count matrix
#' @param cpm_saf the threshold for cmp expression, (default 1)
#' @param sampels_per the threshold for percentage of samples expressing the gene, (default 0.9)
#' @param return_cpm a parameter which determine if to return the cpm or raw data (default F)
#'
#' @return a filtered cpm matrix or the filtered original data
#' @export
#'
filter_lowly_expressed_genes=function(data,cpm_saf=1,sampels_per=0.9,return_cpm=F) {
  cpm_matrix=get_cpm(data)
  #get index for samples to keep
  n_samps=ncol(cpm_matrix)
  ind_keep=apply(cpm_matrix,1,function(x) (sum(x>cpm_saf)/n_samps)>sampels_per)
  if(return_cpm==T) {
    cpm_matrix_f=cpm_matrix[ind_keep,]
    return(cpm_matrix_f)
  }
  data_f=data[ind_keep,]
  return(data_f)
}


#' Get percentage of expression from the top 10 genes
#'
#' This function is a helper function of run_differential_expression
#' The function accepts a cpm_matrix colum and returns the proportion of expression
#' from the top 10 most expressed genes
#'
#' @param expression_vec
#'
#' @return a fraction
#'
top10prop=function(expression_vec) {
  expression_vec_sorted=sort(expression_vec,decreasing = T)
  prop_top_10=sum(expression_vec_sorted[1:10])/sum(expression_vec_sorted)
  return(prop_top_10)
}


#' Run differential expression analysis using edgeR
#'
#' The function accepts an ExpressionSet_data and the name of the group column
#' Lowly expressed genes are filtered out (genes that express less then cpm_saf (default 1) in at least sampels_per (default 0.9) of the smallest group)
#' The function returns differential expression results and if create_plots is True the function also returns diagnostic plots, and go terms enrichmnet
#'
#'
#'
#' @param ExpressionSet_data an ExpressionSet object
#' @param group_col_name a string with the group column name
#' @param cpm_saf the threshold for cmp expression, (default 1)
#' @param sampels_per the threshold for percentage of samples expressing the gene, (default 0.9)
#' @param use_tmm should tmm normalization be used, default is T
#' @param create_plots should plots be created, default is T
#' @param desian an optional design matrix for the model, if none is provided the model is based only on the group
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @return A data frame with gene expression results or a list with the following items
#' dif_exp_results: A data frame with gene expression results
#' norm_factor_diagnostic_plot: A plot of the normalization factor by group, used to determine if to use tmm normalization
#' top_10_exp_plot: A plot of the proportion of expression from the 10 ten most expressed genes, used to determine if to use tmm normalization
#' mds_plot: An mds plot of the results
#' volcano_plot: A volcano plot of the results
#' @export
#'
run_differncial_expression=function(ExpressionSet_data,group_col_name,cpm_saf=1,sampels_per=0.9,use_tmm=T,create_plots=T,design=NULL) {

  if(class(ExpressionSet_data)[1]!='ExpressionSet') {
    stop("Data must be an ExpressionSet data")
  }
  #get count matrix
  count_matrix=Biobase::exprs(ExpressionSet_data)
  #get group vector
  group_vec=ExpressionSet_data[[group_col_name]]
  #get gene annotaton
  gene_ano=ExpressionSet_data@featureData@data
  #create an  edgeR object
  edge_r_counts <- edgeR::DGEList(counts=count_matrix, group=group_vec,genes = gene_ano)

  #remove unexpressed genes
  #get the size of the smallest group
  min_goroup_size=min(table(edge_r_counts$samples$group))
  keep <- rowSums(edgeR::cpm(edge_r_counts)>=cpm_saf) >= min_goroup_size*sampels_per
  print(paste('Retaining',sum(keep),'genes, out of',length(keep),'genes'))
  edge_r_counts=edge_r_counts[keep,]

  #calculate normalization factors using TMM
  if(use_tmm==T) {
    edge_r_counts <- edgeR::calcNormFactors(edge_r_counts)
  }

  #test normalization factors distribution between groups
  if(create_plots==T & use_tmm==T) {
    norm_factor_diagnostic=ggplot2::ggplot(edge_r_counts$samples,ggplot2::aes(x=group,y=norm.factors))+ggplot2::geom_boxplot()+ggplot2::labs(x='',y='Normalization factor')+ggplot2::theme_bw(base_size = 18)
  } else {
    norm_factor_diagnostic=NULL
  }

  #test what proportion of expression arises from the top 10 genes, used to determine the need for TMM normalization
  if(create_plots==T) {
    proortion_of_top_10=apply(edgeR::cpm(edge_r_counts),2, top10prop)
    top_10_exp=ggplot2::qplot(x=edge_r_counts$samples$group,y=proortion_of_top_10,geom='boxplot',xlab = '',ylab ='Proportion of expression\n from top 10 genes')+ggplot2::theme_bw(base_size = 18)
  }
  #Plot mds
  if(create_plots==T) {
    t=limma::plotMDS(edge_r_counts,plot=F)
    mds=t$eigen.vectors
    g_mds=ggplot2::qplot(x=mds[,1],y=mds[,2],col=edge_r_counts$samples$group)+ggplot2::theme_bw(base_size = 18)+ggplot2::labs(col='Group',x='MDS1',y='MDS2')
  }

  #run differential expression
  #if design matrix not provided create one
  if(is.null(design)) {
    design=model.matrix(~edge_r_counts$samples$group)
    colnames(design)[2]=group_col_name
  }

  #estimate dispersion
  edge_r_counts <- edgeR::estimateDisp(edge_r_counts, design, robust=TRUE)

  #fit the model
  fit <- edgeR::glmQLFit(edge_r_counts, design)
  qlf <- edgeR::glmQLFTest(fit, coef=group_col_name)
  dif_exp_results=edgeR::topTags(qlf,n=nrow(edge_r_counts))$table

  #plot volcano
  if(create_plots) {
    #create a sub data for ggrepel
    top10_up=dif_exp_results %>% dplyr::filter(logFC>0 & !is.na(SYMBOL)) %>% dplyr::slice(1:5)
    top10_down=dif_exp_results %>% dplyr::filter(logFC<0 & !is.na(SYMBOL)) %>% dplyr::slice(1:5)
    top=dplyr::bind_rows(top10_down,top10_up)
    #create the plot
    g=ggplot2::ggplot(dif_exp_results,ggplot2::aes(y=-log10(PValue),x=logFC,col=FDR<0.05))+ggplot2::geom_point()+ggrepel::geom_text_repel(data=top,ggplot2::aes(label =SYMBOL),col='blue')+ggplot2::theme_bw(base_size = 18)+ggplot2::scale_color_manual(values = c('red','black'))
    g=g+ggplot2::labs(x='log FC',col='FDR<0.05')
    volcano_plot=g+ggplot2::theme(legend.position="top")
  }


  if(create_plots==T) {
    return(list(dif_exp_results=dif_exp_results,
            norm_factor_diagnostic_plot=norm_factor_diagnostic,top_10_exp_plot=top_10_exp,mds_plot=g_mds,
                volcano_plot=volcano_plot))
  } else {
    return(dif_exp_results)
  }

}

