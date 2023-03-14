


#' Identify file path for cout, and annotations
#'
#' The function identifies the path for the counts or annotation file in
#' the list of files in the provided folder path, if a matching file is not preset
#' or multiple matching files are found the function raises an error
#'
#' @param files_in_folder_path a list of file names
#' @param file_pattern the pattern of the file to search
#' @param file_name the file general name
#'
#' @return
#'
#'
#'
test_if_file_can_be_read=function(files_in_folder_path,file_pattern,file_name) {
  file_path=files_in_folder_path[grep(file_pattern,files_in_folder_path)]
  if(length(file_path)==0) {
    stop(paste('Could not identify a',file_name,'in',folder_path,'please provide full path to the files or make sure the count file has counts in the file name'))
  }
  if(length(file_path)>1) {
    stop(paste('More then 1',file_name, 'is present in',folder_path,'please provide a full path to the files'))
  }
  file_path=paste0(folder_path,file_path)
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
#' @param count_matrix_path
#' @param gene_annotation_path
#' @param sample_annotation_path
#'
#' @return
#' @export
#'
#'
read_expression_data=function(folder_path=NULL,count_matrix_path=NULL,gene_annotation_path=NULL,sample_annotation_path=NULL) {
  #if folder path is provided read all files from the folder
  if(!is.null(folder_path)) {
    files_in_folder_path=list.files(path =folder_path )
    count_matrix_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'count',file_name = 'counts_matrix')
    gene_ano_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'gene.anno',file_name = 'gene annotation')
    sample_ano_path=test_if_file_can_be_read(files_in_folder_path,file_pattern = 'sample.anno',file_name = 'sample annotation')
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
#' @return
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
#' @param return_cpm a parameter which determine if to return the cpm or raw data (default T)
#'
#' @return
#' @export
#'
filter_lowly_expressed_genes=function(data,cpm_saf=1,sampels_per=0.9,return_cpm=T) {
  cpm_matrix=get_cpm(data)
  #get index for samples to keep
  n_samps=ncol(cpm_matrix)
  ind_keep=apply(cpm_matrix,1,function(x) (sum(x>cpm_saf)/n_samps)>sampels_per)
  if(return_cpm==T) {
    cpm_matrix_f=cpm_matrix[ind,]
    return(cpm_matrix_f)
  }
  data_f=data[ind_keep,]
  return(data_f)
}





