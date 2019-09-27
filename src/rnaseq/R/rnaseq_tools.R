# Title     : run_tximport
# Objective : 
# Created by: mdascenzo
# Created on: 2019-05-31

# http://bioconductor.org/packages/release/bioc/html/tximport.html

# ensembldb: an R package to create and use Ensembl-based annotation resources
# https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz031/5301311


# BiocManager::install("tximport")

#/usr/bin/env Rscript
#require(docopt)

library(readr)
library(tximport)
library(rlist)
library(DESeq2)

## tximport.salmon
# example input: 
# sample.ids = c('WB-1044190', 'WB-1044191', 'WB-1044192')
# txt2gene.fp = 'tx2gene.EnsDb.Hsapiens.v86.csv'
# tx_out (default=FALSE): output tx counts instead of gene counts
tximport.salmon = function(salmon_quant_analysis_dir, sample_ids, tx2gene_fp, countsFromAbundance='no', tx_out=FALSE){
  
  # load tx2gene from csv file
  tx2gene = read_csv(tx2gene_fp)
  
  # create file paths to samples
  # possibly this creation step can be omitted if created in pipeline
  files = file.path(salmon_quant_analysis_dir, sample_ids, 'quant.sf')
  names(files) = sample_ids
  
  # basic tximport statement to aggrigate transcript counts at gene
  txi_salmon = tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = countsFromAbundance, txOut = tx_out)
  
  return(txi_salmon)
}

tximport.salmon.alt = function(salmon_dir, tx2gene_fp, countsFromAbundance='no', tx_out=FALSE){
  
  # create txi_fp assuming standard workflow directory structure
  txi_fp = file.path(salmon_dir, 'rnaseq_analysis_txi.rds')
  
  # get sample ids from analysis
  sample_ids = colnames(list.load(txi_fp)$counts)
  
  # load counts for each tx
  txi = tximport.salmon(salmon_dir, sample_ids, tx2gene_fp, countsFromAbundance, tx_out)
  
  return(txi)
}

# save txi to table and list object
write.txi = function(txi, basename, path){
  
  write.table(txi$counts, file=file.path(path, paste(basename, 'counts.txt', sep="_")), sep='\t', quote=FALSE)
  list.save(txi, file=file.path(path, paste(basename, 'txi.rds', sep="_")))
}

# usage:
# sample_table = data.frame(condition =factor(c(rep("N", 10), rep("D", 10))), assay=factor(rep('1', 20)))
load.txi = function(fp, sample_table, filter_min_counts=1){
  
  # load txi data file
  txi = list.load(fp)
  # convert txi into DESeq2 dataset
  dds <- DESeqDataSetFromTximport(txi, sample_table, ~condition)
  # filter by min counts
  keep <- rowSums(counts(dds)) >= filter_min_counts
  dds <- dds[keep,]
  
  vsd <- vst(dds, blind=FALSE)
  
  return(list(txi, dds, vsd))
}



# Title     : mergeFeatureCounts
# Objective : Merge feature count tables from multiple files.
# Created by: jdberndt
# Created on: 2019-09-26

mergeFeatureCounts <- function(files, sample_ids) {
  d <- lapply(files, function(f) read.delim(f, header = T, skip = 1, stringsAsFactors = F))
  ENSG <- gsub("\\.[0-9]{1,2}$", "", d[[1]][,"Geneid"]) #extracts Ensembl geneID without version number
  d <- lapply(d, function(x) x[,7]) #keeps only counts column
  d <- do.call('cbind',d)
  colnames(d) <- sample_ids #renames columns
  rownames(d) <- ENSG
  return(d)
}
