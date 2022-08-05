library(biomaRt)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

load(file.path(annot_dir, "Gencode.v19.annotation.RData"))
expr <- read.csv( file.path(input_dir, "EXPR.txt") , stringsAsFactors=FALSE , sep="\t" )

transcript_gene <- data.frame(matrix(ncol = 3, nrow = length(rownames(features_transcript))))
colnames(transcript_gene) <- c('transcript_id', 'transcript_no_ver', 'gene_id')
transcript_gene$transcript_id <- rownames(features_transcript)
transcript_gene$transcript_no_ver <- str_replace(transcript_gene$transcript_id, '\\.\\d', '')
transcript_gene$gene_id <- features_transcript$gene_id

expr <- expr[rownames(expr) %in% transcript_gene$transcript_no_ver, ]
transcript_gene <- transcript_gene[transcript_gene$transcript_no_ver %in% rownames(expr), ]
expr <- expr[order(rownames(expr)), ]
transcript_gene <- transcript_gene[order(transcript_gene$transcript_no_ver), ]

expr_gene <- data.frame(matrix(ncol = length(colnames(expr)), nrow = length(unique(transcript_gene$gene_id))))
colnames(expr_gene) <- colnames(expr)
rownames(expr_gene) <- unique(transcript_gene$gene_id)

for(gene in rownames(expr_gene)){
  transcripts <- transcript_gene[transcript_gene$gene_id == gene, ]$transcript_no_ver
  if(length(transcripts) > 1){
    df <- expr[transcripts, ]
    expr_gene[rownames(expr_gene) == gene, ] <- unlist(lapply(colnames(df), function(col){
      return(sum(df[, col]))
    }))
  }else{
    expr_gene[rownames(expr_gene) == gene, ] <- expr[transcripts, ]
  }
}

genes <- features_gene[rownames(features_gene) %in% rownames(expr_gene), c('start', 'end', 'gene_id')]
size <- genes$end - genes$start
names(size) <- rownames(genes)

# exp = read.csv( file.path(input_dir, "EXPR.txt") , stringsAsFactors=FALSE , sep="\t" )
# exp$patient_id = paste( "P" , exp$patient_id , sep="" )
# 
# gene = sort( unique(  exp$gene_name ) )
# patient = sort( unique(  exp$patient_id ) )
# 
# expr =  matrix( nrow=length(gene) , ncol= length(patient) )
# colnames(expr) = patient
# rownames(expr) = gene
# 
# for( i in 1:length(patient) ){
# 	p = exp[ exp$patient_id %in% patient[i] , ]
# 	expr[ exp$gene_name , patient[i] ] = p$est_counts
# }
# 
# 
# ##################################
# ## Get Gene Length
# 
# genes <- rownames(expr)
# human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# 
# gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol",values=genes, mart=human)
# 
# size=gene_coords$end_position - gene_coords$start_position
# names(size) = gene_coords[,"hgnc_symbol"]
# 
# gene = intersect( names(size) , rownames(expr) )
# 
# size = size[ gene ]
# expr = expr[ gene , ]

##################################
##################################
## Compute TPM data

GetTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

tpm = log2( GetTPM(expr_gene,size) + 0.001 )


##################################
##################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( tpm , file=file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
