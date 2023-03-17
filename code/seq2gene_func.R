seq2gene_func <- function (seq, tssRegion, flankDistance, TxDb, sameStrand = FALSE)
{
  ChIPseekerEnv <- get("ChIPseekerEnv", envir = .GlobalEnv)
  if (exists("exonList", envir = ChIPseekerEnv, inherits = FALSE)) {
    exonList <- get("exonList", envir = ChIPseekerEnv)
  }
  else {
    exonList <- exonsBy(TxDb)
    assign("exonList", exonList, envir = ChIPseekerEnv)
  }
  exons <- ChIPseeker:::getGenomicAnnotation.internal(seq, exonList, type = "Exon",
                                         sameStrand = sameStrand)
  if (exists("intronList", envir = ChIPseekerEnv, inherits = FALSE)) {
    intronList <- get("intronList", envir = ChIPseekerEnv)
  }
  else {
    intronList <- intronsByTranscript(TxDb)
    assign("intronList", intronList, envir = ChIPseekerEnv)
  }
  introns <- ChIPseeker:::getGenomicAnnotation.internal(seq, intronList,
                                           type = "Intron", sameStrand = sameStrand)

  if (all(!is.na(exons))){
  genes <- c(exons$gene, introns$gene)} else {
    genes <- c(introns$gene)
  }

  genes <- gsub("\\w+\\.*\\d*/(\\d+)", "\\1", genes)
  features <- ChIPseeker:::getGene(TxDb, by = "gene")
  idx.dist <- ChIPseeker:::getNearestFeatureIndicesAndDistances(seq, features,
                                                   sameStrand = sameStrand)
  nearestFeatures <- features[idx.dist$index]
  distance <- idx.dist$distance
  pi <- distance > tssRegion[1] & distance < tssRegion[2]
  promoters <- mcols(nearestFeatures[pi])[["gene_id"]]
  nearest_genes <- mcols(nearestFeatures[!pi][abs(distance[!pi]) <
                                                flankDistance])[["gene_id"]]
  genes <- c(genes, promoters, nearest_genes)
  return(unique(genes))
}
