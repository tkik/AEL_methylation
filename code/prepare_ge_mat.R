


mat <- assay(vsd)
mat <- mat[-grep("PAR_Y", rownames(mat)),]

rownames(mat) <- substr(rownames(mat),1,15)

anno_df <- data.frame(
  ensembl = rownames(mat),
  gene_id =   mapIds(org.Hs.eg.db,
                     keys= rownames(mat),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first"),
  gene_name =  mapIds(org.Hs.eg.db,
                      keys= rownames(mat),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first"),
  stringsAsFactors = FALSE,
  row.names =  rownames(mat)
)



mat <- cbind(mat, anno_df)

write.table(mat, file = "data/gene_expression_mat.txt", sep="\t", row.names=F, na="", quote = F)
