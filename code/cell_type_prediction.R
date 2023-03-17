
cell_types_blood <- list_cell_types()[grep("Blood|Bone_marrow", list_cell_types())]

meth <- readRDS("P:/22_AEL/data/no_snps_methrix_2023.RDS")
res <- housman_cell_type(m = meth, included_cell_types=cell_types_blood)

HSC_mat <- get_matrix(meth[,grep("normal_C010_HSC_HSC", colnames(meth))])
CMP_mat <- get_matrix(meth[,grep("normal_C010_HSC_CMP", colnames(meth))])
MPP_mat <- get_matrix(meth[,grep("normal_C010_HSC_MPP", colnames(meth))])
keep <- complete.cases(HSC_mat) & complete.cases(CMP_mat) & complete.cases(MPP_mat)


meth_hypo <- meth[keep,]
HSC_mat <- HSC_mat[keep,]
CMP_mat <- CMP_mat[keep,]
MPP_mat <- MPP_mat[keep,]
meth_hypo <- meth_hypo[(apply(HSC_mat, 1, function(x) all(x<0.2)) +
                    apply(CMP_mat, 1, function(x) all(x<0.2)) +
                    apply(MPP_mat, 1, function(x) all(x<0.2)))==1,]

library(RnBeads.hg19)
library(RnBeads)
anno <- rnb.get.annotation(type = "probesEPIC")
anno <- unlist(anno)

meth_hypo_f <- subset_methrix(meth_hypo, regions = anno)
meth_hypo_f <- meth_hypo_f[!is.nan(rowMeans2(get_matrix(meth_hypo_f[,-grep("normal", colnames(meth))]), na.rm=T)),]

meth_hypo_f2 <- meth_hypo_f[rowMeans2(get_matrix(meth_hypo_f[,-grep("normal", colnames(meth_hypo_f))]), na.rm=T)>0.6,]
plot <- get_matrix(meth_hypo_f2)[complete.cases(get_matrix(meth_hypo_f2)),]
