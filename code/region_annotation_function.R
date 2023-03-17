##################functions####################
label_func <- function(x){
  breaks <- x
  breaks[breaks==100] <- ">=100"
  breaks
}

region_annotation_function <- function(dmrs, annots_gr = annots_gr_orig, genome="hg19"){

  all_data <- data.frame(data_type=c(rep("Data", length(unique(annots_gr_orig$type))),
                                     rep("Random Regions", length(unique(annots_gr_orig$type)))),
                         annot.type=rep(unique(annots_gr_orig$type), 2))



  #browser()
  plots <- list()

  dmrs <- as.list(dmrs)
  if (length(dmrs)==1){
    names(dmrs) <- "Comparison"
  }

  factor_levels <- gsub(paste0(genome, "_[[:alpha:]]*_"), "", unique(annots_gr$type))
  factor_levels <- factor_levels[order(factor_levels, decreasing = T)]
  cat('\n')

  p_annots_data <- list()
  dmrs_split <- list()
  for (comp in names(dmrs)){
    dmrs_split[[comp]] <- list()
    dmrs_split[[comp]][["up"]] <- makeGRangesFromDataFrame(dmrs[[comp]][dmrs[[comp]]$diff.Methy<0,], keep.extra.columns = T)
    dmrs_split[[comp]][["down"]] <- makeGRangesFromDataFrame(dmrs[[comp]][dmrs[[comp]]$diff.Methy>0,], keep.extra.columns = T)
    p_annots_data[[comp]] <- list()


    for (dataset in c("up", "down")){
      direction_1 <- "Up"
    direction_2 <- "Down"
    direction <- ifelse(dataset=="up", direction_1, direction_2)

    g <- regioneR::getGenomeAndMask(genome)$genome
    g <- regioneR::filterChromosomes(g, chr.type="canonical", organism=genome)


    rnd_annots = annotate_regions(regions = regioneR::randomizeRegions(rep(dmrs_split[[comp]][[dataset]],
                                                                           ifelse(length(dmrs_split[[comp]][[dataset]]) < 500, 10, 1)), genome = g),
                                  annotations = annots_gr,ignore.strand = TRUE)

    result <- annotate_regions(dmrs_split[[comp]][[dataset]], annotations=annots_gr,
                               minoverlap = 10L, ignore.strand = TRUE, quiet = FALSE)


    p_annots = summarize_annotations(annotated_regions = result, annotated_random = rnd_annots)
    p_annots <- merge(p_annots, all_data, by=c("data_type", "annot.type"), all.y=T)
    p_annots$n[is.na(p_annots$n)] <- 0
    p_annots_data[[comp]][[direction]] <- p_annots

    p_annots_data[[comp]][[direction]]$annot.type <- gsub(paste0(genome, "_[[:alnum:]]*_"), "", p_annots_data[[comp]][[direction]]$annot.type)

    p_annots_data[[comp]][[direction]]$annot.type <-  factor(p_annots_data[[comp]][[direction]]$annot.type, levels = factor_levels)

    p_annots_data[[comp]][[direction]]$direction <- direction

  }

    p_annots_data[[comp]] <- dplyr::bind_rows(p_annots_data[[comp]])

  }

  p_annots_data <- dplyr::bind_rows(p_annots_data, .id = "Comparison")

  p_annots_data <- p_annots_data %>% group_by(Comparison, data_type, direction) %>% mutate(percentage=n/sum(n),sum=sum(n))


  p_annots_data$annot.type <- factor(p_annots_data$annot.type, levels = rev(levels(p_annots_data$annot.type)))

  p <- ggplot(p_annots_data[p_annots_data$data_type=="Data",], aes_string(x = "direction", fill = "annot.type",
                                  y = "percentage"))+facet_wrap(~Comparison, ncol = 1)
  p <- p + geom_bar(stat = "identity") + coord_flip() +
    theme_bw()
  p <- p + ylab("Direction") + xlab("Ratio") + ggtitle("Distribution of DMRs")

  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
            "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99",
            "#b15928")
  cols <- colorRampPalette(col3)(length(unique(p_annots_data$annot.type)))
    p <- p + scale_fill_manual(values = rev(cols),
                               guide = guide_legend(reverse = TRUE))
plots[[1]]<- p

 #


#


  final <- p_annots_data %>% group_by(Comparison, direction, annot.type) %>% summarise(ratio = log2(percentage[data_type=="Data"]/percentage[data_type=="Random Regions"]),
                                                                 chi_p = fisher.test(cbind(n, sum))$p.value)

  final$ratio[!is.finite(final$ratio)] <- NaN

  final$logp <- -log10(final$chi_p)
  final$logp[final$logp > 100] <- 100
  final$significant <- ifelse(final$chi_p<0.05, "yes", "no")
  final$significant <- factor(final$significant, levels=c("no", "yes"))

  final$Comparison <- gsub("in", "\n in", final$Comparison)
  #############ggplot_enrichment##################
  enr <- ggplot(data = final, aes(y=annot.type, x=direction))+ facet_wrap(~Comparison, nrow = 1)+
    geom_point(aes(size=logp, fill=ratio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 0, low="#0072B5FF", high="#BC3C29FF", name = "Fold change")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="Log p", labels = label_func) +
    theme_bw()+
    theme(axis.text.x=element_text(size=11, angle = 90), axis.text.y=element_text(size=11),
          strip.text = element_text(size=5), aspect.ratio = 4)+
    ylab("Annotation type")+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=11),  legend.title=element_text(size=11))




  plots[[2]] <- enr

return(plots)
}
