# lapply(c(1,2), function(i) { 
#   
#   imap(salmon.quant[names(salmon.quant) %in% c('aale', 'hbec')][[i]], function(x, x_name) {
#     
#     x$deseq$de %>%
#       filter(abs(log2FoldChange) > 1, padj < 0.05, gene %in% msig.df$HALLMARK_INTERFERON_ALPHA_RESPONSE) %>% 
#       mutate(id = x_name)
#     
#   })
# }) %>% purrr::flatten() %>% 
#   bind_rows() %>% 
#   filter(!grepl('compare|hbec_|v_wt', id)) %>% 
#   group_by(gene) %>% 
#   summarize(context = list(id)) %>% 
#   ggplot(aes(context)) + 
#   geom_bar() + 
#   ggupset::scale_x_upset() +
#   geom_text(stat='count', aes(label=after_stat(count)), size = txt.mm_to_pts(0.6))
# 
# msig.df[grepl('INTERFERON|INFLAMMATORY', names(msig.df))] %>% 
#   lapply(function(x) as.data.frame(x)) %>% 
#   bind_rows(.id = 'pathway') %>% 
#   rbind(salmon.quant$aale$intra$deseq$de %>% 
#           filter(log2FoldChange > 1, padj < 0.05) %>% 
#           mutate(pathway = 'kras.intra') %>% select(pathway, 'x' = gene)) %>% 
#   group_by(x) %>% 
#   summarize(pathways = list(pathway)) %>% 
#   ggplot(aes(pathways)) + 
#   geom_bar() + 
#   ggupset::scale_x_upset() +
#   geom_text(stat='count', aes(label=after_stat(count)), vjust = -0.5)
# 
# salmon.quant$aale$intra$deseq$de %>% 
#           filter(padj < 0.05) %>%
#   mutate(pathway = ifelse(gene %in% msig.df$HALLMARK_INTERFERON_ALPHA_RESPONSE, 'IFNA', 'other')) %>% 
#   ggplot(aes(log2FoldChange, -log10(padj), color = pathway)) + 
#   geom_point(alpha = 0.3)
# 
# salmon.quant$aale$intra$deseq$de %>% filter(padj< 0.05) %>% 
#   merge(salmon.quant$aale$exo$deseq$de %>% filter(padj < 0.05), by = 'gene') %>% 
#   mutate(pathway = ifelse(gene %in% msig.df$HALLMARK_INTERFERON_ALPHA_RESPONSE, 'IFNA', 'other')) %>%
#   ggplot(aes(log2FoldChange.x, log2FoldChange.y, color = pathway)) + 
#   geom_point(alpha = 0.3)




# RERWIP
# panel.save <- function(figure, name, width = 'single', height = 'single', plot.in = last_plot()){
#   
#   if (width == 'single'){
#     panel.width = panel.figure.width
#   } else if (width == 'medium'){
#     panel.width = panel.figure.width.medium
#   } else{
#     panel.width = panel.figure.width.double
#   }
#   
#   if (height == 'single'){
#     panel.height = panel.figure.height/3
#   } else if (height == 'double'){
#     panel.height = (panel.figure.height/3)*1.5
#   } else if (height == 'medium'){
#     panel.height = (panel.figure.height/3)*2
#   } else{
#     panel.height = panel.figure.height
#   }
#   
#   if(!dir.exists(paste0(figure.dir, '/','fig.',figure))){
#     
#     dir.create(paste0(figure.dir, '/','fig.',figure))
#     
#   }
#   
#   ggsave(paste0(figure.dir,
#                 '/',
#                 'fig.',
#                 figure,
#                 '/', 
#                 name, 
#                 '.',
#                 panel.width, 
#                 'x', 
#                 panel.height, 
#                 '.pdf'),
#        width = panel.width, 
#        height = panel.height, 
#        units = 'mm', plot = plot.in, device = cairo_pdf)
# }
# 
# save_pheatmap_pdf <- function(x, filename, width=panel.figure.width, height=3) {
#    stopifnot(!missing(x))
#    stopifnot(!missing(filename))
#    pdf(filename, width=unit(width, 'mm'), height=unit(height, 'mm'), family = 'Helvetica')
#    grid::grid.newpage()
#    grid::grid.draw(x$gtable)
#    dev.off()
# }
# 
# calc_ht_size = function(ht, unit = "mm") {
#     pdf(NULL)
#     ht = ComplexHeatmap::draw(ht)
#     w = ComplexHeatmap:::width(ht)
#     w = grid::convertX(w, unit, valueOnly = TRUE)
#     h = ComplexHeatmap:::height(ht)
#     h = grid::convertY(h, unit, valueOnly = TRUE)
#     dev.off()
# 
#     c(w, h)
# }



# ```{r}
# msig.df <- 
#   msigdbr(species = 'Homo sapiens', category = 'H') %>% 
#   dplyr::select(gene_symbol, gs_name) %>% 
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# msig.df$HALLMARK_KRAS_G12D_ZNF_DN <- salmon.quant$aale$intra$deseq$de %>% 
#   filter(gene %in% filter(meta.data$trono.krab.znfs, ZNF == 'KRAB-ZNFs')$gene,
#          log2FoldChange <= -1) %>% 
#   dplyr::select(gene) %>% unlist() %>% unname()
# 
# 
# ```

## bam count promoters NOT WORKING
# ```{r}
# getSeq(BSgenome.Hsapiens.UCSC.hg38, ifn_uniq_peak_ranges) -> ifn_uniq_seqs
# getSeq(BSgenome.Hsapiens.UCSC.hg38, other_uniq_peak_ranges) -> other_uniq_seqs
# memes::runAme(input = ifn_uniq_seqs, control = other_uniq_seqs) -> ifn_uniq_motifs
# 
# salmon.quant$aale$intra$deseq$de %>% 
#   filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
#   select(gene_id) %>% unlist() %>% unname() -> de_list
# 
# de_list_ranges <- gr[names(gr) %in% de_list]
# 
# names(de_list_ranges) <- paste0(seqnames(de_list_ranges)) #, '_', names(de_list_ranges))
# 
# de_list_promoter_seqs <-
#   GenomicFeatures::getPromoterSeq(subject = BSgenome.Hsapiens.UCSC.hg38, 
#                                   query = de_list_ranges, upstream = 1000, downstream = 500)
# 
# lapply(ifn_uniq_motifs$motif_id, function(motif) {
#   
#   print(motif)
#   
#   this_motif <- MotifDb::query(MotifDb::MotifDb, motif) %>% 
#                   universalmotif::convert_motifs()
#   
#   memes::runFimo(de_list_promoter_seqs, motifs = this_motif) -> de_list_promoter_regions
#   
#   bamsignals::bamCount(kras_bam_path,
#                           de_list_promoter_regions) -> kras_de_prom_profile
# 
#   bamsignals::bamCount(ctrl_bam_path,
#                           de_list_promoter_regions) -> ctrl_de_prom_profile
# 
#   lst('kras' = kras_de_prom_profile,
#       'ctrl' = ctrl_de_prom_profile)
#   
# 
# }) -> de_list_fimo_out
# 
# de_list_fimo_out
# ```

# ```{r}
# salmon.quant$aale$intra$deseq$de %>% 
#   filter(padj < 0.05, log2FoldChange > 0.5, gene %in% meta.data$gencode.35.lnc.gene.names$gene) %>% 
#   select(gene_id) %>% unlist() %>% unname() -> de_list
# 
# salmon.quant$aale$intra$deseq$de %>% 
#   filter(padj < 0.05, log2FoldChange < -0.5, gene %in% meta.data$gencode.35.lnc.gene.names$gene) %>% 
#   select(gene_id) %>% unlist() %>% unname() -> de_list_w_hbec
# 
# gr <- rowRanges(salmon.quant$gxi)
# 
# ifn_de_gr <- gr[names(gr) %in% de_list]
# ifn_de_gr_hbec <- gr[names(gr) %in% de_list_w_hbec]
# 
# promoter_seqs <- 
#   GenomicFeatures::getPromoterSeq(subject = BSgenome.Hsapiens.UCSC.hg38, 
#                                   query = ifn_de_gr, upstream = 1000, downstream = 500)
# 
# promoter_seqs_hbec <- 
#   GenomicFeatures::getPromoterSeq(subject = BSgenome.Hsapiens.UCSC.hg38, 
#                                    query = ifn_de_gr_hbec, upstream = 1000, downstream = 500)


# dreme_results <- memes::runDreme(promoter_seqs, control = "shuffle")
# 
# tomTom_res <- memes::runAme(control = 'shuffle',
#                             input = promoter_seqs)


# ```

## Figure S2_
# ```{r}

# salmon.quant$aale$intra$deseq$counts %>%
#   filter(grepl('^IFN', gene),
#          !grepl('P', gene)) %>% 
#   select(gene, contains('.')) %>% 
#   gather(sample, count, -gene) %>% 
#   mutate(count = scale(log1p(count), center = T)) %>%
#   select(sample, gene, count) %>%
#   pivot_wider(names_from = gene, values_from = count) %>%
#   distinct() %>% 
#   column_to_rownames('sample') %>% 
#   t() -> IFN_hm_counts
# 
# IFN_hm_counts %>% 
#   colnames() %>% 
#   as.data.frame() %>% 
#   rename('.' = 'sample') %>% 
#   mutate(condition = ifelse(grepl('ctrl', sample), 'ctrl', 'kras')) %>% 
#   column_to_rownames('sample') %>% 
#   ComplexHeatmap::columnAnnotation(df = . , 
#                                    col = list(condition = c('ctrl' = viridis(3)[1], 
#                                                             'kras' = viridis(3)[2])),
#                                     show_legend = T,
#                                     annotation_legend_param = list(condition = 
#                                                                      list(nrow = 1, title_position = 'lefttop',
#                                                                           title_gp = grid::gpar(fontsize = 7, face = 'bold'))),
#                                     annotation_name_gp = grid::gpar(fontsize = 8),
#                                     simple_anno_size = unit(4.5, "mm")) -> IFN_sample_anno
# 
# IFN_hm_counts %>% 
#   rownames() %>% 
#   as.data.frame() %>% 
#   rename('.' = 'gene') %>% 
#   merge(salmon.quant$aale$intra$deseq$de %>% select(gene, log2FoldChange), by = 'gene', all.x = T) %>% 
#   column_to_rownames('gene') -> IFN_de
#   
# IFN_de[order(match(rownames(IFN_de), rownames(IFN_hm_counts))), , drop =F] -> IFN_de
#   
# toil_de_list <- setNames(IFN_de$log2FoldChange, rownames(IFN_de))
#   
# ComplexHeatmap::rowAnnotation(
#   `AALE log2FC` = ComplexHeatmap::anno_barplot(x = toil_de_list,
#       bar_width = 1, 
#       which = 'row',
#       gp = grid::gpar(col = "black", fill = viridis(3)[2]), 
#       border = FALSE,
#       annotation_name_gp = grid::gpar(fontsize = 8),
#       width = unit(1.25, "cm"),
#       height = unit(0.75, "cm")), 
#   show_annotation_name = TRUE,
#   annotation_name_gp = grid::gpar(fontsize = 8)) -> log2fc_ha
# 
# ComplexHeatmap::Heatmap(IFN_hm_counts, 
#                       name = 'z-score', 
#                       show_row_names = T,
#                       row_dend_side = 'left',
#                       show_column_names = F,
#                       row_split = NULL,
#                       row_title = NULL,
#                       row_title_gp = grid::gpar(fontsize = 7, face = 'bold'),
#                       row_names_side = 'right',
#                       row_dend_width = unit(3, 'mm'),
#                       top_annotation = IFN_sample_anno,
#                       right_annotation = log2fc_ha,
#                       # top_annotation = toil_sample_anno,
#                       heatmap_legend_param = 
#                         list(labels_gp = grid::gpar(fontsize = 7),
#                              legend_gp = grid::gpar(fontsize = 7),
#                              legend_direction = 'vertical',
#                              legend_height = unit(3.5, "cm"),
#                              title_gp = grid::gpar(fontsize = 7, face = 'bold')),
#                       row_names_gp = grid::gpar(fontsize = 7),
#                       heatmap_width = patch_fig_width*0.375) -> IFN_expression_hm.plt
# 
# plot_big_wig_coverage(loci = 'chr19:39,294,877-39,297,773', 
#                       gene_name = 'IFNL1', y_lab = 'scaled read depth') -> plt_S2A.plt
# 
# fig_S2_patch <- "
# AB
# CB
# DD
# "
# 
# 
# 
# wrap_plots(wrap_elements(full = plt_S2A.plt),
#            patchwork::wrap_elements(full = grid::grid.grabExpr(ComplexHeatmap::draw(IFN_expression_hm.plt,
#                                                                          heatmap_legend_side = 'right',
#                                                                          annotation_legend_side = 'top'), 
#                                                                          width = patch_fig_width*0.375)),
#            patchwork::plot_spacer(),
#            patchwork::plot_spacer(),
#            widths = c(1, 1),
#            design = fig_S2_patch) +
#   plot_annotation(tag_levels = 'A') &
#   theme(plot.tag = element_text(size = rel(0.95), face = 'bold'),
#         plot.tag.position = c(0.01, 0.99)) -> figure_S2_patchwork.fig
# 
# 
# ggsave(plot = figure_S2_patchwork.fig,
#        filename = file.path(figure.dir, 'fig.2', 'figure_S2_patchwork.pdf'),
#        width = patch_fig_width, height = patch_fig_height, device = cairo_pdf, units = 'mm')

# ```


# plt_4C_patch <- "
# AAA
# BCD
# EEE
# "
# 
# p4 <- ggplot(data.frame(l = 'DE TF RNA log2FoldChange', x = 1, y = 1)) +
#       geom_text(aes(x, y, label = l), angle = 0, size = txt.mm_to_pts(0.8)) + 
#       theme_void() +
#       theme(plot.margin = margin(t = -6, unit = 'mm')) +
#       coord_cartesian(clip = "off")
# 
# wrap_elements(full = wrap_plots(guide_area(),
#                               intra.aale.te_ranges_list$plts$volcano_plt + theme(legend.direction = 'horizontal') + ggtitle('AALE in.') + xlab(''),
#                               exo.aale.te_ranges_list$plts$volcano_plt + theme(legend.direction = 'horizontal') + ggtitle('AALE ex.') + ylab('') + xlab(''),
#                               hbec.aale.te_ranges_list$plts$volcano_plt + theme(legend.direction = 'horizontal') + ggtitle('HBEC in.') + ylab('') + xlab(''), 
#                               p4,
#                               guides = 'collect', heights = c(0.05, 1, 0.005), design = plt_4C_patch)) -> plt_4C.plt


# base::do.call(rbind, uniq_ctrl_tb) %>%
#   group_by(Ontology) %>%
#   top_n(10, -HyperFdrQ) %>%
#   filter(HyperFdrQ < 0.05,
#          NumFgGenesHit > 1) %>%
#   ggplot(aes(log2(RegionFoldEnrich), reorder(Desc, RegionFoldEnrich), label = Desc, color = -log10(HyperFdrQ))) +
#   geom_point(aes(size = NumFgGenesHit)) +
#   geom_text_repel(aes(label = Desc), 
#                 size = txt.mm_to_pts(0.8), show.legend = F, point.padding = unit(2, 'mm'), color = 'black') + 
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.title = element_text(size = rel(0.7), face = 'plain'),
#         legend.position = 'top', legend.direction = 'horizontal') +
#   scale_x_continuous() +
#   ylab('GO Term') + xlab('log2(FoldEnrichment)') +
#   scale_color_viridis_c() -> plt_6E.plt

# salmon.quant$aale$intra$deseq$de %>% 
#   filter(padj < 0.05, log2FoldChange < -0.5, grepl('^HOX|^KRT', gene) | gene %in% meta.data$trono.krab.znfs$gene) %>% 
#   filter(gene %in% filter(salmon.quant$hbec$lacz_v_v12$deseq$de, padj < 0.05, log2FoldChange > 0)$gene) %>% 
#   select(gene_id) %>% unlist() %>% unname() -> de_list
