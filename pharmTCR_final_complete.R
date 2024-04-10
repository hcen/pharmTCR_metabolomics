
#install.packages("tidyverse")
library(tidyverse)

#install.packages("readxl")
library(readxl)
library(writexl)

#install.packages("devtools") # most recent version of complexheatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) 

#install.packages("gridtext")
library(gridtext)
library(scales)

#install.packages("ggfortify")
library(ggfortify) # for making PCA

library(ggrepel) # for label text in ggplot2 # https://ggrepel.slowkow.com/articles/examples.html

#install.packages("ggbeeswarm")
library(ggbeeswarm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ##Set working directory to where this file is.
getwd()


# load source file of metabolomics results =======

xcms <- read_excel("input/metabolomics_results_final.xlsx",sheet = 1) 
dim(xcms) #936 detected
dim(xcms[xcms$pvalue<0.05,]) #114 sig
dim(xcms[xcms$pvalue<0.05 & xcms$direction.median_logFC=="up",]) #21 up
dim(xcms[xcms$pvalue<0.05 & xcms$direction.median_logFC=="down",]) #93 down

View(xcms)
colnames(xcms)

maria.new <- read_excel("input/Identified metabolites included in Fig X A_ppm search_HCedits.xlsx",sheet = 1) 
View(maria.new)
colnames(maria.new)
maria.select <- maria.new %>% 
  select(all_of(c("name","HMDB most confident ID_HC",
                "HMDB most confident ID_new", "HMDB Link_new", "ppm" )))

xcms.new <- xcms %>%
  left_join(maria.select, by="name") %>%
  relocate("HMDB most confident ID_HC":"ppm", .before = "mzmed")
View(xcms.new)

xcms.30 <- xcms.new %>%
  mutate(final.label = ifelse(ppm<=30 , `HMDB most confident ID_HC`,NA)) %>%
  relocate(final.label, .before='first.label')
xcms.20 <- xcms.new %>%
  mutate(final.label = ifelse(ppm<=20 , `HMDB most confident ID_HC`,NA)) %>%
  relocate(final.label, .before='first.label')
# 86 to 43 (ppm 20) or 52 (ppm 30)
df <- xcms.30$final.label[!xcms.30$final.label %in% xcms.20$final.label]
View(xcms.30[xcms.30$final.label %in% df,]) 
# between ppm30 and ppm20, none of the additional compounds are in bar graph or volcano (q<0.05), so used ppm20
View(xcms.20)
write_xlsx(xcms.20,"output/xcms_20ppm.xlsx")
xcms.20 <- read_excel("input/xcms_20ppm.xlsx",sheet = 1) 

# heatmaps with compound labels for signigicant ones p<0.05 ==========

m.anno <- xcms.20 %>% filter(pvalue < 0.05) # 114
dim(m.anno)
colnames(m.anno)

m <- as.matrix(as.data.frame(lapply(m.anno[,36:73], as.numeric),check.names=F))

m.z <- t(scale(t(m))) 
View(m.z)
colnames(m.z)

# remove outliers (cutoff is 3 SD) from z score calculation
m.z[m.z>3|m.z<(-3)] <- NA
m.z <- t(scale(t(m.z)))

ceiling(max(abs(m.z)))

# set metabolite labels
xcms.labeled <- xcms.20[!is.na(xcms.20$final.label),]
at <- match(xcms.labeled$name,m.anno$name)
labels <- xcms.labeled$final.label


number_of_pre <- length(grep("Pre", colnames(m.z)))
number_of_post <- length(grep("Post", colnames(m.z)))

end_index_pre <- grep("Pre", colnames(m.z))[number_of_pre]
end_index_post <- grep("Post", colnames(m.z))[number_of_post]

number_of_post
end_index_post

###
heatmap.sig.label <- Heatmap(m.z, #matrix_Z_score_total,
                             name = "Z score",
                             
                             show_column_names = FALSE,
                             show_row_dend = TRUE,
                             
                             show_row_names = FALSE,
                             # show_row_names = T,
                             
                             #    row_labels = gt_render(m.anno$first.label),
                             
                             #row_names_gp = gpar(fontsize = 10),
                             column_names_side = "top",
                             row_dend_side = "left",
                             column_dend_side = "bottom",
                             layer_fun = function(j, i, x, y, width, height, fill) {
                               mat = restore_matrix(j, i, x, y)
                               ind = unique(c(mat[, c(end_index_pre
                               )]))
                               grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                         y = y[ind], 
                                         width = unit(0.05, "inches"), 
                                         height = unit(1/nrow(m.z), "npc"),
                                         gp = gpar(col = "white")
                               )
                             },
                             #col = colorRamp2(c(-ceiling(max(abs(m.z))), 0, ceiling(max(abs(m.z)))), c("blue", "white", "red")),
                             col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
                             
                             top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                                  , height = unit(12, "mm")
                             )),
                             
                             column_order = 1:ncol(m.z),
                             height = 
                               
                               unit(270, 
                                    #60,
                                    "mm"), 
                             
                             width = ncol(m.z)*unit(2, "mm"),
                             border_gp = gpar(col = "black"),
                             show_heatmap_legend = TRUE,
                             heatmap_legend_param = list(
                               title = "Z-score",
                               title_position = "topleft",
                               at = c(-3,-2,-1, 0, 1,2, 3),
                               legend_height = unit(4, "cm"))) +
  rowAnnotation(label = anno_mark(at = at, labels = labels,
                                  labels_gp = gpar(col = "black", fontsize = 11))
  )
draw(heatmap.sig.label)

# save heatmap
png(file = "figures/heatmap_metabolites_p0.05_label_outlierNA_ppm20.png",
    width = 9, 
    height = 12,
    units = "in", res = 600)

draw(heatmap.sig.label)
#
seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

####Condition labels

#Condition label 1
grid.rect(x = (loc2$x - loc1$x)*(end_index_pre)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_pre)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = "lightgrey",
                    col = "lightgrey"
          )
)
grid.text("Pharm-TCR pre", 
          x = (loc2$x - loc1$x)*(end_index_pre)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2
grid.rect(x = (loc2$x - loc1$x)*(end_index_pre + 
                                   end_index_post)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_post)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = "moccasin",
                    col = "moccasin"
          )
)

grid.text("Pharm-TCR post", 
          x = (loc2$x - loc1$x)*(end_index_pre + 
                                   end_index_post)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

## Top label gaps
#Vertical lines
grid.rect(x = end_index_pre/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.03, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

dev.off()


# volcano plot =============
xcms.20$sig.change <- factor(xcms.20$sig.change,
                      levels = c("up","down","ns"))
xcms.20 <- xcms.20 %>%
  mutate(
    volcano_label.q=ifelse(qvalue<0.05, final.label, '')
  )

ggplot(data=xcms.20, aes(x=median_logFC, y=-log10(pvalue), 
                    col=sig.change, 
                    label=volcano_label.q)) +
  scale_color_manual(values=c("red", # nothing positively correlate with leucine or Glu+FA
                              "blue","grey"), #c("#F8766D", "#00BFC4", "grey"),
                     guide = "none"
  ) +
  ylab("-log10 (p value)"#expression(-log[10](p~value))
  ) +
  xlab("log2 (post/pre)"#expression(log[2](fold~change))
  ) +
  geom_point(size=1.2) + 
  
  geom_text_repel(#aes(label=FC),
    size=4,#3.5,
    #color="grey2",
    box.padding   = 0.5,
    point.padding = 0,
    #force=10,
    #force_pull=5,
    max.overlaps = Inf, # always show all label, regardless of overlap
    min.segment.length = 0.01, # always draw line
    segment.color = 'grey70'
  ) +
  
  labs(color="")+
  
  #xlim(-4.9,4.9) +
  
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=20),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "top",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.key=element_rect(fill="white"), # remove grey square around dots
        aspect.ratio = 1/1.2, 
        panel.grid.major = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(filename="figures/volcano_qvalue_median.logFC_ppm20.png",width=6,height=6,units="in",dpi=600)


# bar graphs fold change =================================

xcms.bar <- xcms.20 %>%
  arrange(median_logFC) %>%
  # mutate(order=1:nrow(dat.median)) %>%
  filter(bar.graph=='Yes',ppm<=20) 
View(xcms.bar)
colnames(xcms.bar)

bar <- xcms.bar[,c("final.label","sig.change",
                   "014_GAK_Post_0_Pre_logFC",
                   "016_IAY_Post_0_Pre_logFC",
                   "024_GBK_Post_0_Pre_logFC",
                   "028_MAD_Post_0_Pre_logFC",
                   "032_IAM_Post_0_Pre_logFC",
                   "036_LAA_Post_0_Pre_logFC",
                   "040_GAJ_Post_0_Pre_logFC",
                   "042_IBD_Post_0_Pre_logFC",
                   "044_LAE_Post_0_Pre_logFC",
                   "046_GAU_Post_0_Pre_logFC",
                   "048_GAQ_Post_0_Pre_logFC",
                   "050_IBA_Post_0_Pre_logFC",
                   "056_GAS_Post_0_Pre_logFC",
                   "058_IAZ_Post_0_Pre_logFC",
                   "060_IAC_Post_0_Pre_logFC",
                   "062_FAH_Post_0_Pre_logFC",
                   "066_IBE_Post_0_Pre_logFC",
                   "070_FAD_Post_0_Pre_logFC",
                   "072_GBH_Post_0_Pre_logFC")] %>%
  pivot_longer(cols = -(1:2),names_to = 'subjects',values_to = 'log2FC')
head(bar)

bar$sig.change <- factor(bar$sig.change, levels = c("up","down"))
bar$final.label <- factor(bar$final.label, levels = xcms.bar$final.label)

ggplot(bar, aes(x=final.label, y=log2FC)) +
  geom_hline(yintercept = 0,linetype="dashed", color = "grey60")+
  
  geom_quasirandom(method = "pseudorandom", 
                   aes( colour= sig.change),
                   cex=1.2, width=0.3, alpha=0.5) +
  geom_boxplot(#color='grey40',
    fill=NA,outlier.shape = NA, 
    width=0.3 #, alpha=0.1
  ) +
  coord_flip()+
  
  # set limits or not
  scale_y_continuous(
    #    limits = c(-4,4),
    breaks=c(-8,-6,-4,-2,0,2,4,6)
  ) +
  
  
  #scale_fill_manual(values = c("red","blue") ) +
  scale_color_manual(values = c("red","blue") ) +
  #scale_shape_manual(values=c(23,21,24,22))+ 
  labs(
    y= "log2 (post/pre)" #expression(log[2]~(post/pre))
    #y="Normalized counts"
  )+
  theme_bw()+
  theme(
    axis.title = element_text(#family = "Roboto Light", 
      color="black", size=16),
    axis.title.y = element_blank(),
    axis.text.x = element_text(#family = "Roboto Light",
      color="black",size=12),
    axis.text.y = element_text(#family = "Roboto Light",
      color="black",size=16),
    #axis.ticks = element_line(colour="black"),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(), #,hjust=1,angle = 45
    legend.text = element_text(#family = "Roboto Light",
      color="black",size=16), 
    legend.title = element_blank(),
    legend.position = "top", #c(0.15,0.9)
    strip.text = element_text(#family="Roboto Light",
      color="black",size=14),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    #axis.line = element_line(colour = "black"),
    #aspect.ratio = 1/0.5
  ) + 
  guides(colour = guide_legend(override.aes = list(size=2)))


ggsave("figures/boxplot_metabolites_v3_ppm20_deoxyglucosone.png",width = 8, height = 4, units = "in",dpi=600)


# Prepare a csv for submitting to MetaboAnalyst========

xcms.sub <- xcms %>% select(mzmed, pvalue, tstat, rtmed,log2fold, qvalue)
head(xcms.sub)

colnames(xcms.sub)[c(1:4)] <- c('m.z', 'p.value',	't.score',	'r.t')

xcms.sub <- xcms.sub[order(xcms.sub$t.score, decreasing = T),]

write.csv(xcms.sub[,c(1:4)], file = "output/xcms_metabo.csv")
# note: retention time is in minutes not seconds. Negative Ion mode.


# Organize and plot MetaboAnalyst results ========================

mintg.50.05 <- read_excel("input/MetaboAnalyst_enrichment.xlsx",sheet = 1) 
mgsea.50.05 <- read_excel("input/MetaboAnalyst_enrichment.xlsx",sheet = 2) 
mpath.50.05 <- read_excel("input/MetaboAnalyst_enrichment.xlsx",sheet = 3) 


mpath.50.05 <- read.csv(file="Download_50ppm_p0.05/mummichog_pathway_enrichment.csv")
View(mpath.50.05)
mgsea.50.05 <- read.csv(file="Download_50ppm_p0.05/mummichog_fgsea_pathway_enrichment.csv")
View(mgsea.50.05)
mintg.50.05 <- read.csv(file="Download_50ppm_p0.05/mummichog_integ_pathway_enrichment.csv") # used integrated 2 methods in the plot

colnames(mpath.50.05)
head(mpath.50.05)
head(mgsea.50.05)
mintg.50.05 <- mintg.50.05 %>% left_join(mgsea.50.05[,c(1,5:12)],by="...1") %>% 
  left_join(mpath.50.05[,c(1,4:6)],by="...1") 
View(mintg.50.05)

mintg.sig <- mintg.50.05[mintg.50.05$Combined_Pvals<0.05,] 
View(mintg.sig)
colnames(mintg.sig)

mintg.sig <- mintg.sig %>% mutate(NES.sig=ifelse(P_val<0.05, NES, NA))

# plot pathways

p <- ggplot(mintg.sig, 
            aes(y = fct_reorder(`...1`, -Combined_Pvals))) + 
  geom_point(aes(size = Sig_Hits, color = NES.sig, x=Combined_Pvals))+ 
  scale_x_reverse(limits = c(0.05, 0),breaks=c(0.05,0.03,0.01))+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_gradient2(midpoint = 0, low = "blue3", mid = "white",
                        high = "red3", space = "Lab" )+
  ylab(NULL)+
  xlab("p value")+
  
  theme(#axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black",size=14),
        axis.text.x = element_text(colour = "black",size=12),
        legend.title = element_text(color = "black", size = 14))+
  labs(color = "NES",size="Sig. Hits") #6.75x3.5
p

png(file = "figures/Pathways_Intg_50ppm_p05.png",
     width = 7.5, 
     height = 4, 
     units = "in", res = 600)
p
dev.off()


