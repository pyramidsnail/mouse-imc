setwd("~/projects/WC")
options(stringsAsFactors = F)
library(imcRtools)
library(paletteer)
library(ggplot2)
library(tidyverse)
library(data.table)
load("speWithUndefined.RData")
load("sub_spes.Rdata")
##### check for information of ZXJ
info <- read.csv("info(2).csv", header = T)
markers <- read.csv("data/marker.csv")


spe <- spe[,spe$Annotation!="Undefined"]

spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster3")] <- "ProliferatingMφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster5")] <- "LAMP3+DC"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster7")] <- "PD-L1+Mφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster9")] <- "AntigenPresentingMφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster11")] <- "AntigenPresentingMφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster15")] <- "InfiltratingMφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster1","Cluster2","Cluster6","Cluster14")] <-
  "C1QC+Mφ" 
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster4","Cluster8")] <- "CD163+Mφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster12")] <- "SPP1+IDO+Mφ"
spe.mac$sub.anno[spe.mac$sub.cluster %in% c("Cluster10","Cluster13")] <- "SPP1+Mφ"
table(spe.mac$sub.anno)

spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster15","Cluster23","Cluster18")] <- "CD4T&B"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster1","Cluster4")] <- "GZMB+CD8Teff"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster2")] <- "CD45RO+ActivatedCD8T"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster19")] <- "CD38+CD8T"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster3")] <- "CD38+Treg"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster20","Cluster11","Cluster16")] <- "NK"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster12")] <- "ProliferatingB"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster14","Cluster17","Cluster22")] <- "B"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster24")] <- "CXCL13+B"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster6","Cluster8")] <- "Treg"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster5")] <- "CD8T&B"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster9")] <- "CD4&CD8T"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster10")] <- "CD45RO+CD4T"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster7")] <- "NK&CD8T"
spe.lym$sub.anno[spe.lym$sub.cluster %in% c("Cluster13","Cluster21")] <- "outlier"
# spe.lym <- spe.lym[,spe.lym$sub.anno!="outlier"]

spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster16","Cluster20","Cluster21","Cluster24","Cluster25","Cluster11")] <- 
                   "ProliferatingEpithelial"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster22","Cluster19","Cluster15","Cluster18")] <- "EMT"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster2","Cluster5")] <- "IDO1+Epithelial"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster9","Cluster26","Cluster4")] <- "PD-L1+Epithelial"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster23","Cluster1","Cluster14")] <- "Epithelial"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster17","Cluster13","Cluster3","Cluster6","Cluster8")] <- 
  "IDO1+S100A9+Epithelial"
spe.epi$sub.anno[spe.epi$sub.cluster %in% 
                   c("Cluster12","Cluster10","Cluster7")] <- "S100A9+Epithelial"


spe$sub.anno <- spe$Annotation
spe$id <- paste0(spe$ImageNb,"-",spe$z)
spe$sub.anno[match(paste0(spe.epi$ImageNb,"-",spe.epi$z), spe$id)] <- 
  spe.epi$sub.anno
spe$sub.anno[match(paste0(spe.lym$ImageNb,"-",spe.lym$z), spe$id)] <- 
  spe.lym$sub.anno
spe$sub.anno[match(paste0(spe.mac$ImageNb,"-",spe.mac$z), spe$id)] <- 
  spe.mac$sub.anno
spe <- spe[,spe$sub.anno!="outlier"]
spe.lym <- spe.lym[,spe.lym$sub.anno!="outlier"]

main.markers <- markers[markers$lineage==1,"Antibody"]
epi.markers <- c("S100A9","HLAI","E-cadherin","IDO1","PD-L1","Ki-67",
                 "Pan-cytokeratin","HIF1a")
lym.markers <- c("FoxP3","CD4","CD8","CXCL13","CD127","CD45","CD38","CCR7",
                 "CD20","GranzymeB",
                 "PD-1","Ki-67","CD45RA","CD3","CD45RO","CD57","HIF1a")
mac.markers <- c("S100A9","CD16","CD14","CD45","IDO1","PD-L1","Ki-67",
                   "CD68","HIF1a","CD11b","CD11c","SPP1","HLA-DR",
                   "CD208","C1QC","CD163")

##### cluster heatmaps
library(pheatmap)
main.order <- c("Epithelial","Endothelial","CollagenI+Fibroblast","FAP+Fibroblast",
                "αSMA+Fibroblast","Granulocyte","Lymphocytes","Mono&Mac","Lineage-")
main.markers.order <- c("Pan-cytokeratin","E-cadherin","CD31","CollagenI","FAP","aSMA",
                        "CD15","CD45","CD3","CD4","CD8","CD57","CD20","CD68","CD14","CD16",
                        "CD163","CD11b","CD11c")
epi.order <- c("EMT","Epithelial","ProliferatingEpithelial","IDO1+S100A9+Epithelial",
               "IDO1+Epithelial","S100A9+Epithelial","PD-L1+Epithelial")
epi.markers.order <- c("Pan-cytokeratin","E-cadherin","Ki-67","IDO1","S100A9","HLAI",
                       "PD-L1","HIF1a")
lym.order <- c("CD45RO+CD4T","Treg","CD38+Treg","CD4T&B","CD4&CD8T",
               "GZMB+CD8Teff","CD45RO+ActivatedCD8T","CD38+CD8T","CD8T&B",
               "NK&CD8T","NK","B","CXCL13+B","ProliferatingB")
lym.markers.order <- c("CD45","CD45RA","CD45RO","CD3","CD4","CD8","FoxP3",
                       "CCR7","CD127","CD38","CXCL13","GranzymeB","PD-1","HIF1a",
                       "CD57","CD20",
                       "Ki-67")

breaksList = seq(0, 1, by = 0.01)

exprs.m <- aggregate(t(assay(spe, "exprs")[main.markers,]),
                     by=list(spe$Annotation), 
                     mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
grDevices::cairo_pdf("figures/heatmap_all.pdf", height = 8, width = 12)
pheatmap(exprs.m[main.order, main.markers.order], border_color = "white", scale = "none", 
         cellwidth = 30,
         cellheight = 30*.618,fontsize=15,
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
dev.off() 

exprs.m <- aggregate(t(assay(spe.epi, "exprs")[epi.markers,]),
                     by=list(spe.epi$sub.anno), 
                     mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
grDevices::cairo_pdf("figures/heatmap_epi.pdf", height = 8, width = 12)
pheatmap(exprs.m[epi.order, epi.markers.order], border_color = "white", scale = "none", 
         cellwidth = 30,
         cellheight = 30*.618,fontsize=15,
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
dev.off() 

exprs.m <- aggregate(t(assay(spe.lym, "exprs")[lym.markers,]),
                     by=list(spe.lym$sub.anno), 
                     mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
grDevices::cairo_pdf("figures/heatmap_lym.pdf", height = 8, width = 12)
pheatmap(exprs.m[lym.order, lym.markers.order], border_color = "white", scale = "none", 
         cellwidth = 30,
         cellheight = 30*.618,fontsize=15,
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
dev.off() 

exprs.m <- aggregate(t(assay(spe.mac, "exprs")[mac.markers,]),
                     by=list(spe.mac$sub.anno), 
                     mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
mac.order <- c("CD163+Mφ","InfiltratingMφ","C1QC+Mφ","PD-L1+Mφ",
               "AntigenPresentingMφ","SPP1+Mφ","SPP1+IDO+Mφ",
               "LAMP3+DC","ProliferatingMφ")
mac.markers.order <- c("CD45","CD68","CD16","CD14","CD163","CD11b","CD11c",
                       "SPP1","IDO1","S100A9","HIF1a","HLA-DR",
                       "C1QC","PD-L1","CD208","Ki-67")
grDevices::cairo_pdf("figures/heatmap_mac.pdf", height = 8, width = 12)
pheatmap(exprs.m[mac.order, mac.markers.order], border_color = "white", scale = "none", 
         cellwidth = 30,
         cellheight = 30*.618,fontsize=15,
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
dev.off() 
save(spe, spe.epi, spe.lym, spe.mac, file="spes.RData")

##### composition analysis, assign over-all order and colors 
stats.ct <- colData(spe) %>% as_tibble() %>% count(sub.anno)
stats.ct$frequncy <- stats.ct$n/sum(stats.ct$n)

main.cols <- setNames(c("#b7ded2","#e56b6f","#ffecb8","#665687"),c("Epithelial","Stromal","Lymphocyte",
                                                                   "Myloid"))
ct.order <- c("ProliferatingEpithelial","EMT","IDO1+S100A9+Epithelial","PD-L1+Epithelial",
              "Epithelial","S100A9+Epithelial","IDO1+Epithelial",
              "Endothelial","αSMA+Fibroblast","FAP+Fibroblast","CollagenI+Fibroblast",
              "B","CD4T&B","Treg","ProliferatingB","CD38+Treg","CD45RO+CD4T","GZMB+CD8Teff","NK",
              "CD45RO+ActivatedCD8T","CD8T&B","NK&CD8T","CD4&CD8T","CD38+CD8T","CXCL13+B",
              "Granulocyte","C1QC+Mφ", "AntigenPresentingMφ","SPP1+Mφ","CD163+Mφ","ProliferatingMφ",
              "SPP1+IDO+Mφ","InfiltratingMφ","PD-L1+Mφ","LAMP3+DC")
ct.cols <- setNames(c("#809bce",  "#84a6d0", "#8db1d1", "#91bdd3","#9acbd5","#acdad3", "#b8e0d2",
                      "#9D2A36","#C44E54","#E56B6F","#FFCFCE",
                      "#ffb59f", "#ffbc9e", "#ffc39e", "#ffcb9f", "#ffd2a1", "#ffd7a4",
                      "#ffdba7", "#ffe0aa", "#ffe3ad", "#ffe6b1", "#ffe9b4", "#ffecb8",
                      "#F2FEDC", "#FFF9E8",
                      "#665687", "#705f91","#79689a", "#8d7bae","#9987b8", 
                      "#a693c2", "#b2a0cc", "#d6c5e5", "#f8ebff","#f7f3ff"),ct.order)
     

colData(spe) <- left_join(colData(spe) %>% as.data.frame(), info, 
                          by=c("ImageNb"="roi_id")) %>% DataFrame()
stats <- colData(spe) %>% as_tibble() %>% group_by(ImageNb) %>% count(sub.anno)
stats.total <- colData(spe) %>% as_tibble() %>% count(ImageNb) 
stats <- left_join(expand_grid(ImageNb = unique(stats$ImageNb), 
                               sub.anno=unique(stats$sub.anno)),stats, 
                   by=c("ImageNb","sub.anno")) %>% replace(is.na(.), 0)
stats$total <- stats.total$n[match(stats$ImageNb, stats.total$ImageNb)]
stats$frequency <- stats$n/stats$total
write.csv(stats[,c(1,2,5)],"ml/ct_freq.csv", col.names = T,
          row.names = F, quote = F)
stats <- left_join(stats, info, by=c("ImageNb"="roi_id"))

ttest.res<- Reduce(rbind, list(do.call(rbind,lapply(split(stats,factor(stats$sub.anno)), function(x){
  ttest <- t.test(x[x$group1=="pre(0-1)","frequency"],
                  x[x$group1=="pre(2-3)","frequency"])
  return(c(unique(x$sub.anno), ttest$estimate, ttest$p.value,"pre(0-1)\nv.s.\npre(2-3)"))
})) %>% as.data.frame(),
do.call(rbind,lapply(split(stats,factor(stats$sub.anno)), function(x){
  ttest <- t.test(x[x$group2=="post(0-1)","frequency"],
                  x[x$group2=="post(2-3)","frequency"])
  return(c(unique(x$sub.anno), ttest$estimate, ttest$p.value,"post(0-1)\nv.s.\npost(2-3)"))
})) %>% as.data.frame(),
do.call(rbind,lapply(split(stats,factor(stats$sub.anno)), function(x){
  ttest <- t.test(x[x$group3=="post(0-1)","frequency"],
                  x[x$group3=="pre(0-1)","frequency"])
  return(c(unique(x$sub.anno), ttest$estimate, ttest$p.value,"post(0-1)\nv.s.\npre(0-1)"))
})) %>% as.data.frame(),
do.call(rbind,lapply(split(stats,factor(stats$sub.anno)), function(x){
  ttest <- t.test(x[x$group4=="post(2-3)","frequency"],
                  x[x$group4=="pre(2-3)","frequency"])
  return(c(unique(x$sub.anno), ttest$estimate, ttest$p.value,"post(2-3)\nv.s.\npre(2-3)"))
})) %>% as.data.frame()
))
colnames(ttest.res)[2:3] <- c("x","y")
ttest.res <- ttest.res[ttest.res$V1!="Lineage-",]
ttest.res <- ttest.res %>% mutate(x=as.numeric(x)) %>% 
  mutate(y=as.numeric(y)) %>% 
  mutate(V4=as.numeric(V4)) %>% 
  mutate(fc = x/y) %>%
  mutate(V5 = factor(V5, levels = c("pre(0-1)\nv.s.\npre(2-3)","post(0-1)\nv.s.\npost(2-3)",
                                    "post(0-1)\nv.s.\npre(0-1)","post(2-3)\nv.s.\npre(2-3)"))) %>%
  mutate(V1=factor(V1, levels = rev(ct.order)))
range(log2(ttest.res$fc))

# library(hrbrthemes)
p <- ggplot(ttest.res, aes(x=V5,y=V1))+
  geom_point(aes(color=log2(fc), size=abs(log2(fc))>1 & as.numeric(V4)<.05), shape=15)+
  scale_size_manual(values = c(.8,3.5),name="|log2(FoldChange)|>1&\np-value<0.05")+
  scale_color_gradientn(colours=as.character(rev(paletteer_d("RColorBrewer::RdYlBu"))),
                        values =c(0, seq(.2,.8, length.out=9), 1),
                        name="log2(FoldChange)", limits=c(-5, 5))+
  theme_classic() +xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90,hjust=1, vjust = 0.5),legend.margin=margin(t=-5))+
  coord_fixed(ratio = 1)+
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=10))
p
grDevices::cairo_pdf("figures/cluster_ratio_compare.pdf", width = 6, height = 6)
print(p)
dev.off()

stats.ct <- stats.ct %>% 
  mutate(sub.anno=factor(sub.anno,levels = ct.order))
stats.ct <- stats.ct[complete.cases(stats.ct),]
p <- ggplot(stats.ct, aes(x=sub.anno,y=frequncy))+geom_bar(aes(fill=sub.anno), stat="identity")+
  scale_fill_manual(values = ct.cols)+ylab("Frequency")+xlab("")+
  theme_classic()+theme(axis.text.x = element_blank(), legend.position = "none")
grDevices::cairo_pdf("figures/cluster_ratio_bar.pdf", width = 6, height = 2.5)
print(p)
dev.off()
spe$meta <- spe$Annotation
spe$meta[spe$Annotation %in% c("FAP+Fibroblast","αSMA+Fibroblast","CollagenI+Fibroblast","Endothelial")] <- 
  "Stromal"
spe$meta[spe$Annotation %in% c("Granulocyte","Mono&Mac")] <- "Myloid"
save(spe, file="speWithAnnotation.RData")

##### region analysis
# spe.test <- spe[, colData(spe)$ImageNb=="T001_ROI001"]
source("utils.R")

# spe.test <- spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")]
# spe.test<- buildSpatialGraph(spe.test, img_id = "ImageNb",type = "expansion", 
#                         threshold = 30)
# 
# plotSpatial(spe.test,
#             img_id = "ImageNb",
#             node_color_by = "meta",
#             node_shape_by = "ImageNb")
#             # draw_edges = TRUE,
#             # colPairName = "expansion_interaction_graph",
#             # directed = FALSE,
#             # scales = "free")
# immune_region <- spe.test$sub.anno %in% c(lym.order, mac.order,"Granulocyte")
# spe.test <- patchDetection(spe.test, patch_cells = immune_region,
#                      colPairName = "expansion_interaction_graph",
#                      expand_by = 0.1, min_patch_size = 20,
#                      convex=F,
#                      img_id = "ImageNb")
# plotSpatial(spe.test, node_color_by="patch_id",img_id = "ImageNb")
# epi_region <- spe.test$sub.anno %in% epi.order
# spe.test <- patchDetection(spe.test, patch_cells = epi_region,
#                            colPairName = "expansion_interaction_graph",
#                            expand_by = 0.1, min_patch_size = 20,
#                            convex=F,
#                            img_id = "ImageNb")
# plotSpatial(spe.test, node_color_by="patch_id",img_id = "ImageNb")
# spe.test$label <- immune_region
# plotSpatial(spe.test, node_color_by="label",img_id = "ImageNb")
# 
# 
# node_id <- seq_len(sum(immune_region))
# cur_graph <- graph_from_data_frame(colPair(spe[,epi_region], "expansion_interaction_graph"),
#                                    vertices = data.frame(row.names = node_id, id = node_id))
# cur_components <- components(cur_graph)
# cur_clusters <- cur_components$membership
# cur_clusters[!(cur_clusters %in%
#                  which(cur_components$csize >= min_patch_size))] <- NA
# 
# cur_out <- vector(mode = "character", length = ncol(dat))
# cur_out[!patch_cells] <- NA
# cur_out[patch_cells] <- cur_clusters
# colData(dat)[[name]] <- cur_out
# plotSpatial(dat,img_id = "ImageNb",node_color_by = "region",flip_y=F)
# 
# cur_coords <- spatialCoords(cur_obj)[,c("Pos_X", "Pos_Y")]
# cells <- st_multipoint(as.matrix(cur_coords))
# cells_sfc <- st_cast(st_sfc(cells), "POINT")
# plot(cells_sfc)
# 
# x <- cbind(colData(dat), cur_coords) %>% as_tibble %>%
#   filter(!is.na(!!sym("region"))) %>%
#   nest_by(!!sym("region"))
# x <- x[1,2][[1]][[1]]
# cur_coords <- as.matrix(cbind(x[["Pos_X"]], x[["Pos_Y"]]))
# hull <- data.frame(concaveman(cur_coords, concavity = 2.3))
# polygon <- st_polygon(list(as.matrix(hull)))
# plot(polygon)
# 
# test <- st_intersects(polygon, cells_sfc)
# plot(cells_sfc[test[[1]]])
# 
# 
# 

spe<- buildSpatialGraph(spe, img_id = "ImageNb",type = "expansion", 
                        threshold = 30)
stromal_region <- spe$sub.anno %in% c("FAP+Fibroblast","αSMA+Fibroblast","CollagenI+Fibroblast","Endothelial")
epi_region <- spe$sub.anno %in% epi.order
immune_region <- spe$sub.anno %in% c(lym.order, mac.order,"Granulocyte")
fib_region <- spe$sub.anno %in% c("FAP+Fibroblast","αSMA+Fibroblast",
                                  "CollagenI+Fibroblast")

# spe <- patchDetection(spe, patch_cells = stromal_region,
#                               colPairName = "expansion_interaction_graph",
#                               expand_by = 0, min_patch_size = 50,
#                               img_id = "ImageNb")
# colData(spe)$stromal_region <- !is.na(spe$patch_id)
# spe <- patchDetection(spe, patch_cells = stromal_region,
#                       colPairName = "expansion_interaction_graph",
#                       expand_by = 20, min_patch_size = 50,
#                       img_id = "ImageNb")
# colData(spe)$stromal_milieu <- (!is.na(spe$patch_id))&spe$stromal_region==F
colData(spe)$stromal_milieu <- NA
colData(spe)$stromal_region <- NA

spe <- patchDetection(spe, patch_cells = epi_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 0.1, min_patch_size = 20,
                      convex=F,
                      img_id = "ImageNb")
colData(spe)$epi_region <- !is.na(spe$patch_id)
spe <- patchDetection(spe, patch_cells = epi_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 20, min_patch_size = 20,
                      convex = F,
                      img_id = "ImageNb")
colData(spe)$epi_milieu <- (!is.na(spe$patch_id))&spe$epi_region==F

spe <- patchDetection(spe, patch_cells = immune_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 0.1, min_patch_size = 20,
                      convex=F,
                      img_id = "ImageNb")
colData(spe)$immune_region <- !is.na(spe$patch_id)
spe <- patchDetection(spe, patch_cells = immune_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 20, min_patch_size = 20,
                      convex = F,
                      img_id = "ImageNb")
colData(spe)$immune_milieu <- (!is.na(spe$patch_id))&spe$immune_region==F

spe <- patchDetection(spe, patch_cells = fib_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 0.1, min_patch_size = 20,
                      convex=F,
                      img_id = "ImageNb")
colData(spe)$fib_region <- !is.na(spe$patch_id)
spe <- patchDetection(spe, patch_cells = fib_region,
                      colPairName = "expansion_interaction_graph",
                      expand_by = 20, min_patch_size = 20,
                      convex = F,
                      img_id = "ImageNb")
colData(spe)$fib_milieu <- (!is.na(spe$patch_id))&spe$fib_region==F

plotSpatial(spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")], img_id = "ImageNb", 
            node_color_by = "meta",
            scales = "free")

plotSpatial(spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")], img_id = "ImageNb", 
            node_color_by = "immune_region",
            scales = "free")
plotSpatial(spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")], 
            img_id = "ImageNb", node_color_by = "fib_region",
            scales = "free")
plotSpatial(spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")], 
            img_id = "ImageNb", node_color_by = "fib_milieu",
            scales = "free")

plotSpatial(spe[,spe$ImageNb %in% c("T001_ROI001","T001_ROI002")], img_id = "ImageNb", 
            node_color_by = "Annotation",scales = "free")

ggplot(colData(spe.test) %>% as.data.frame(), aes(x,y,color=Annotation))+geom_point()+
  theme_classic()+scale_y_reverse()

colData(spe)$epi_imm <- colData(spe)$epi_region&colData(spe)$immune_region
colData(spe)$epi_fib <- colData(spe)$epi_region&colData(spe)$fib_region
save(spe, file="speWithRegion.RData")

##### first regions area
library(reshape2)
info <- info %>% mutate(group=group1)
info$group[info$group==""] <- info$group2[info$group==""]
group.cols <- setNames(c("#ff8e7f","#ff5d49","#7fb9b3","#4c8f7e"),
                       c("pre(0-1)","post(0-1)","pre(2-3)","post(2-3)"))
region.cols <- setNames(c("#adcbe3","#f6a6b2","#ffae58","#b7ded2","#c2b7db"),
                        c("Epithelial","Fibrosis","Immune","Epi&Immune","Epi&Fibrosis"))
ct.order <- c("Epithelial","EMT","PD-L1+Epithelial",
              "S100A9+Epithelial","IDO1+Epithelial","IDO1+S100A9+Epithelial","ProliferatingEpithelial",
              "Endothelial","αSMA+Fibroblast","FAP+Fibroblast","CollagenI+Fibroblast",
              "CD45RO+CD4T","Treg","CD38+Treg","CD4T&B","CD4&CD8T","GZMB+CD8Teff",
              "CD45RO+ActivatedCD8T","CD38+CD8T","NK&CD8T","CD8T&B","B","CXCL13+B","ProliferatingB",
              "NK","Granulocyte","LAMP3+DC","C1QC+Mφ", "AntigenPresentingMφ","SPP1+Mφ","CD163+Mφ",
              "ProliferatingMφ","SPP1+IDO+Mφ","InfiltratingMφ","PD-L1+Mφ")


stats.total <- colData(spe) %>% as.data.frame() %>% dplyr::count(ImageNb)
stats <- colData(spe) %>% as.data.frame() %>% group_by(ImageNb) %>%
  summarise_at(c("fib_region","epi_region","immune_region","epi_imm","epi_fib"), sum)
stats <- inner_join(stats.total, stats, by="ImageNb")
stats <- left_join(stats, info, by=c("ImageNb"="roi_id"))
comparisons.p=list(c("pre(0-1)","post(0-1)"),c("pre(2-3)","post(2-3)"))
comparisons.u=list(c("pre(0-1)","pre(2-3)"),c("post(0-1)","post(2-3)"))
stats[,c("fib_region","immune_region","epi_region","epi_imm","epi_fib")] <- 
  stats[,c("fib_region","immune_region","epi_region","epi_imm","epi_fib")]/stats$n
write.csv(stats[,c(1,3:7)], "ml/region_freq.csv", quote=F, col.names=T,
          row.names=F)

stats.m <- stats %>% group_by(group,patient_name) %>% 
  summarise_at(c("fib_region","epi_region","immune_region","epi_imm","epi_fib"), mean) %>% 
  mutate(group=factor(group, levels = c("pre(0-1)","post(0-1)","pre(2-3)","post(2-3)")))
stats.m <- stats.m %>% filter(!(group=="pre(2-3)" & patient_name=="ZXJ"))
stats.m <- melt(stats.m, id.var=c("group","patient_name"))
upper <- max(stats.m$value)
stats.m$variable  <- factor(stats.m$variable, levels = c("epi_region","immune_region","fib_region",
                                                         "epi_imm","epi_fib"))

p <- ggpaired(stats.m, x = "group", y = "value",
              id="patient_name",
              color = "group", 
              line.color = "gray", line.size = 0.4
)+scale_color_manual(values=group.cols)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  stat_compare_means(comparisons=comparisons.p,step.increase=c(0,0),
                     vjust = -.3, tip.length = 0.01,aes(label = ..p.signif..),
                     method = "t.test", paired=T)+
  stat_compare_means(comparisons=comparisons.u,
                     # step.increase=c(.1,.3),
                     vjust = -.3, tip.length = 0.01,
                     aes(label = ..p.signif..),
                     label.y=c(1.3*upper,1.6*upper),
                     method = "t.test")+
  facet_wrap(.~variable, ncol=5, labeller = 
               labeller(variable=setNames(c("Fibrosis","Epithelial","Immune","Epi&Immune","Epi&Fibrosis"), 
                                 c("fib_region","epi_region","immune_region","epi_imm","epi_fib"))))+
  theme_classic()+xlab("")+ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black",vjust = 0.5), 
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold",color = "black"))
ggsave("figures/region_area_compare.pdf", p, height = 3, width = 10)
# save(spe, file="speWithAnnotation.RData")

##### ct's proportion for each area
stats.epi.total <- colData(spe) %>% as.data.frame() %>% filter(epi_region==T) %>% 
  dplyr::count(ImageNb)
stats.epi <- colData(spe) %>% as.data.frame() %>% filter(epi_region==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.epi$total <- stats.epi.total[match(stats.epi$ImageNb, stats.epi.total$ImageNb),"n"]
stats.epi$freq <- stats.epi$n/stats.epi$total
stats.epi$group <- info$group[match(stats.epi$ImageNb, info$roi_id)]
stats.epi <- stats.epi %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))
write.csv(stats.epi[,c(1,2,5)], file="ml/epi_region_freq.csv", col.names = T,
          row.names = F, quote = F)

epi.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.epi, factor(stats.epi$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.epi, factor(stats.epi$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.epi, factor(stats.epi$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.epi, factor(stats.epi$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>% 
  mutate(region="Epithelial")

##### for immune region
stats.immune.total <- colData(spe) %>% as.data.frame() %>% filter(immune_region==T) %>% 
  dplyr::count(ImageNb)
stats.immune <- colData(spe) %>% as.data.frame() %>% filter(immune_region==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.immune$total <- stats.immune.total[match(stats.immune$ImageNb, stats.immune.total$ImageNb),"n"]
stats.immune$freq <- stats.immune$n/stats.immune$total
stats.immune$group <- info$group[match(stats.immune$ImageNb, info$roi_id)]
stats.immune <- stats.immune %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))
write.csv(stats.immune[,c(1,2,5)], file="ml/immune_region_freq.csv", col.names = T,
          row.names = F, quote = F)


immune.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>%
  mutate(region="Immune")

##### for immune region
stats.immune.total <- colData(spe) %>% as.data.frame() %>% filter(immune_region==T) %>% 
  dplyr::count(ImageNb)
stats.immune <- colData(spe) %>% as.data.frame() %>% filter(immune_region==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.immune$total <- stats.immune.total[match(stats.immune$ImageNb, stats.immune.total$ImageNb),"n"]
stats.immune$freq <- stats.immune$n/stats.immune$total
stats.immune$group <- info$group[match(stats.immune$ImageNb, info$roi_id)]
stats.immune <- stats.immune %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))

immune.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.immune, factor(stats.immune$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>%
  mutate(region="Immune")

##### for fibroblast region
stats.fib.total <- colData(spe) %>% as.data.frame() %>% filter(fib_region==T) %>% 
  dplyr::count(ImageNb)
stats.fib <- colData(spe) %>% as.data.frame() %>% filter(fib_region==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.fib$total <- stats.fib.total[match(stats.fib$ImageNb, stats.fib.total$ImageNb),"n"]
stats.fib$freq <- stats.fib$n/stats.fib$total
stats.fib$group <- info$group[match(stats.fib$ImageNb, info$roi_id)]
stats.fib <- stats.fib %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))
write.csv(stats.fib[,c(1,2,5)], file="ml/fib_region_freq.csv", col.names = T,
          row.names = F, quote = F)

fib.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.fib, factor(stats.fib$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.fib, factor(stats.fib$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.fib, factor(stats.fib$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.fib, factor(stats.fib$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>%
  mutate(region="Fibrosis")

##### for epi&immune region
stats.ei.total <- colData(spe) %>% as.data.frame() %>% filter(epi_imm==T) %>% 
  dplyr::count(ImageNb)
stats.ei <- colData(spe) %>% as.data.frame() %>% filter(epi_imm==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.ei$total <- stats.ei.total[match(stats.ei$ImageNb, stats.ei.total$ImageNb),"n"]
stats.ei$freq <- stats.ei$n/stats.ei$total
stats.ei$group <- info$group[match(stats.ei$ImageNb, info$roi_id)]
stats.ei <- stats.ei %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))
write.csv(stats.ei[,c(1,2,5)], file="ml/ei_region_freq.csv", col.names = T,
          row.names = F, quote = F)

ei.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.ei, factor(stats.ei$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ei, factor(stats.ei$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ei, factor(stats.ei$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ei, factor(stats.ei$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>%
  mutate(region="Epi&Immune")

##### for epi&fib region
stats.ef.total <- colData(spe) %>% as.data.frame() %>% filter(epi_fib==T) %>% 
  dplyr::count(ImageNb)
stats.ef <- colData(spe) %>% as.data.frame() %>% filter(epi_fib==T) %>% 
  mutate(sub.anno=factor(sub.anno, levels = ct.order)) %>% 
  dplyr::count(., ImageNb, sub.anno, .drop=F)
stats.ef$total <- stats.ef.total[match(stats.ef$ImageNb, stats.ef.total$ImageNb),"n"]
stats.ef$freq <- stats.ef$n/stats.ef$total
stats.ef$group <- info$group[match(stats.ef$ImageNb, info$roi_id)]

stats.ef <- colData(spe) %>% as.data.frame() %>% filter(epi_fib==T) %>%
 dplyr::count(., ImageNb, sub.anno, .drop=F) %>% group_by(ImageNb) %>%
    mutate(freq=prop.table(n))
stats.ef <- stats.ef %>% group_by(sub.anno) %>% mutate(s.freq=scale(freq))
write.csv(stats.ef[,c(1,2,5)], file="ml/ef_region_freq.csv", col.names = T,
          row.names = F, quote = F)

ef.region.stats <- Reduce(rbind,list(Reduce(rbind,lapply(split(stats.ef, factor(stats.ef$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","pre(2-3)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="pre(0-1)v.s.pre(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ef, factor(stats.ef$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("pre(0-1)","post(0-1)")) %>%
                    mutate(group=recode(group,"pre(0-1)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order, group="post(0-1)v.s.pre(0-1)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ef, factor(stats.ef$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","post(0-1)")) %>%
                    mutate(group=recode(group,"post(2-3)"=0, "post(0-1)"=1)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(0-1)v.s.post(2-3)") %>% as.data.frame(),

Reduce(rbind,lapply(split(stats.ef, factor(stats.ef$sub.anno)), function(x){
  modelOut <- glm(group ~ s.freq, dat = x %>% filter(group %in% c("post(2-3)","pre(2-3)")) %>%
                    mutate(group=recode(group,"post(2-3)"=1, "pre(2-3)"=0)), 
                  family = binomial)
  model_tidy <- as.data.table(broom::tidy(modelOut)) 
  estimate <- model_tidy[2, .(estimate, p.value, std.error)]
  return(estimate)
})) %>% mutate(ct=ct.order,group="post(2-3)v.s.pre(2-3)") %>% as.data.frame())) %>%
  mutate(region="Epi&Fibrosis")
ef.region.stats$estimate[is.na(ef.region.stats$estimate)] <- 0
ef.region.stats$p.value[is.na(ef.region.stats$p.value)] <- 1


