library(fgsea)
library(purrr)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(purrr)
library(magrittr)
library(ComplexHeatmap)
library(tibble)
library(stringr)

#GSEA to understand the contibution of cell types to the tissue

#first gather top 100 DGE from each cell type in vectors, the refer to your old code

Ex_N<-DEG_10x.big_fil_N_EX%>% filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)
Inh_N<-DEG_10x.big_fil_N_IN%>%  filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)
Astro<-DEG_10x.big_fil_Astro%>%  filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)
Oligo<-DEG_10x.big_fil_Oligo%>%  filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)
Micro<-DEG_10x.big_fil_Micro_PVM%>%  filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)
OPC<-DEG_10x.big_fil_OPC%>%  filter(abs(avg_log2FC)>0.26 & p_val_adj<0.05)%>% slice_min(p_val_adj, n=100)


Cell_type_GO <- list(Ex_N=Ex_N$Genes, Inh_N=Inh_N$Genes,Astro=Astro$Genes,
                       Oligo=Oligo$Genes, Micro=Micro$Genes, OPC=OPC$Genes )

GSEA <- function(DE_file) {
  #DE file is the ouput of DE from a Seurat object
  library(dplyr)
  S4table = DE_file
  gene_list = S4table$avg_log2FC
  names(gene_list) = S4table$Genes
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  head(gene_list)
  gene_list
  myGO = Cell_type_GO
  fgRes<- fgsea::fgsea(pathways =myGO, 
                       stats = gene_list,
                       minSize=5,
                       maxSize=500, nPermSimple = 10000) %>% 
    as.data.frame() 
  fgRes$Enrichment = ifelse(fgRes$NES > 0,"Up-regulated","Down-regulated")
  fgRes$FDR <- p.adjust(fgRes$pval, method = "fdr")
  fgRes <- fgRes %>% mutate(minuslog_padj= -log10(FDR))
  fgRes = fgRes[order(fgRes$minuslog_padj, decreasing = TRUE),]
  fgRes = fgRes[order(fgRes$NES, decreasing = TRUE),]
}



conditonDElist <- list(Border=DEG_10x.big_fil_Border, MAC=DEG_10x.big_fil_MAC,L1=DEG_10x.big_fil_L1,
                       L2_3=DEG_10x.big_fil_L2_3, L5_6=DEG_10x.big_fil_L5_6, VAC=DEG_10x.big_fil_VAC, WM=DEG_10x.big_fil_WM  )


GSEA_conditionDElist <- conditonDElist %>% map (GSEA)




L1 <- ggplot(GSEA_conditionDElist$L1, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="L1") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

L1


ggsave("Cell_type_RS_GSEA_L1.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_L1.eps", width = 2 , height = 2.5)



L2 <- ggplot(GSEA_conditionDElist$L2_3, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="L2/3") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

L2


ggsave("Cell_type_RS_GSEA_L2.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_L2.eps", width = 2 , height = 2.5)


L5_6 <- ggplot(GSEA_conditionDElist$L5_6, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="L5/6") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

L5_6


ggsave("Cell_type_RS_GSEA_L5_6.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_L5_6.eps", width = 2 , height = 2.5)



WM <- ggplot(GSEA_conditionDElist$WM, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="WM") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

WM


ggsave("Cell_type_RS_GSEA_WM.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_WM.eps", width = 2 , height = 2.5)


Border <- ggplot(GSEA_conditionDElist$Border, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="Border") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

Border


ggsave("Cell_type_RS_GSEA_Border.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_Border.eps", width = 2 , height = 2.5)


MAC<- ggplot(GSEA_conditionDElist$MAC, aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_segment(aes(reorder(pathway, abs(NES)), xend=pathway, y=0, yend=abs(NES))) +
  geom_point(aes(size=0.02, fill=pathway), shape=21, stroke=0.2, show.legend = FALSE) +
  coord_flip() +
  labs(x="", y="", title="MAC") +
  theme_bw() +
  guides(fill=FALSE)+ theme(plot.title = element_text(hjust = 0.5)) # This will remove the fill legend

MAC

ggsave("Cell_type_RS_GSEA_MAC.png", width = 2 , height = 2.5)
ggsave("Cell_type_RS_GSEA_MAC.eps", width = 2 , height = 2.5)



(L1|L2|L5_6)/(WM|Border|MAC)+ plot_layout(widths = c(0.5, 0.5), heights = c(0.8, 0.8))

ggsave("Cell_type_RS_GSEA_Layers.png", width = 7 , height = 4, units = "in", dpi = 450)
ggsave("Cell_type_RS_GSEA_Layers.eps", width = 7 , height =4,units = "in", dpi = 450)


