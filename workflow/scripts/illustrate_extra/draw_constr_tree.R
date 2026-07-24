# --------------------------------------------------------------
# R script to reproduce the phylogenetic tree style in your image
# Packages: ape, ggtree, ggplot2, treeio (optional)
# --------------------------------------------------------------

# Install if needed (uncomment first time)
# install.packages(c("ape", "ggtree", "ggplot2"))
# BiocManager::install("ggtree")   # if ggtree is not installed yet

library(ape)
#library(ggtree)
#library(ggplot2)
#library(dplyr)
#library(phytools)
#library(phangorn)

plotConstTree <- function(infile, outfile) {
    tree <- read.tree(infile)
    pdf(outfile)
    plot(tree)
    dev.off()
}

plotConstTree(
    infile = "/data/pals2/old_maf_zenodo/configs/constraints/constr_nrl_cmaf_mafb_mafa_tj.newick",
    outfile = "FigS1a_MAFL_constr_nrl_cmaf_mafb_mafa_tj.pdf"
)

plotConstTree(
    infile = "/data/pals2/old_maf_zenodo/configs/constraints/constr_maff_mafg_mafk.newick",
    outfile = "FigS2a_MAFS_constr_maff_mafg_mafk.pdf"
)

plotConstTree(
    infile = "/data/pals2/old_maf_zenodo/configs/constraints/constr_allmafs.newick",
    outfile = "FigS3a_ALLMAFS_constr_allmafs.pdf"
)

plotConstTree(
    infile = "/data/pals2/old_maf_zenodo/configs/constraints/constr_otx1_otx2_crx_oc.newick",
    outfile = "FigS4a_OTX_constr_otx1_otx2_crx_oc.pdf"
)

quit()

pdf("Rplots.pdf", height=50, width=15)

# 1. Read your tree file (newick, nexus, phylip, etc.)
# Replace "tree.newick" with your actual file
tree <- read.tree("tree.newick")   # or read.nexus(), etc.

# 2. If your tree has internal node labels (bootstrap/support values)
#    ggtree will automatically detect them
#    Otherwise you can add them manually later
#
#
#    
#
add_support_labels<-function(node=NULL,support,
  cols=c("white","black"),cex=1.05){
  scale<-c(min(support,na.rm=TRUE),1)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]))
      points(x[i],y[i],pch=21,,cex=cex,
        bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

bs<-as.numeric(tree$node.label)

plotTree(ladderize(tree),fsize=0.7)
minmax<-add_support_labels(support=bs)
add.color.bar(leg=0.3*max(nodeHeights(tree)),
  cols=colorRampPalette(c("white","black"))(100),
  lims=minmax*100,title="bootstrap %",
  subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
  y=0.02*par()$usr[4])

dev.off()
quit()

# 3. Basic circular/unrooted plot (your tree is unrooted)
p <- ggtree(tree, layout = "rectangular") #, branch.length = "none") + 
     geom_tiplab(align = TRUE, linesize = 0.5, size = 3.5, offset = 0.5) 

# 4. Color the tip points according to support values
#    Your legend: >90% red, >80% orange, >70% yellow, >60% green, ≤50% white
p <- p + geom_tippoint(aes(color = as.numeric(label)), size = 3.5, na.rm = TRUE)
plot(p)

# Define the color scale exactly like in the figure
p <- p + scale_color_gradientn(
       colours = c("red", "#FF8C00", "yellow", "#32CD32", "white"),  # orange, limegreen
       values  = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
       breaks  = c(50, 60, 70, 80, 90),
       labels  = c("≤50%", ">60%", ">70%", ">80%", ">90%"),
       limits  = c(50, 100),
       na.value = "grey50",
       name    = "Bootstrap support"
     ) +
     guides(color = guide_legend(override.aes = list(size = 5)))

 plot(p)

# 5. (Optional) Highlight major clades with vertical labels on the right
#    Example for the two big clades you have ("MafA" and "MafB")
p <- p + geom_cladelabel(node = 45, label = "MafB", offset = 8, 
                         barsize = 1.2, fontsize = 5, extend = 0.7, 
                         align = TRUE, color = "black") +
         geom_cladelabel(node = 68, label = "MafA", offset = 8, 
                         barsize = 1.2, fontsize = 5, extend = 0.7, 
                         align = TRUE, color = "black")

# You need to find the correct node numbers for your clades:
# Use: viewClade(ggtree(tree) + geom_text2(aes(label=node), hjust=-.3))

# 6. Clean theme and expand limits so labels fit nicely
p <- p + theme(legend.position = c(0.1, 0.85),
               legend.background = element_rect(fill = "white", colour = "black"),
               legend.title = element_text(size = 11),
               legend.text  = element_text(size = 10)) +
         xlim(0, max(node.depth.edgelength(tree)) * 1.6)   # give space on the right

# 7. Show the plot
print(p)

# 8. Save high-resolution figure (exactly like the published one)
ggsave("phylogenetic_tree_like_figure.pdf", p, width = 12, height = 9, dpi = 600)
ggsave("phylogenetic_tree_like_figure.png",  p, width = 12, height = 9, dpi = 600)

# --------------------------------------------------------------
# Tips to get it identical to your figure:
#   • Use branch.length="none" → collapses branches (as in your tree)
#   • Adjust 'offset' in geom_tiplab and geom_cladelabel until alignment is perfect
#   • Find correct node numbers for MafA/MafB with:
#        ggtree(tree) + geom_text2(aes(label=node), color="red")
#   • If your support values are on internal nodes (not tips), use:
#        geom_nodelab(aes(label=label), color="black", size=2.5, hjust=1.3)
#        and then color nodes instead of tips.
# --------------------------------------------------------------
