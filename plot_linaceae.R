
# plot_linaceae.R

library(ape)
library(phytools)
library(OutbreakTools)
library(strap)

# read in tree
nexus = "MAP_ages_fixed.tree"
tree = read.annotated.nexus(nexus)

# turn on pdf writer
pdf(paste0(nexus, ".pdf"), height=10.0, width=8.5)

# ladderize the tree
#tree = ladderize(tree, right=FALSE)

# get root height
tree$root.time = max(nodeHeights(tree))

#tree = drop.tip(tree, "Hugonia_jenkinsii")

# get posterior annotations
posteriors = vector()
for (i in 1:length(tree$annotations)) {

    if (!is.null(tree$annotations[[i]]$posterior)) 
        posteriors[i] = tree$annotations[[i]]$posterior
}

hpd_min = vector()
hpd_max = vector()

# add values into vectors
for(i in 1:length(tree$annotations)){

    if (!is.null(tree$annotations[[i]]$`height_95%_HPD`)) {

        min_i = tree$annotations[[i]]$`height_95%_HPD`[[1]]
        max_i = tree$annotations[[i]]$`height_95%_HPD`[[2]]

        if (min_i < tree$root.time) 
            min_i = tree$root.time - min_i
        else
            min_i = -1 * (min_i - tree$root.time)

        if (max_i < tree$root.time) 
            max_i = tree$root.time - max_i
        else
            max_i = -1 * (max_i - tree$root.time)
        

        hpd_min = append(hpd_min, min_i)
        hpd_max = append(hpd_max, max_i)

    } 
}


# plot the tree with geological scale
geoscalePhylo(tree, units=c("Period", "Epoch"), boxes="Epoch", cex.tip=0.5, cex.age=0.5, cex.ts=0.5, label.offset=0.5, x.lim=c(-40,tree$root.time)+7, width=1, quat.rm=TRUE)

# add posterior probs to edges
edgelabels(round(posteriors, 2), frame="n", cex=0.3, adj=c(0,-0.6))

# add node bars
node_bar_col = rgb(red=0,green=0,blue=255,alpha=100,max=255)
bar_width = 4
# get the coordinates of the last tree plot
lastPP = get("last_plot.phylo", envir = .PlotPhyloEnv)
internal_nodes = (lastPP$Ntip + 1):length(lastPP$xx)
XX = lastPP$xx[internal_nodes]
YY = lastPP$yy[internal_nodes]

for(i in 1:length(hpd_min)) {

    if (i == length(hpd_min))
        segments(x0=hpd_min[i], x1=hpd_max[i], y0=YY[1], y1=YY[1], col=node_bar_col, lwd=bar_width)
    else
        segments(x0=hpd_min[i], x1=hpd_max[i], y0=YY[i+1], y1=YY[i+1], col=node_bar_col, lwd=bar_width)
}


# add clade bars

bar_text_size = 0.6
line_x = 97.5

segments(line_x, 96, line_x, 100, lwd=2)
text(line_x + 1, 97.5, "Ixonanthaceae", offset=0, pos=4, cex=bar_text_size)

segments(line_x, 76, line_x, 95, lwd=2)
text(line_x + 1, 86, "Hugonioideae", offset=0, pos=4, cex=bar_text_size)

segments(line_x, 70, line_x, 75, lwd=2)
text(line_x + 1, 72, "S. Asian genera", offset=0, pos=4, cex=bar_text_size)

segments(line_x, 66, line_x, 69, lwd=2)
text(line_x + 1, 67.5, "Linum sect. Dasylinum", offset=0, pos=4, cex=bar_text_size)

segments(line_x, 46, line_x, 65, lwd=2)
text(line_x + 1, 60, "Linum sect. Linum", offset=0, pos=4, cex=bar_text_size)


# turn off pdf writer
dev.off()

