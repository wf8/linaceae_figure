
# plot_linaceae.R

library(ape)
library(phytools)
library(OutbreakTools)
library(strap)


# helper function to get coordinates from phylo plot
get_plot_coords <- function()
{
    zz = get(x="last_plot.phylo", envir=.PlotPhyloEnv)

    x = zz$xx
    y = zz$yy

    xy = as.data.frame(cbind(x,y))
    return(xy)
}


nexus = "MAP_ages_fixed.tree"
tree = read.annotated.nexus(nexus)

# ladderize the tree
#tree = ladderize(tree, right=FALSE)
#tree = ladderize(tree, right=TRUE)

posterior = vector()
internal_edges = vector()

# get annotations
for (i in 1:(length(tree$annotations)-1)) {

    if (!is.null(tree$annotations[[i]]$posterior)) {

        posterior[i] = round(tree$annotations[[i]]$posterior, 2)
        internal_edges = append(internal_edges, i)

    } else { 
        posterior[i] = 1.0
    }

}


# plot the tree with geological scale
tree$root.time = max(nodeHeights(tree))
geoscalePhylo(tree, units=c("Period", "Epoch"), boxes="Epoch", cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0.5, x.lim=c(-10,tree$root.time), lwd=1, width=1, quat.rm=TRUE)

# add posterior probs to edges
#edgelabels(posterior, edge=internal_edges, frame="n", cex=0.3, adj=c(0,-0.7))

# add bars at nodes
node_bar_col = rgb(red=0,green=0,blue=255,alpha=100,max=255)
bar_width = 4
xy = get_plot_coords()
x0 = vector()
x1 = vector()
for (i in 1:(length(tree$annotations)-1)) {
    
    if (!is.null(tree$annotations[[i]]$`height_95%_HPD`)) {

        #print("")
        #print(tree$annotations[[i]]$`height_95%_HPD`[[1]])
        #print(tree$annotations[[i]]$`height_95%_HPD`[[2]])
        #print(tree$root.time - tree$annotations[[i]]$`height_95%_HPD`[[1]])
        #print(tree$root.time - tree$annotations[[i]]$`height_95%_HPD`[[2]])

        x0 = append(x0, tree$root.time - tree$annotations[[i]]$`height_95%_HPD`[[1]])
        x1 = append(x1, tree$root.time - tree$annotations[[i]]$`height_95%_HPD`[[2]])
    } 

}
#max_x0 = max(x0)
#max_x1 = max(x1)
#x0 = max_x0 + tree$root.time - x0
#x1 = max_x1 + tree$root.time - x1
segments(x0=x0, x1=x1, y0=xy$y, y1=xy$y, col=node_bar_col, lwd=bar_width)

# make pdf
#pdf(nexus + ".pdf", height=11.0, width=8.5)


