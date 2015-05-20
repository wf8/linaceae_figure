
# plot_linaceae.R

library(ape)
library(phytools)
library(OutbreakTools)
library(strap)
library(plotrix)

# read in tree
nexus = "MAP_ages_fixed_abr.tree"
tree = read.annotated.nexus(nexus)

# turn on pdf writer
pdf(paste0(nexus, ".pdf"), height=10.0, width=8.5)

# ladderize the tree
tree = ladderize(tree, right=FALSE)

# get root height
tree$root.time = max(nodeHeights(tree))

#tree = drop.tip(tree, "Hugonia_jenkinsii")

# get posterior annotations
posteriors = vector()
for (i in 1:length(tree$annotations)) {

    if (!is.null(tree$annotations[[i]]$posterior)) 
        #posteriors[i] = tree$annotations[[i]]$posterior
        posteriors = append(posteriors, tree$annotations[[i]]$posterior)
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
geoscalePhylo(tree, units=c("Period", "Epoch"), boxes="Epoch", cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=0.5, x.lim=c(-40,tree$root.time)+7, width=1, quat.rm=TRUE)

#############################
# add posterior probs to edges
#edgelabelstPP$edges(round(posteriors, 2), frame="n", cex=0.3, adj=c(0,-0.6))
lastPP = get("last_plot.phylo", envir = .PlotPhyloEnv)
internal_nodes = (lastPP$Ntip + 1):length(lastPP$xx)
XX = lastPP$xx[internal_nodes]
YY = lastPP$yy[internal_nodes]
        
for (i in 1:length(posteriors)) {

    if (posteriors[i] > 0.80) {
        if (i == length(posteriors))
            text(x=XX[i], y=YY[i], round(posteriors[i], 2), cex=0.3, offset=0)
#           text(x=XX[i], y=YY[i], i, cex=0.3, adj=c(0,-0.6))
                
        else {
#           x = (lastPP$xx[lastPP$edge[, 1][i]] + lastPP$xx[lastPP$edge[, 2][i]])/2
#           y = lastPP$yy[lastPP$edge[, 2][i]]
#           text(x=x, y=y, round(posteriors[i], 2), cex=0.3, adj=c(0,-0.6))

            text(x=(XX[i+1] - 0.9), y=(YY[i+1] + 0.8), round(posteriors[i], 2), cex=0.3, offset=0)
        }
    }
}


##############################
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


############################
# family labels
text(1, 3.5, "Ixonanthaceae", offset=0, pos=4, cex=0.8)
text(1, 24.5, "Linaceae", offset=0, pos=4, cex=0.8)

# subfamily labels
text(27, 30, "Linoideae", offset=0, pos=4, cex=0.7)
text(23.8, 16, "Hugonioideae", offset=0, pos=4, cex=0.7)

# fossil labels
draw.circle(0,12.73,radius=1,nv=100,border=NA,col=rgb(red=0,green=100,blue=25,alpha=100,max=255),lty=3,lwd=1)
draw.circle(40.4,47.8,radius=1,nv=100,border=NA,col=rgb(red=0,green=100,blue=25,alpha=100,max=255),lty=3,lwd=1)

############################
# add clade bars

bar_text_size = 0.6
line_x = 92.2

segments(line_x, 1, line_x, 5, lwd=2)
text(line_x + 1, 2, substitute(italic("Ochthocosmus")), offset=0, pos=4, cex=bar_text_size)
text(line_x + 1, 3.5, substitute(italic("Cyrillopsis / Ixonanthes")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 5.7, line_x, 7.2, lwd=2)
text(line_x + 1, 6.2, substitute(italic("Hebepetalum")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 7.9, line_x, 10, lwd=2)
text(line_x + 1, 8.7, substitute(italic("Roucheria")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 10.6, line_x, 11.2, lwd=2)
text(line_x + 1, 10.6, substitute(paste(italic("Hugonia"), " sect. ", italic("Durandea"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 11.8, line_x, 13.2, lwd=2)
text(line_x + 1, 12.5, substitute(italic("Philbornea / Indorouchera")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 13.9, line_x, 16.3, lwd=2)
text(line_x + 1, 14.8, substitute(paste(italic("Hugonia"), " sect. ", italic("Durandea"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 16.9, line_x, 17.5, lwd=2)
text(line_x + 1, 16.9, substitute(paste(italic("Hugonia"), " sect. ", italic("Hugoniopsis"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 18.2, line_x, 25.2, lwd=2)
text(line_x + 1, 21, substitute(paste(italic("Hugonia"), " sect. ", italic("Hugonia"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 25.9, line_x, 31.1, lwd=2)
text(line_x + 1, 29, substitute(italic("Anisadenia / Tirpitzia")), offset=0, pos=4, cex=bar_text_size)
text(line_x + 1, 27.5, substitute(italic("Reinwardtia")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 31.8, line_x, 35.2, lwd=2)
text(line_x + 1, 33, substitute(paste(italic("Linum"), " sect. ", italic("Dasylinum"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 35.9, line_x, 53.8, lwd=2)
text(line_x + 1, 45, substitute(paste(italic("Linum"), " sect. ", italic("Linum"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 54.5, line_x, 54.9, lwd=2)
text(line_x + 1, 54.5, substitute(italic("Radiola")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 55.6, line_x, 56, lwd=2)
text(line_x + 1, 55.7, substitute(paste(italic("Linum"), " sect. ", italic("Cathartolinum"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 56.7, line_x, 64.3, lwd=2)
text(line_x + 1, 60, substitute(paste(italic("Linum"), " sect. ", italic("Linopsis"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 65, line_x, 69.5, lwd=2)
text(line_x + 1, 67, substitute(paste(italic("Linum"), " sect. ", italic("Syllinum"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 70.2, line_x, 71, lwd=2)
text(line_x + 1, 70.2, substitute(italic("Sclerolinon")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 71.7, line_x, 84.3, lwd=2)
text(line_x + 1, 77, substitute(italic("Hesperolinon")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 85, line_x, 95, lwd=2)
text(line_x + 1, 90, substitute(paste(italic("Linum"), " sect. ", italic("Linopsis"))), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 95.7, line_x, 96.5, lwd=2)
text(line_x + 1, 95.9, substitute(italic("Cliococca")), offset=0, pos=4, cex=bar_text_size)

segments(line_x, 97.2, line_x, 100, lwd=2)
text(line_x + 1, 98, substitute(paste(italic("Linum"), " sect. ", italic("Linopsis"))), offset=0, pos=4, cex=bar_text_size)

# turn off pdf writer
dev.off()

