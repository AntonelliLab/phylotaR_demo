# Make figures for publication

# LIBS
library(ape)
library(ggplot2)
library(treemapify)
library(grid)
library(gridExtra)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

whitetext <- function(p) {
  p$layers <- p$layers[1]
  p + geom_treemap_text(fontface="italic", colour="white",
                        place="centre", grow=TRUE)
}
  
# Figure S2, tx treemap
palms <- whitetext(readRDS(file=file.path('figures', 'palms_tx_treemap.RData')))
primates <- whitetext(readRDS(file=file.path('figures', 'primates_tx_treemap.RData')))
caption <- paste0('Figure S2. Relative distribution of ',
                  'number of sequences (box size) and clusters (colour)',
                  ' across tribes/genera for palms (a) and primates (b).')
mylegend <- g_legend(palms + theme(legend.position = 'top') +
                       guides(fill=guide_legend(title="No. clusters")))
dev.off()
pdf(file.path('figures', "figure_S2.pdf"), width=14, height=7)
grid.arrange(palms + theme(legend.position='none') + labs(x='a.'),
             primates + theme(legend.position='none') +
               labs(x='b.'), ncol=2,
             bottom=textGrob(caption, just="centre",
                             gp=gpar(fontsize=10, fontface='bold')),
             top=mylegend)
dev.off()

# Figure S3, cl treemap
palms <- whitetext(readRDS(file=file.path('figures', 'palms_cl_treemap.RData')))
primates <- whitetext(readRDS(file=file.path('figures', 'primates_cl_treemap.RData')))
caption <- paste0('Figure S3. Relative distribution of ',
                  'number of sequences (box size) and taxa (colour)',
                  ' across clusters for palms (a) and primates (b).')
mylegend <- g_legend(palms + theme(legend.position = 'top') +
                       guides(fill=guide_legend(title="No. taxa")))
dev.off()
pdf(file.path('figures', "figure_S3.pdf"), width=14, height=7)
grid.arrange(palms + theme(legend.position='none') + labs(x='a.'),
             primates + theme(legend.position='none') + labs(x='b.'),
             ncol=2,
             bottom=textGrob(caption, just="centre",
                             gp=gpar(fontsize=10, fontface='bold')),
             top=mylegend)
dev.off()

# Figure 3, tx pa
palms <- readRDS(file=file.path('figures', 'palms_cl_pamap.RData'))
primates <- readRDS(file=file.path('figures', 'primates_cl_pamap.RData'))
tiff(file.path('figures', "figure_3.tiff"), width=7, height=14, units="cm",
     res=2500, compression='lzw', pointsize=1)
grid.arrange(palms + theme(legend.position='none',
                           text=element_text(size=8),
                           axis.text.y=element_text(size=4)) +
               labs(x='a.'),
             primates + theme(legend.position='none',
                              text=element_text(size=8),
                              axis.text.y=element_text(size=3)) +
               labs(x='b.'), heights=c(1,2))
dev.off()

# Figure 4, trees
palms <- readRDS(file=file.path('figures', 'palms_tree.RData'))
primates <- readRDS(file=file.path('figures', 'primates_tree.RData'))
tiff(file.path('figures', "figure_4.tiff"), width=14, height=14, units="cm",
     res=2000, compression='lzw')
split.screen(figs=c(1,2))
split.screen(figs=c(2,1), screen=1)
screen(3)
par(mar=c(.5,0,0,0))
plot(ladderize(palms), edge.width=.6, label.offset=2, srt=0, cex=.5,
     font = 1)
nodelabels(text=palms$node.label, frame='none', cex=.35,
           adj=c(-.1, .5))
mtext(text='a.', side=1, cex=.75)
screen(2)
par(mar=c(.5,0,0,0))
plot(ladderize(primates), edge.width=.5, label.offset=1, srt=0,
     cex=.3, font = 3)
nodelabels(text=primates$node.label, frame='none', cex=.25,
           adj=c(-.1, .5))
mtext(text='b.', side=1, cex=.75, line = -1)
close.screen(all.screens = TRUE)
dev.off()
