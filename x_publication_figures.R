# Make figures for publication

# LIBS
library(ggplot2)
library(gridExtra)

# Figure 3, tx treemap
palms <- readRDS(file=file.path('figures', 'palms_tx_treemap.RData'))
primates <- readRDS(file=file.path('figures', 'primates_tx_treemap.RData'))
tiff(file.path('figures', "figure_3.tiff"), width=14, height=7, units="cm",
     res=2750, compression='lzw')
grid.arrange(palms + theme(legend.position='none') + labs(x='a.'),
             primates + theme(legend.position='none') + labs(x='b.'),
             ncol=2)
dev.off()

# Figure 4, cl treemap
palms <- readRDS(file=file.path('figures', 'palms_cl_treemap.RData'))
primates <- readRDS(file=file.path('figures', 'primates_cl_treemap.RData'))
tiff(file.path('figures', "figure_4.tiff"), width=14, height=7, units="cm",
     res=2000, compression='lzw')
grid.arrange(palms + theme(legend.position='none') + labs(x='a.'),
             primates + theme(legend.position='none') + labs(x='b.'),
             ncol=2)
dev.off()

# Figure 5, tx pa
palms <- readRDS(file=file.path('figures', 'palms_cl_pamap.RData'))
primates <- readRDS(file=file.path('figures', 'primates_cl_pamap.RData'))
tiff(file.path('figures', "figure_5.tiff"), width=14, height=7, units="cm",
     res=2500, compression='lzw', pointsize=1)
grid.arrange(palms + theme(legend.position='none',
                           text=element_text(size=8)) + labs(x='a.'),
             primates + theme(legend.position='none',
                              text=element_text(size=8),
                              axis.text.y=element_text(size=3)) +
               labs(x='b.'),
             ncol=2)
dev.off()

# Figure 6, trees
palms <- readRDS(file=file.path('figures', 'palms_tree.RData'))
primates <- readRDS(file=file.path('figures', 'primates_tree.RData'))
tiff(file.path('figures', "figure_6.tiff"), width=14, height=14, units="cm",
     res=2000, compression='lzw')
split.screen(figs=c(1,2))
split.screen(figs=c(2,1), screen=1)
screen(3)
par(mar=c(.5,0,0,0))
plot(palms, edge.width=.6, label.offset=2, srt=-10, cex=.5)
nodelabels(text=palms$node.label, frame='none', cex=.3,
           adj=c(-.1, .5))
mtext(text='a.', side=1, cex=.75)
screen(2)
par(mar=c(.5,0,0,0))
plot(primates, edge.width=.5, label.offset=1, srt=-10, cex=.3)
nodelabels(text=primates$node.label, frame='none', cex=.2,
           adj=c(-.1, .5))
mtext(text='b.', side=1, cex=.75, line = -1)
close.screen(all.screens = TRUE)
dev.off()
