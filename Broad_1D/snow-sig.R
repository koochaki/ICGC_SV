#!/usr/bin/env Rscript

##################
#### DEAL WITH OPTIONS
##################
library(optparse)

option_list = list(
    make_option(c("-o", "--outdir"),       type = "character", default = "./",  help = "Folder to place output files in"),
    make_option(c("-c", "--cores"),        type = "numeric", default = 1,     help = "Number of cores to use"),
    make_option(c("-w", "--binwidth"),     type = "numeric", default = 1e6,   help = "Width, in bps, of bins"),
    make_option(c("-a", "--bedA"),    type = "character", default = "",   help = "BED file (A) to test A-B hypothesis with"),
    make_option(c("-b", "--bedB"),    type = "character", default = "",   help = "BED file (B) to test A-B hypothesis with"),
    make_option(c("-r", "--raRDS"),        type = "character", default = "",  help = "RDS file containing all events in a data.table (if not provided, calcualted instead)"),
    make_option(c("-p", "--hitplots"),        type = "logical", default = FALSE,  help="Make the hit plots"),
    make_option(c("-f", "--fsRDS"),        type = "character", default = "",  help = "RDS file containing F-score background"),
    make_option(c("-s", "--rRDS"),        type = "character", default = "",  help = "RDS file containing F-score background"),
    make_option(c("-n", "--numpermute"),   type = "numeric", default = 1e4,   help = "RDS file containing all events in a data.table (if not provided, calcualted instead)"),
    make_option(c("-v", "--VCFlist"),      type = "character", default = "",  help = "Number of times to permute. Default: 1e4"),
    make_option(c("-l", "--lengthCutoff"), type = "numeric", default = 1000,  help = "event Length cutoff (minimum bases). Default: 1000"))

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

# make output directory if necessary
if (!file.exists(opt$outdir))
  if (!dir.create(opt$outdir,recursive=TRUE))
  stop(paste("Cannot create ",opt$outdir,sep=" "))


#source(file.path(Sys.getenv('GIT_HOME'), 'isva', 'TaigaSig', 'sig.setup.R'))
  
suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))
ra = data.table()
if (file.exists(opt$raRDS)) {
  print(paste("importing events from", opt$raRDS))
  ra <- readRDS(opt$raRDS)
}

if (!file.exists(opt$VCFlist) && nrow(ra) == 0) {
  print(print_help(parseobj))
  stop("Need to supply either file with a list of VCFfiles or pre-loaded ra rds")
}

bedA.exists <- file.exists(opt$bedA)
if (!bedA.exists && nchar(opt$bedA)) 
  stop(paste("BED file supplied for A does not exist:", opt$bedA))
bedB.exists <- file.exists(opt$bedB)
if (!bedB.exists && nchar(opt$bedB)) 
  stop(paste("BED file supplied for B does not exist:", opt$bedB))
if (!(nchar(opt$bedA) == 0) ==(nchar(opt$bedB)==0))
  stop(print("Can supply zero or two BED files, but not one (for A-B test)"))

print(paste("...setting output directory to", opt$outdir))
setwd(opt$outdir)

print("...sourcing all R libraries")
source.all()

print(paste("Output directory:", getwd()))
print(paste("CPU cores:", opt$cores))
print(paste("Bin size:", opt$binwidth))
print(paste("Number of permutations:", opt$numpermute))
print(paste("Minimum event length:", opt$lengthCutoff,"bases"))
if (file.exists(opt$VCFlist))
  print(paste("VCF list:", opt$VCFlist))
if (nrow(ra) > 0)
  print(paste("Loaded", nrow(ra), "breakpoints from file", opt$raRDS))


chr_colors = structure(c("#0B9599","#05D21F","#33CAC9","#81535E","#000000","#23A7E9","#2B75D2","#9B9B51",
  "#D04AE0","#07B7CC","#C90B9C","#05AC12","#B22F6C","#D05958","#911FC0","#7A1C38",
  "#294C4A","#17D65B","#29FDA1","#076FE4","#B534BE","#02F600","#DBA809","#704C91",
  "#6A1E70"), names=c(1:22, 'X', 'Y', 'M'))

plot.cancer.genes <- function(ra, opt, max.rank = NULL) {

  ## data table to GRanges to do overlaps
  gr <- dt2gr(ra)
  uid = mcols(gr)$uid
  mcols(gr) = NULL
  grl <- split(gr, uid)
  
  ###############################################
  ## check if anything overlaps with cancer genes
  ###############################################
  cgc.genes <- track.load('cgc')
  fo <- gr.findoverlaps(gr, gr.pad(cgc.genes, 5e4))
  tab <- table(fo$subject.id)
  gene.nums = as.numeric(names(tab))
  names(tab) <- cgc.genes$gene[gene.nums]
  tab.lenscale = tab / width(cgc.genes[gene.nums])
  
  ## get the genes track data
  suppressWarnings(td.rg.cgc <- track.refgene(genes = cgc.genes$gene, height=3))


  ord = order(tab.lenscale, decreasing=T)
  dir.create(file.path(opt$outdir, "CGCgenes"), showWarnings=FALSE)

  if (is.null(max.rank))
    max.rank = length(ord)
  #assert("max.rank must be <= number of cancer genes (440)", max.rank <= length(ord))
  
  ###################
  ## MAKE PLOTS FOR ALL CANCER GENES
  ###################
  dum <- sapply(seq_along(ord)[1:max.rank], function(x) 
                {
                  gene = names(tab.lenscale[ord[x]])
                  ## find the partners
                  gene.window = gr.pad(cgc.genes[cgc.genes$gene == gene], 5e4)
                  rar.hits = grl[ceiling(gr.findoverlaps(gr, gene.window)$query.id/2)]
                  suppressWarnings(windows <- streduce(grbind(gr.pad(streduce(rar.hits), 5e4), gene.window)))
                  nindiv = length(unique(mcols(rar.hits)$individual))
                  
                  print(paste("plotting CGC gene", gene))
                  pdf(file.path(opt$outdir, "CGCgenes", paste("rank_",sprintf("%04d",x),"_",names(tab.lenscale[ord[x]]), "_NIndiv_", nindiv, ".pdf", sep="")))
                                        #pdf(file.path(opt$outdir, "ROS1_special.pdf"))  
                  td.rg.cgc$xaxis.newline = T
                  td.rg.cgc$xaxis.cex = 0.5
                  td.rg.cgc$xaxis.cex.label = 0.5
                  td.rg.cgc$xaxis.nticks = 2
                  display(td.rg.cgc, links=rar.hits, window=windows)
                  dev.off()
                })
  
}


plot.matrix <- function(opt, ra, bin=5e6) {
  
  grt  = gr.tile(gr.all(), w=bin)
  mat <- sig.tri.matrix(ra, grt, log=TRUE)
  mat.raw <- sig.tri.matrix(ra, grt, log=FALSE)
  allopts = list()
  td.mat  <- do.call('trackData', c(allopts, list(grt, mdata=mat, triangle=TRUE,
                                                  cmap.min=min(mat), cmap.max=max(mat)+0.1, track.name='Triangle',
                                                  height=25, sep.lwd=0.5, m.bg.col='white',
                                                  track.name='Breakpoint-pair Heatmap', islog=TRUE, xaxis.nticks=0,
                                                  xaxis.prefix="", xaxis.chronly=TRUE)))
  td.mat.raw  <- do.call('trackData', c(allopts, list(grt, mdata=mat.raw, triangle=TRUE,
                                                      cmap.min=min(mat.raw), cmap.max=max(mat.raw)+1, track.name='Triangle',
                                                      height=25, sep.lwd=0.5, m.bg.col='white',
                                                      track.name='Breakpoint-pair Heatmap', islog=FALSE, xaxis.nticks=0,
                                                      xaxis.prefix="", xaxis.chronly=TRUE)))
  
  pdf("overview_matrix.pdf", width=12, height=15)
  display(td.mat, windows=streduce(gr.all()))
  dev.off()
  
  pdf("overview_matrix_raw.pdf", width=12, height=12)
  display(td.mat.raw, windows=streduce(gr.all()))
  dev.off()
  
  ## per chrom
  grt  = gr.tile(gr.all(), w=bin)
  td.mat  <- do.call('trackData', c(allopts, list(grt, mdata=mat, triangle=TRUE,
                                                  cmap.min=min(mat)+0.1, cmap.max=max(mat)+0.1, track.name='Triangle',
                                                  height=25, sep.lwd=0.5, m.bg.col='white',
                                                  track.name='Breakpoint-pair Heatmap', islog=TRUE, xaxis.nticks=0,
                                                  xaxis.prefix="", xaxis.chronly=TRUE)))
  if (FALSE)
  dum <- lapply(seq(23), function(x) {
    print(paste("plotting for chr", x))
    pdf( paste("chr",x,"_matrix.pdf",sep=""), width=12, height=12)
    display(td.mat, windows=streduce(gr.all()[x]))
    dev.off()
  })
  
  
}

#################
## power.law.plot
#################

round.n <- function(x, n) n * round(x / n)

lm_eqn <- function(m){
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                              list(a = format(coef(m)[1], digits = 2),
                                                 b = format(coef(m)[2], digits = 2),
                                                r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

power.law.plot <- function(pspan, plotname) {
  print(paste('...making power-law plot for', plotname))
  ex <- ecdf(pspan[pspan > 0])
  #grid <- c(seq(0,10000,5), seq(10000,10^7,100))
  df = data.frame(pspan = pspan, cdf = ex(pspan), pspan.log = log(pspan,10), cdf.log = log(ex(pspan), 10))
  df <- df[df$cdf.log > -1000,] ## get rid of -Inf
  lm.s <- lm(df$cdf.log ~ df$pspan.log)
  lm.eq = lm_eqn(lm.s)

  yb = c(0.01,0.02,0.03,0.04,0.05,0.10,0.2, 0.3, 0.4,0.5,0.75,1)

  pdf(opt$outdir, plotname, width=5,height=5)
  print(g <- ggplot() + geom_point(data=df, aes(x=pspan.log, y=cdf.log), size=1, color="blue") +
    geom_smooth(data=df, aes(x=pspan.log, y=cdf.log), method="lm",se=F) +
        geom_text(aes(x=4.5,y=log(0.01,10), label=lm.eq), label=lm.eq, parse=T) +
    theme_bw() + ylab("CDF") + xlab("Span (bp)") +
    scale_y_continuous(breaks=log(yb,10), labels=format(yb,digits=2), limits=c(log(yb[1],10), 0)) +
    scale_x_continuous(breaks=seq(0,7), labels=parse(text=paste('10', seq(0,7), sep='^'))))
  dev.off()
}

##############################
#### MAIN SIGNIFICANCE PIPELINE
##############################

## load the data
if (nrow(ra) == 0){
  print('...importing VCF files')
  ra <- import.vcf(opt$VCFlist, cache.file = "", mc.cores=opt$cores)
}

## remove samples with too many
print(paste('...removing samples with too many calls and small rearrangements. Original Number of breakpoints:', nrow(ra)))
ra[, sample.count := nrow(.SD), by='individual']
ra <- ra[ra$sample.count < 2000, ]
print(paste("...after removal:", nrow(ra)))      

## get the spans
#@!!! should move span processing to import.vcf
print('...calculating event spans for all events')
ra[, span := abs(end[1]-end[2]) * (seqnames[1] == seqnames[2]), by=uid]
ra$span[ra$span == 0] <- NA

## length filter
print(paste("removing events less than",opt$lengthCutoff,"bases"))
ra <- ra[is.na(ra$span) | ra$span >= opt$lengthCutoff]
print(paste("...after sanger filter:", nrow(ra)))

## cache it
print("...saving the imported breaks to ra.rds")
saveRDS(ra, "ra.rds")

#######################
#### CODE FOR METHOD COMPARISON
#######################
## set the logical for methods
##   s  = mcols(ra.bound)$CALLER == "S"
##   d = mcols(ra.bound)$CALLER == "D"
##   ds = mcols(ra.bound)$CALLER == "DS"

## ## make the overlap pie chart
## cols = c(ovlp.hue, dran.hue, snow.hue)
## levs = c("Both", "dRanger", "SnowmanSV")
## df <- data.frame(Caller=c("SnowmanSV", "dRanger", "Both"), counts=c(sum(s, na.rm=T), sum(d, na.rm=T), sum(ds, na.rm=T)), stringsAsFactors=FALSE)
## g <- ggplot(df, aes(x=factor(1), y=counts, fill=Caller)) + geom_bar(stat='identity', position="fill") + coord_polar(theta="y") + xlab("") + ylab("") +
##   scale_fill_manual(values=cols, breaks=levs)
## pdf(file.path(opt$outdir, "caller_overlap_pie.pdf")); print(g); dev.off()
###############################

## make the 1/L distribution
#power.law.plot(span[which(d)], "power_law_D.pdf")

if (!opt$hitplots || TRUE) {

  ## matrix plot
  print('...plotting event 2D matrix')
  plot.matrix(opt, ra, bin=5e6)

################################################
## Start the significance testing ##
################################################

  ## load data for 1e6 binned mapping
  print('...loading 1d mappability scores')
  dload <- sig.load.bins('1e6')
  drt <- dload$grt
  
  ## calculate the 2d background from 1/L and 1d mappablities
  print('...calculating 2d background with 1/L, using 1d mappabilities')
  mat.background = sig.2d.background.from.1d.probs(drt)

  ## place events into bins
  fo <- sig.bin.counts(ra, drt)
                                        #foA <- sig.bin.counts(ra, gr2dt(gr.bedA))

## 
#foA.dd <- foA; setkey(foA.dd, tile.bin, event.count.per.bin); foA.dd <- foA.dd[!duplicated(foA.dd)]; foA.dd[, c("query.id", "sample", "sample.counts") := NULL]
#foA.dd[, event.density := event.count.per.bin / width(gr.bedA[tile.bin])]


## normalize background probablity matrix by total intra/inter chrom rate
print('...generating empirical background probability matrix')
empirical.background.prob <- sig.generate.normalized.matrices(ra, mat.background)

## get the total vector of background F-scores
if (!file.exists(opt$fsRDS)) {
  print('...permuting events')
  system.time(f.back <- sig.2d.permute(empirical.background.prob, ra, num.permute = opt$numpermute, method='1d', opt$cores))
  print("...saving the background F distribution to fback.rds")
  saveRDS(f.back, "fback.rds")
} else {
  print(paste('...loading F score background distribution from', opt$fsRDS))
  f.back <- readRDS(opt$fsRDS)
}

## plot the distribution of f scores
#f.back2 <- f.back[f.back > 0]
#pdf("fscore_distribution.pdf", height=3.5, width=10)
#hist(f.back2, col='black', breaks=300, xlab='F-score', ylab='Count', main='F-score distribution under NULL')
#dev.off()

## get a trackData for the bin probabilities
print('...making bin probabilities trackdata and plotting')
all.bin.probs <- colSums(empirical.background.prob)
drt$binprobs <- all.bin.probs
drt$col <- drt$border <- chr_colors[as.character(seqnames(dt2gr(drt)))]
drt$binprobs <- -log(drt$binprobs, 10)
td.binprob <- trackData(dt2gr(drt), y.field='binprobs', xaxis.newline=FALSE, xaxis.chronly=TRUE, track.name='-Log(Bin Background Probability)',
                        y0=3.3, y1=3.9, circles=TRUE, lwd.border=0.5, xaxis.nticks = 0, xaxis.prefix="")
pdf("td_bin_background_probabiltiies.pdf", width=10)
display(td.binprob, windows=gr.all())
dev.off()

## get the actual f-score
print('...calculating real f-scores')
f.real <- sig.fscore(ra, drt, empirical.background.prob)

## get the p and q values
print('...calculating pvalues and qvalues')
pq = sig.fscore2pval(f.back, f.real)
drt$pval <- pq$pval
drt$qval <- pq$qval

## make the qq-plot
print("...making qq plot")
pdf("qqplot.pdf", width=5, height=5)
qq_pval(drt$pval)
dev.off()

## make the manhattan plot
drt$n.log10.q <- -log(drt$qval, 10)
#drt$col <- drt$border <- chr_colors[as.character(seqnames(dt2gr(drt)))]
td.sig <- trackData(dt2gr(drt), y.field='n.log10.q', circles=TRUE, track.name='-log10(Q-value)', lwd.border=0.5, xaxis.chronly=T, y0=0, y1 = max(drt$n.log10.q)+0.5, xaxis.prefix="", xaxis.nticks=0)
pdf("sig.pdf", width=30, height=15)
display(td.sig, windows=gr.all())
dev.off()

## overlap with genes
print("...getting hits overlaps withs genes")
cgc.genes <- track.load('cgc')
gr.genes <- track.load('genes')
genes <- cgc.genes

fo <- gr.findoverlaps(dt2gr(drt), gr.genes)
drt$genes <- drt$cancer.genes <- ""
drt$genes[unique(fo$query.id)] <- sapply(unique(fo$query.id), function(x) {
  gn <- paste(gr.genes$gene[fo$subject.id[fo$query.id==x]], collapse=",")
  return(gn)
  })
fo <- gr.findoverlaps(dt2gr(drt), cgc.genes)
drt$cancer.genes[unique(fo$query.id)] <- sapply(unique(fo$query.id), function(x) {
  gn <- paste(cgc.genes$gene[fo$subject.id[fo$query.id==x]], collapse=",")
  return(gn)
  })


 ## need this function for below
interleave <- function(v1,v2)  {   ord1 <- 2*(1:length(v1))-1;    ord2 <- 2*(1:length(v2));    c(v1,v2)[order(c(ord1,ord2))] }

## count how many breakpoints and uniq sample bps in each bi
fo <- gr.findoverlaps(dt2gr(drt), dt2gr(ra))
fo <- gr2dt(fo)
fo[, num_bps := nrow(.SD), by='query.id']
fo[, individual := ra$individual[subject.id]]
fo[, num_uniq_samples := length(unique(individual)), by='query.id']

## read in the SIFS file
print("...reading SIFS file and annotating regions")
sifs <- data.table(read.delim("/cga/fh/pancan_data/pcawg_calls/pcawg_train2_metadata/PCAWG_Data_Freeze_Train_2.0_Paired_Download_table_2014_11_18-2171_Tumour-Normal_pairs_2014_11_18.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE))
sifs[, short_code := gsub("([A-Z]+).*", '\\1',Project.Code) ]
setkey(sifs, Tumour.Analyzed.Sample.Aliquot.GUUID)  
fo[, id := gsub("(.*)?.svcp.*", '\\1', individual)]
fo$study <- sifs[fo$id]$Study
fo$code <- sifs[fo$id]$short_code
setkey(fo, query.id)


## need this for below
num.vcf = length(unique(ra$individual))
sifs.tab <- table(sifs$short_code)
format.codes <- function(code)
  {
    tab <- table(sort(code))
    nam <- as.character(names(tab))
    tots <- sifs.tab[nam]
    paste(interleave(nam, sapply(seq_along(tots), function(x) return(paste(tab[x], tots[x], sep="/")))), collapse=" - ")
  }
fo[, code_summary := format.codes(code), by="query.id"]

top.code <- function(code)
  {
    tab <- table(sort(as.character(code)))
    nam <- as.character(names(tab))
    tots <- sifs.tab[nam]
    pval <- -log(sapply(seq_along(tots), function(x) 1-sum(dbinom(0:tab[x], size=sum(tab), prob=tots[x]/num.vcf))),10)
    paste(nam[order(pval, decreasing=T)], pval[order(pval, decreasing=T)], collapse=" - ")
 }
fo[, code_pval := top.code(code), by='query.id']

drt$code_pval = ""
drt$code_pval[fo$query.id] <- fo$code_pval
drt$code_summary=""
drt$code_summary[fo$query.id] <- fo$code_summary
drt$num_bps <- 0
drt$num_uniq_samples <- 0
drt$num_bps[fo$query.id] <- fo$num_bps
drt$num_uniq_samples[fo$query.id] <- fo$num_uniq_samples
  
drt$num_uniq_samples_NOT_rearranged <- length(unique(ra$individual)) - drt$num_uniq_samples
setkey(drt, qval, num_uniq_samples_NOT_rearranged)
write.htab(drt, file.path(getwd(), "results.html"))
print(paste("Saving results HTML to results.rds / results.html"))
saveRDS(drt, "results.rds")

print("#### SUCCESSFULLY COMPLETED ####")

}
if (opt$hitplots) {

suppressWarnings(td.refgenes <- track.refgene())
dir.create("topsighits", showWarnings=FALSE)
### plot all the hits above -log10(q) = 1
drt[, n.log10.q :=  -log(qval,10)]
which.hits <-  which(drt$n.log10.q > 1)
hitnum <- 0 # counter: used?

cgc.genes <- track.load('cgc')
td.rg.cgc <- trackData(cgc.genes)
##
print("...Reformattting data into GRanges for plotting")
radf <- data.frame(seqnames=ra$seqnames, start=ra$start, end=ra$end, uid=ra$uid, strand=ra$strand, col=ra$col, border=ra$border)
ra.ul = gr.fix(dt2gr(data.table(radf)),hg19_len)
ra.bound = split(ra.ul, ra$uid)

mcols(ra.bound)$col <- ra.ul$col[seq(1,length(ra.ul), by=2)]
mcols(ra.bound)$border <- ra.ul$border[seq(1,length(ra.ul), by=2)]

## loop through hits and plot
dum <- lapply(seq_along(which.hits), function(x) {

  y = which.hits[x]
  window = dt2gr(drt)[y]
  genes <- unique(gr.genes$gene[gr.findoverlaps(gr.genes, window)$query.id])
  genes.string <- paste(genes, collapse = " ")
  print(paste("Plotting top significance bin", x, "of", sum(drt$n.log10.q > 1), "which is chr: ", drt$seqnames[y],":", drt$start[y], "-", drt$end[y], "which contains genes", genes.string))

  rar.hits <- ra.bound[ceiling(gr.findoverlaps(grl.unlist(ra.bound), window)$query.id/2)]
  disp.window <- streduce(grbind(gr.pad(streduce(grl.unlist(rar.hits)),5e4),window))

  ## get overlaps of events with window
  fo <- gr2dt(gr.findoverlaps(disp.window, ra.ul))
  fo[, num_hits := nrow(.SD), by='query.id']
  disp.window$num_hits = 0
  disp.window$num_hits[fo$query.id] <- fo$num_hits

  this.window <- gr.genes[gr.genes$gene=='HLA-DRB5']
  ##this.window <- disp.window[which.max(disp.window$num_hits)]
  show.ra <- grl.allin(rar.hits, this.window, only=T) ## display only these rearrangments (e.g. ones where BOTH ends are in plot window)
  
  pdf(fp <- file.path("topsighits", paste("hit_", x, "q_",sprintf("%04f",drt$n.log10.q[y]),".pdf", sep="")), height=15, width=30)
  print(fp)
  td.rg.cgc$xaxis.newline = T
  td.rg.cgc$xaxis.cex = 0.5
  td.rg.cgc$xaxis.cex.label = 0.5
  td.rg.cgc$xaxis.nticks = 2
  td.rg.cgc$cex.tick = 0.5
  td.rg.cgc$track.name = "Cancer Genes"
  td.refgenes$track.name = "All Genes"
  td.rg.cgc$height = 0.3
  td.rg.cgc$xaxis.cex.label = 3
  td.rg.cgc$xaxis.cex = 3
  #td.rg.cgc$xaxis.unit=1e3
  td.rg.cgc$xaxis.suffix="Kb"
  #display(c(td.rg.cgc, td.refgenes), links=rar.hits, window=disp.window[which.max(disp.window$num_hits)]-1e5)
  display(c(td.rg.cgc, td.refgenes), links=ra.bound[grl.allin(ra.bound, this.window, only=T)], window=this.window)
  dev.off()
  hitnum = hitnum + 1
})

}


if (FALSE) {




################
###### generalized linear model
gr <- dt2gr(ra)
bg <- grl.unlist(track.load('grl.85k')) 
dt.bg <- gr2dt(bg)
  
## cancer genes
fo <- gr.findoverlaps(gr, cgc.genes)
ra$cgc <- 0
ra$cgc[fo$query.id] <- 1
fo <- gr.findoverlaps(bg, cgc.genes)
dt.bg$cgc <- 0
dt.bg$cgc[fo$query.id] <- 1

## replication timing
gr.rep <- track.load("reptime")
fo <- gr.findoverlaps(gr, gr.rep)
ra$rep <- NA
ra$rep[fo$query.id] <- gr.rep$timing[fo$subject.id]

## all genes
gr.gen <- track.load("genes")
fo <- gr.findoverlaps(gr, gr.gen)
ra$gene <- 0
ra$gene[fo$query.id] <- 1
fo <- gr.findoverlaps(bg, gr.gen)
dt.bg$gene <- 0
dt.bg$gene[fo$query.id] <- 1

## fragile sites
gr.frag <- import.ucsc("/xchip/gistic/Jeremiah/tracks/common_fragile_sites_fixed.bed")
fo <- gr.findoverlaps(gr, gr.frag)
ra$frag <- 0
ra$frag[fo$query.id] <- 1
fo <- gr.findoverlaps(bg, gr.frag)
dt.bg$frag <- 0
dt.bg$frag[fo$query.id] <- 1

## mappability
gr.100 <- readRDS('/xchip/gistic/Jeremiah/tracks/100map.rds')
fo <- gr.findoverlaps(gr, gr.100)
ra$map <- 0
ra$map[fo$query.id] <- 1
fo <- gr.findoverlaps(bg, gr.100)
dt.bg$map <- 0
dt.bg$map[fo$query.id] <- 1

fm <- data.table(response=c(rep(1, nrow(ra)), rep(0, nrow(dt.bg))), gene=c(ra$gene, dt.bg$gene), cgc=c(ra$cgc, dt.bg$cgc), frag=c(ra$frag, dt.bg$frag))

gm <- glm(response ~ gene + cgc + frag, data=fm, family=binomial)



###### fish hook
### collecting all target sets

## grab covered sites track
mappability = readRDS('/home/unix/marcin/DB/Tracks/wgEncodeCrgMapabilityAlign100mer.gr.rds')
blacklist = readRDS('/home/unix/marcin/Projects/LUAD/db/wgs.mut.blacklist.v2.rds')
covered = setdiff(gr.fix(mappability[mappability$score>0.9], hg_seqlengths(), drop = TRUE), blacklist)

## gather tracks / track sources
hg = readRDS('/home/unix/marcin/DB/ffTracks/hg19.rds')
hgcontext = readRDS('/home/unix/marcin/DB/ffTracks/hg19.context64.rds')
dnase = readRDS("/home/unix/marcin/DB/Tracks//A549_DNAse_narrowpeak.rds")
rept = readRDS('/home/unix/marcin/DB/Pubs/Koren2014/Smoothed_consensus_profile_hg19.rds')
chromhmm = readRDS('/home/unix/marcin/DB/Tracks/A549_ChromHMM.rds')
chromhmm.active = chromhmm[grepl('Act', chromhmm$type)]
chromhmm.het = chromhmm[grepl('(Het)|(Repress)', chromhmm$type)]
bigtile = readRDS('/home/unix/marcin/Projects/LUAD/mut/wgs/tracks/bigtile.density.rds')
gencode = read_gencode()
cds = gencode[gencode$type == 'CDS']; cds = split(cds, cds$gene_name)


## targets
gencode = read_gencode()
tiles500 = gr.tile(hg_seqlengths(), 5e2)
all.targets = list(
    tiles500 = tiles500,
    tfbs = '/home/unix/marcin/DB/Tracks/TFBS_UCI/HUMAN_hg19_BBLS_1_00_FDR_0_10.rds',
    chromhmm = '/home/unix/marcin/DB/Tracks/A549_ChromHMM.rds',
    cage = '/home/unix/marcin/DB/FANTOM/CAGE.gencodev10.gr.rds',
    fantom_enhancers = '/home/unix/marcin/DB/FANTOM/downloads/enhancers/human_permissive_enhancers_phase_1_and_2.bed',
    super_enhancers = '/home/unix/marcin/DB/Tracks/Hnisz2013.se.rds',
    cds_uncovered = '/home/unix/marcin/Projects/LUAD/mut/wgs/targets/cds.uncov.rds',
    lincRNA = gencode[gencode$gene_type == 'lincRNA']
    )



ra2 <- data.table(seqnames=ra$seqnames, start=ra$start, end=ra$end)
maf.wgs.gr <- dt2gr(ra2)
for (tname in names(all.targets))
    {
        WID = 1e5 ## window
        cat('Starting', tname, '\n')
         system(sprintf('mkdir -p mut/wgs/analyses/%s/', tname))
        out.path = sprintf('mut/wgs/analyses/%s/targets.rds', tname)
        if (!file.exists(out.path))
            targets = annotate.targets(all.targets[[tname]],
                covered = covered,
                events = maf.wgs.gr,
                base.context = list(track = hg, signature = list(GC = c('G', 'C'), AT = c('A', 'T')), type = 'sequence'),
                tn.context = list(track = hgcontext, signature = list(TpC = '(TC.)|(.GA)', CpG = '(.CG)|(CG.)'), grep = TRUE, type = 'sequence'),
                dnase = list(track = dnase, type = 'interval', pad = WID),
                activechrom = list(track = chromhmm.active, type = 'interval', pad = WID),
                hetchrom = list(track = chromhmm.het, type = 'interval', pad = WID),
                rept = list(track = rept, field = 'val', type = 'numeric', pad = WID),
                local.mut.density = list(track = bigtile, field = 'density', type = 'numeric', pad = 5e5),
                max.chunk = Inf, 
                mc.cores = 10, max.slice = 1e5, ff.chunk = 1e7, out.path = out.path)
    }


targ500 = readRDS('mut/wgs/analyses/tiles500/targets.rds')
clock({targ10K = aggregate.targets(targ500, rolling = 20, verbose = TRUE)})



}

