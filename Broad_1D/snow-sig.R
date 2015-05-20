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

## load the breakpoints raRDS if it exists
suppressMessages(suppressWarnings(library(data.table, quietly=TRUE)))
ra = data.table()
if (file.exists(opt$raRDS)) {
  print(paste("importing events from", opt$raRDS))
  ra <- readRDS(opt$raRDS)
  print(paste("imported", nrow(ra), "breakpoints"))
}

## check that there is some data provided
if (!file.exists(opt$VCFlist) && nrow(ra) == 0) {
  print(print_help(parseobj))
  stop("Need to supply either file with a list of VCFfiles or pre-loaded ra rds")
}

print(paste("...setting output directory to", opt$outdir))
setwd(opt$outdir)
#print("...sourcing all R libraries")
##source.all()

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

################################################
#### DATA LOADING AND FILTERING
################################################

## load the data from VCF files is we didn't already load data
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


## length filter
print(paste("removing events less than",opt$lengthCutoff,"bases"))
ra <- ra[is.na(ra$span) | ra$span >= opt$lengthCutoff]
print(paste("...after sanger filter:", nrow(ra)))

## add in the sifs data
print("...reading SIFS file and annotating regions")
sifs <- read.sifs()
if (!"analysis.id" %in% colnames(ra))
  ra$analysis.id <- gsub("(.*?).svcp_.*", "\\1", ra$uid)
if (!"short_code" %in% colnames(ra))
  ra$short_code <- sifs$short_code[match(ra$analysis.id, sifs$Tumour.Analyzed.Sample.Aliquot.GUUID)]

## cache it
print("...saving the imported breaks to ra.rds")
saveRDS(ra, "ra.rds")

##########################################
### BASIC PLOTTING OF THE RAW DATA
##########################################
print("...plotting counts by tumor type")
g <- sig.type.plot(ra)
pdf("counts_by_type.pdf", width=15); print(g); dev.off()

## plot top hits bar graph
g <- sig.plot.gene.hits(ra)
pdf("raw_gene_hits.pdf", width=25, height=7); print(g[[1]]); dev.off()
g <- sig.plot.gene.hits(ra, exclude.fragile=T, exclude.large=T)
pdf("raw_gene_hits_nofrag.pdf", width=25, height=7); print(g[[1]]); dev.off()

## plot top hits bar graph per tumor type
setkey(ra, short_code)
mclapply(unique(ra$short_code), function(x)
         {
           g <- sig.plot.gene.hits(ra, x)
           pdf(paste(x,"raw_gene_hits.pdf", sep=""), width=25, height=7); print(g[[1]]); dev.off()
           pdf(paste(x,"enrichment.pdf", sep=""), width=25, height=7); print(g[[3]]); dev.off()
         }, mc.cores=opt$cores)

## plot top hits bar graph per tumor type, excluding large and fragile genes
setkey(ra, short_code)
mclapply(unique(ra$short_code), function(x)
         {
           g <- sig.plot.gene.hits(ra, x, exclude.fragile=T, exclude.large=T)
           pdf(paste(x,"raw_gene_hits_nofrag.pdf", sep=""), width=25, height=7); print(g[[1]]); dev.off()
         }, mc.cores=opt$cores)


## matrix plot
#print('...plotting event 2D matrix')
#plot.matrix(opt, ra, bin=5e6)

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
sifs <- read.sifs()

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

print("...Reformattting data into GRanges for plotting")
ra.bound <- breaks.dt2grl(ra)

### plot all the hits above -log10(q) = 1
drt[, n.log10.q :=  -log(qval,10)]
which.hits <-  which(drt$n.log10.q > 1)
hitnum <- 0 # counter: used?

td.rg.cgc <- gTrack(cgc.genes)

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

