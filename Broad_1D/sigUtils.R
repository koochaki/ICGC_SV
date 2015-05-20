###
###  Rearrangement Signifigance Utility Functions
###

#' #' Import data from list of VCF files
#'
#' Reads rearrangement a list of VCF files into a data table. If a cache file is provided
#' and exists, the table is loaded directly from the cache file. If the cache file
#' does not exist, it is created and the table is written to it.
#' @param vcf.file Text file containing paths to VCF files to import, one per line
#' @param mc.cores Number of threads to use
#' @return Data table containing rearrangement information from the vcf. segnames
import.vcf <- function(vcf.file="",cache.file="",mc.cores=1) {

    # if cache.file supplied and file exists read it
    if (nchar(cache.file)>0 && file.exists(cache.file)) {
        dt.all <- readRDS(cache.file)
    }
    # otherwise, import from list (file)
    else {
        v <- read.delim(vcf.file, header=FALSE, stringsAsFactors=FALSE)[,1]
        v <- v[grepl("vcf", v)]
        print(paste('...importing VCF files. Total of', length(v)))
        cols = sample(colors(), length(v), replace=TRUE)
        b = mclapply(v, function(x) {
            print(paste('   ', basename(x), " is ", match(x, v), " of ", length(v), sep=""))
            suppressWarnings(ab <-ra_breaks(x))
            if (length(ab) == 0)
                return (data.table())
            mcols(ab)$individual = basename(x)
            mcols(ab)$border = mcols(ab)$col = cols[match(x, v)]
            if (!inherits(ab, "GRangesList")) {
                print(paste('error for case', match(x,v)))
                return (data.table)
            }
            dt <- grl.unlist(ab)
            dt <- gr2dt(dt)
 
            ## remove some stuff to make lighter (SANGER)
            if ("TRDS" %in% colnames(dt))
                dt[, TRDS := NULL]
          
            return (dt)
        }, mc.cores=mc.cores)
        dt.all <- rbindlist(b,fill=TRUE)
        dt.all$uid = paste(dt.all$individual, dt.all$grl.ix, sep="_")
        #assert("Non-unique ids in import.vcf. Duplicated files?", length(table(table(dt.all$uid))) == 1)
        # if cache.file argument supplied, save to cache file
        if (nchar(cache.file)>0) {
            saveRDS(dt.all,cache.file)
        }
    }

    print(paste("SV data.table size is", object.size(dt.all)/1e6, "Mb"))
    return (dt.all)
}

#' Make  an overview plot given a gene set
#' 
sig.oplot <- function(grl, grt, genes=NULL, file=NULL, gene.pad=5e6, cov.pad=1e3, windows=NULL, plot.events=FALSE, display.pad=1e3, plot.ref=TRUE,
                      skip.matrix = FALSE, height=10, width=10) {

  if (is.null(file))
    file <- '~/public_html/plot.pdf'

  if (length(grl) == 0) {
    print('sig.oplot: Input data was empty')
    return()
  }

  print('...setting up')

  nullw = is.null(windows)
  if (is.null(windows))
    windows <- gr.all()

  ## load the ref genes
  if (!exists('trg'))
    assign('trg', track.refgene(), envir=.GlobalEnv)

  ## set some xaxis options
  allopts <- list()

  ## get the matrix data
  if (!skip.matrix)  {
    print('...getting matrix')
    logthemat = TRUE
    if (!"individual" %in% colnames(mcols(grl)))
      stop("To plot matrix, must have individual mcols filled in")
    mat <- sig.tri.matrix(grl, grt, log=logthemat)
    td.mat  <- do.call('trackData', c(allopts, list(grt, mdata=mat, triangle=TRUE, cmap.min=min(mat), cmap.max=max(mat)+0.1, track.name='Triangle', height=12, sep.lwd=0.5, m.bg.col='white',
                                                    track.name='Breakpoint-pair Heatmap', islog=islog, xaxis.nticks=0, xaxis.prefix="", xaxis.chronly=TRUE)))
    ppdf(display(td.mat, windows=gr.all()[1]))
  }

  print('...getting background probabilities')
  ## get the background probabilities
  dload <- sig.load.bins('1e6') ## save, with MAP information
  grt = dt2gr(dload$grt)
  M <- table(mcols(grl)$individual)
  all.mats <- sig.generate.normalized.matrices(grl, dload$mat, single.mat=TRUE)
  #f.back <- sig.2d.permute(all.mats, M, num.permute = 1000)
              
  pback.vec <- colSums(all.mats)
  ab <- t(replicate(length(M), pback.vec))
  ab <- sweep(ab, MARGIN=1, M, '*')
  fo = sig.binom(grl.dra, grt, ab)

  p.back <- sig.background(grl, grt)
  gr.w <- sig.binom(grl.dra, grt, p.back)

  gr.w$col <-  gr.w$border <- chr_colors[as.character(seqnames(gr.w))] ## color it, chr_colors in .Rprofile
  td.w <- do.call('trackData', c(allopts, list(gr.w, y.field='pval.fdr.log', circles=TRUE, track.name='-Log q', height=5, lwd.border=0.4)))
  
  ## get the gene heatmap
  if (nullw) {
    grt <- gr.tile(gr.all(), w=40e6)
    tt <- sig.gene.heat(grl, grt, windows=windows)
  }
  
  ## get the hotspot
  th <- sig.hotspot(grl, cov.pad=1e5)
  
  td.w@formatting$xaxis.chronly=TRUE
  td.w@formatting$xaxis.nticks=0
  td.w@formatting$xaxis.newline = TRUE

  if (nullw){
    tt@formatting$xaxis.chronly=TRUE
    tt@formatting$xaxis.nticks=0
  }

  if (!skip.matrix && nullw)
    td.all <- c(tt, td.w, th, td.mat)
  else if (!skip.matrix && !nullw)
    td.all <- c(td.w, th, td.mat)
  else if (skip.matrix && nullw)
    td.all <- c(tt, td.w, th)
  else if (skip.matrix && !nullw)
    td.all <- c(td.w, th)
  
  ## plot the plot
  print(paste('...sending to plot at', file))
  pdf(file, height=height, width=width)
  par(xpd=NA)
  par(mar=c(4,4,4,4))
  par(oma=c(0,0,0,0))
  display(td.all, windows=windows, y.gap=0.1)
  dev.off()

}

###############
#' sig.gene.heat
#'
#' Returns a trackData object for
#' a heatmap showing breakpoints that have
#' other end in specified genes
#' @param grl Input GRangesList of breakpoint pairs
#' @param grt GRanges specifying the tiling
#' @param genes character list of genes to calculate
#' @param windows GRanges specificy the windows you want to calc for (exclude events outside)
#' @param gene.pad integer for how far event can be from gene to count as overlap. Default 1e5
###############
sig.gene.heat <- function (grl, grt, genes = NULL, windows=NULL, gene.pad=1e5) {

  if (is.null(genes))
    genes <- c('TP53', 'ARID1A', 'GATA3', 'PIK3CA', 'BRCA1', 'BRCA2',
               'HRAS', 'NRAS', 'ALK', 'EML4', 'CHEK2')
  gr.genes <- track.load('genes')
  gr.genes <- gr.genes[gr.genes$gene %in% genes]
  if (length(gr.genes) == 0) {
    print('sig.plot: Check the names of the genes')
    return()
  }

  if (is.null(windows))
    windows = gr.all()
  
  gr <- grl.unlist(grl)
  
  ## limit to only windows events, for speed
  if (!is.null(windows)) {
    fo <- gr.findoverlaps(gr, windows)
    grl <- grl[gr$grl.ix[fo$query.id]]
    gr <- grl.unlist(grl)

    fo <- gr.findoverlaps(grt, windows)
    grt <- grt[fo$query.id]
  }

  gr.dt <- gr2dt(gr)
  gr.dt$gr.event <- seq(nrow(gr.dt))

  print('...getting overlaps')  
  ## get the overlap of tiles to events
  fo <- gr2dt(gr.findoverlaps(grt, gr))
  setnames(fo, c('query.id', 'subject.id'), c('tile.bin', 'gr.event'))
  fo$grl.ix <- gr$grl.ix[fo$gr.event]
  
  ## get the overlap of genes to tiles
  fot <- gr2dt(gr.findoverlaps(gr, gr.pad(gr.genes, gene.pad)))
  fot$gene <- gr.genes$gene[fot$subject.id]
  fot$grl.ix <- gr.dt$grl.ix[fot$query.id]

  ## if BOTH ends of a break lie in same gene, just keep one (don't double count)
  setkey(fot, grl.ix, gene)
  fot <- fot[!duplicated(fot)]
  if (nrow(fot) > 0) {
    setnames(fot, 'query.id', 'gr.event')
    fot[, c('subject.id', 'seqnames', 'start', 'end', 'strand', 'grl.ix') := NULL]

    ## combine
    fo2 <- merge(fo, fot, by='gr.event')
    setkey(fo2, gene, tile.bin)
    fo2[, count := nrow(.SD), by='tile.bin'] ## count the events
    fo2 <- fo2[!duplicated(fo2)]
    fo2[, c('gr.event', 'seqnames', 'start', 'end','strand', 'grl.ix') := NULL] # remove meaningless
    setkey(fo2, gene)
    
    cmap.max = max(fo2$count)
    cmap.min = 0
    cs <- colorRampPalette(c("light green", "yellow", "orange", "red"))(length(seq(cmap.min, cmap.max))+1)        
    #cs <- col.scale(seq(cmap.min, cmap.max), val.range=c(cmap.min, cmap.max), col.min=col.min, col.max=col.max)
  } 

  ## put the counts into the list
  outlist <- list()
  for (i in genes) {
    grt.this <- grt
    grt.this$count <- 0
    if (nrow(fot) > 0) {
      grt.this$count[fo2[i]$tile.bin] <- fo2[i]$count
      grt.this$col <- cs[grt.this$count + 1]
    } else {
      grt.this$col <- col.min
    }
    grt.this$border <- 'black'

    ## color the overlapping border differently
    grt.this$border[gr.findoverlaps(grt.this, gr.genes[gr.genes$gene == i])$query.id] <- 'purple'
    
    outlist[i] <- grt.this
  }

  ## make the list into a GRangesList
  grl.out <- GRangesList(outlist)

  print('...making tracks and getting matrix')

  ## make a trackData
  hei = 5
  tt <- trackData(grl.out, angle=0, lwd.border=2, col=NA, height=hei, ywid=(hei+1)/length(grl.out), gr.cex.label=0.5, xaxis.prefix="",
                  track.name='Events by Genes', xaxis.nticks=2, xaxis.newline=FALSE)
  return(tt)
  
}

###############
#' sig.hotspot
##############
sig.hotspot <- function(gr, cov.pad = 1e5) {

  if (inherits(gr, 'GRangesList'))
    gr <- grl.unlist(gr)

  ## make the hotspot coverage track
  gr.cov <- as(coverage(gr.pad(gr, cov.pad)), 'GRanges')

  gr.cov$col <-  gr.cov$border <- chr_colors[as.character(seqnames(gr.cov))] ## color it, chr_colors in .Rprofile
  
  #gr.cov$col <- 'black'
  #gr.cov$col[gr.cov$score > max(gr.cov$score)*0.75] <- 'orange'  
  #gr.cov$col[gr.cov$score > max(gr.cov$score)*0.95] <- 'red'
  #gr.cov$border <- 'black'

  allopts = list()
  td.cov <- do.call('trackData', c(allopts, list(gr.cov, y.field='score', circles=TRUE, track.name='Hotspots', height=5, lwd.border=0.4)))
  
  return(td.cov)
}

###############
#' sig.plot.events
###############
## gr <- gr.pad(gr, display.pad)
## set.seed(80)

## col.vec <- unique(mcols(grl)$individual)
## col.vec <- structure(sample(colours(), length(col.vec)), names=col.vec)
## gr$col <- col.vec[rep(mcols(grl)$individual, each=2)]
## gr$border <- NA
## mc <- mcols(grl)

## grl <- split(gr, gr$grl.ix)
## mcols(grl) <- mc

## #indr <- mcols(grl)$individual[gr$grl.ix[fo$query.id]]
## #tmp <- grl[gr$grl.ix[fo$query.id]]

## td.ev <- trackData(grl, draw.paths=TRUE, track.name='Individual Events', col='purple', labels.suppress=TRUE)

###############
#' sig.tri.matrix
#'
#'
#'
###############
sig.tri.matrix <- function(ra.dt, grt, log=FALSE, inter.only=FALSE) {

  set1 <- ra.dt$grl.iix == 1
  set2 <- ra.dt$grl.iix == 2
  #set1 <- seq(from=1, to=length(ra.dt)*2, by=2)
  #set2 <- seq(from=2, to=length(grl)*2, by=2)
  #indivs <- mcols(grl)$individual
  indivs <- ra.dt$individual[set1]
  #if (is.null(indivs))
  #  stop('sig.tri.matrix: need to have "individual" field as part of grl') 
  
  gr <- dt2gr(ra.dt)
  fo <- gr.findoverlaps(gr, grt)
  gr1 = gr[set1]
  gr2 = gr[set2]
  fo1 = data.table(as.data.frame(gr.findoverlaps(gr1, grt))[,6:7])
  fo2 = data.table(as.data.frame(gr.findoverlaps(gr2, grt))[,6:7])
  setnames(fo1, c('subject.id'), 'subject.id1')
  setnames(fo2, c('subject.id'), 'subject.id2')
  sn <- as.character(seqnames(grt))
  
  fo <- merge(fo1, fo2, by='query.id') # only keep events where both ends overlap somewhere in grt
  fo[, ind1 := (subject.id1-1) * length(grt) + subject.id2]
  fo[, ind2 := (subject.id2-1) * length(grt) + subject.id1]  
  setkey(fo, ind1, ind2)
  fo[, sample := indivs[query.id]]    
  fo[, num := length(unique(sample)), by=ind1]

  ## make the interchromosomal matrix
  fo.int <- data.table(s1=rep(seq_along(grt), length(grt)), s2=rep(seq_along(grt), each=length(grt)))
  fo.int[, chr1 := sn[s1]]
  fo.int[, chr2 := sn[s2]]  
  fo.int[, inter := chr1 != chr2]
  fo.int[, ind1 := (s1-1) * length(grt) + s2]
  fo.int[, ind2 := (s2-1) * length(grt) + s1]  
  
  ## dedupe samples in bin
  #setkey(fo, sample, subject.id1, subject.id2)
  #fo <- fo[!duplicated(fo)]
  #fo[, new.count := length(query.id), by=c('subject.id1', 'subject.id1')]
  
  mat = mat.ic <- matrix(nrow=length(grt), ncol=length(grt), 0)
  mat[fo$ind1] <- fo$num
  mat[fo$ind2] <- fo$num

  mat.ic[fo.int$ind1] <- fo.int$inter
  mat.ic[fo.int$ind2] <- fo.int$inter

  if (inter.only)
    mat <- mat * mat.ic
  
  if (log) {
    mat <- log(mat, 10)
    mat[mat < 0] <- 0
  }

  return(mat)
  
}


##############
#' sig.mapq.pair
#'
#' 
#'
##############

##############
#' sig.method.compare
#'
#'
#'
##############
xsig.method.compare <- function(g1, g2, pad = 100, ignore.strand=TRUE, name1='g1', name2='g2') {

  if (length(g1) != length(g2)) {
    warning('sig.method.compare: Lengths of comparators should be equal')
  }

  if (is.null(names(g1)))
    names <- as.character(seq_along(g1))
  else
    names <- names(g1)

  suppressWarnings(ix <- !sapply(g1, is.na) & !sapply(g2, is.na))
  comp <- num.overlap <- num.private.g1 <- num.private.g2 <- rep(NA, length(g1))
  
  comp[ix] <- lapply(seq_along(g1[ix]), function(x) {
    print(paste('overlapping', x, 'of', sum(ix)))
    ra.overlaps(g1[ix][[x]], g2[ix][[x]], ignore.strand=ignore.strand, pad=pad)
  })
  
  len.g1 <- sapply(g1, length)
  len.g2 <- sapply(g2, length)
  len.g1[!ix] <- NA
  len.g2[!ix] <- NA
  
  num.overlap[ix] <- sapply(comp[ix], function(x)
                            if (is.na(x[1]))
                                return(0)
                            else
                               return(length(unique(x[, 'ra1.ix'])))
                            )
  num.private.g1[ix] <- sapply(seq_along(comp[ix]), function(x) return(len.g1[ix][x] - num.overlap[ix][x]))
  num.private.g2[ix] <- sapply(seq_along(comp[ix]), function(x) return(len.g2[ix][x] - num.overlap[ix][x]))  
  
  dt <- data.table(overlap=num.overlap, g1=len.g1, g2=len.g2, name=names, private.g1= num.private.g1, private.g2 = num.private.g2)

  setnames(dt, c('g1', 'g2', 'private.g1', 'private.g2'), c(name1, name2, paste('private', c(name1, name2), sep='.')))
  return(list(dt=dt, ra=comp))
  
}

############
#' sig.load.dranger
#'
#'
#'
############
sig.load.dranger <- function(files) {

  ix <- file.exists(files)
  grl.dran <- as.list(rep(NA, length(files)))
  grl.dran[ix] <- lapply(files[ix], function(x) ra_breaks(read.delim(x, strings = F), keep.features = T))
  for (i in seq_along(grl.dran))
    if (is.null(grl.dran[[i]]))
      grl.dran[[i]] <- NA
  names(grl.dran) <- names(files)
  return(grl.dran)
  
}

############
#' sig.load.sno
#'
#'
############
sig.load.snow <- function(filer, verbose=TRUE, mapq=0, tsplit=0) {

   out <- lapply(filer, function(x) {
    print(x)
    if (!file.exists(x)) {
      warning(paste('sig.load.snow: file does not exist: ', x))
      return(NA)
    }

    abo <- read.table(x, header=T, sep=',', stringsAsFactors=FALSE)
    gr = GRanges()
    ab1 <- abo[,2:9]
    ab2 = ab1[c(4,5,6,1,2,3,8,7)]
    colnames(ab2) <- colnames(ab1)
    ab = rbind(ab1, ab2)

    if (grepl('germline', x))
      samplename = basename(gsub("/breakpoints/breakpoints.germline.txt", "", x))
    else 
      samplename = basename(gsub("/breakpoints/breakpoints.somatic.txt", "", x))
    abo$sample <- samplename
    
    ## fix x and y
    ab$chr1[ab$chr1 == '23'] <- 'X'
    ab$chr1[ab$chr1 == '24'] <- 'Y'
    ab$chr2[ab$chr2 == '23'] <- 'X'
    ab$chr2[ab$chr2 == '24'] <- 'Y'

    
    gr = with(ab, GRanges(chr1, IRanges(pos1, width=1), strand=strand1, mapq=mapq1))
    if (length(gr) == 0)
      return(GRangesList())
    gr$ra.index <- rep(seq(nrow(ab1)),2)
    strand(gr) <- ifelse(as.character(strand(gr))=="+", "-", "+") ## set the marcin convention

    ## fix X and Y
    seqnames(gr)
    grl <- split(gr, gr$ra.index)
    mcols(grl) <- abo
    grl <- grl[pmin(abo$mapq1, abo$mapq2) >= mapq & abo$tsplit >= tsplit]
    return(grl)
  })
                
 return(out)
  
}


###############
#' sig.binom
#'
#' @description Provides p and q values that the observed number of disrupted samples
#'   per bin is greater than expected by chance. Does FDR correction with BH.
#'
#' @param grl GRangesList with ra-pairs and "individual" meta-data field
#' @param grt GRanges with the tiles (bins) for the test. 
#' @param p.background 
#' @param indiv.field Field in grl specifying where the individual is stored
#' @param model Compute with either the "binomial" or Poisson Binomial "poibin"
###############
sig.binom <- function(grl, grt, p.background, indiv.field='individual', model='binomial', random=FALSE, force.background=FALSE) {

  require(poibin)
  ## check the input
  if (!indiv.field %in% colnames(mcols(grl)))
    stop('sig.binom: grl (events) must have a field describing which indiviual event is from')
  if (is.matrix(p.background))
    if (ncol(p.background) != length(grt))
      stop('sig.binom: Expecting nrow(p.background) == length(grt)')
  if (is.vector(p.background))
    if (length(p.background) != length(grt))
      stop('sig.binom: Expecting length(p.background) == length(grt)')

  if (inherits(grl, 'GRangesList'))
    gr <- grl.unlist(grl)
  
  ## get the sample ids for each BP
  indivs <- rep(mcols(grl)[, indiv.field], each=2)
  Ns = length(unique(indivs))

  if (is.matrix(p.background) && Ns != nrow(p.background))
    stop('sig.binom: Expecting nrow(p.background) == Number of samples, for p.background is matrix (Poisson Binomial implementation)')
  if (is.null(rownames(p.background)))
    warning('sig.binom: Row names of p.background *should* contain sample labels. If not, trusting that in same order as unique(mcols(grl)[, indiv.field])')
  else if (!identical(unique(indivs) ,rownames(p.background)))
    warning('sig.binom: Row names of p.background not in same order as indiv.field in grl...')

  ## find the event - bin overlaps
  fo <- gr2dt(gr.findoverlaps(gr, grt))
  setnames(fo, 'subject.id', 'tile.bin')
  orig.num <- nrow(fo)

  ## associate indiviuals with bins, de-duplicate
  fo[, sample := indivs[query.id]]
  setkey(fo, tile.bin, sample)
  fo <- fo[!duplicated(fo)]
  post.num <- nrow(fo)

  ## get the probabilities from the distribution
  fo[, count := nrow(.SD), by='tile.bin']
  setkey(fo, tile.bin) ## now that counted, get rid of samples
  fo <- fo[!duplicated(fo)]

  ## random it for debugging
  #randr <- round(rnorm(mean=mean(fo$count), sd=sd(fo$count), n=length(fo$count)))
  if (random) {
    randr <- rbinom(prob=p.background, size=Ns, n=length(p.background))
    print(t.test(randr, p.background*Ns))
    fo[, count := NULL]
    fo$count <- randr[fo$tile.bin]
  }

  if (force.background) {
    p.background <- rep(0, length(grt))
    p.background[fo$tile.bin] <- fo$count / length(M)
    gg <- sapply(p.background, function(x) rbinom(1, prob=x, size=length(M)))
    fo$count <- gg[fo$tile.bin]
  }
  
  if (is.matrix(p.background)) {
    if (model == 'binomial')
      fo[, pval := 1-pbinom(count, Ns, mean(p.background[, tile.bin]))]
    else if (model == 'poibin')
      fo[, pval := 1-ppoibin(pp=p.background[, tile.bin], kk=count, method='DFT-CF')]
 } else {
    fo[, pval := pbinom(count, Ns, p.background[tile.bin], lower.tail=FALSE) + dbinom(count, Ns, p.background[tile.bin])] # pbinom (LT=F) is X > x, we want X >= x
    #fo[, pval := pmin(pbinom(count, Ns, p.background[tile.bin]), pbinom(count, Ns, p.background[tile.bin], lower.tail=FALSE))]    
}

  ## correct for multiple hypotheses
  fo[, pval.fdr := p.adjust(pval, 'BH')]
  fo[, pval.fdr.log := -log(pval.fdr, 10)]
  fo$pval.fdr.log[is.na(fo$pval.fdr.log)] <- 0
  
  ## print some info
  print(paste("sig.binom:", post.num, "unique event/bin pairs of", orig.num, "total events"))
  print(paste("   ", post.num, "unique event/bin pairs of", orig.num, "total events"))
  print(paste("   ", "Max -log(q)", max(fo$pval.fdr.log, na.rm=TRUE), "with count:", max(fo$count)))

  return(dt2gr(fo))
  
}

#####################
#' sig.rainfall
#'
#' Produce a trackData object for a rainfall plot
#####################
sig.rainfall <- function(gr) {

  if (inherits(gr, 'GRangesList'))
    gr <- grl.unlist(gr)

  ## tmp destroy strandedness, due to strand sorting
  strand(gr) <- '*'
  gr <- sort(gr)

  sn <- as.character(seqnames(gr))
  first = sapply(unique(sn), function(x) which(sn==x)[1])

  diffr <- c(NA,diff(start(gr)))

  ## make all the first ones on chromosomes NA
  diffr[first] <- NA
  ix <- !is.na(diffr)
  diffr <- diffr[ix]
  if (any(diffr < 0))
    warning('sig.rainfall: not expecting neg distances')
  gr <- gr[ix]
  gr$diff = diffr
  
  td = trackData(gr, y.field='diffr', col='black', circles=TRUE, lwd.border=1)
  
}

pbinom.pvec <- function(p, k, log.p=FALSE) {

  mysum = prod(1-p); # the zero term
  len = length(p)
  numperm <- 1000
  if (k > 0)
    for (i in seq(1,k)) {
      sm <- mean(sapply(1:numperm, function(x) {
         ps = sample(p, length(p))
         prod(ps[1:i])*prod(1-ps[(i+1):len])*choose(length(p), i)
      }))
      mysum = mysum + sm
  }
  if (log.p)
    mysum <- log(mysum, 10)
  return (mysum)
}


pchernoff <- function(p, k) {

  mu <- mean(p)*length(p) ## mu is the <X>

  ## seems to only work if k is greater than exp.val?
  delta = k/mu - 1
  pd = 1 + delta

  tailr <- (exp(delta)/(pd^pd))^mu
  tailr[k <= mu] <- NA
  
  return(tailr)
  
}

##########
#' sig.pbin
#'
#' Find the probability for a given sample and bin
##########
sig.pbin <- function() {

  
  
}


########
#' sig.prep.tiles
#'
#'
########
sig.prep.tiles <- function(width=1e6, cutoff=NULL) {

 ## load the data
 dt.100 <- readRDS('/xchip/gistic/Jeremiah/tracks/100map.dt.rds') # keyed at seqnames, start
 
 ## make the tiles
 dt.tile <- gr2dt(gr.tile(gr.all(), width))
 dt.tile$bin <- seq(nrow(dt.tile))

 if (!is.null(cutoff))
   ix <- dt.100$score >= cutoff
 dt.100[, width := end - start + 1]

 print('...setting key and getting subs')
 setkey(dt.100, seqnames)
 dt.subs <- lapply(unique(dt.tile$seqnames), function(x) dt.100[x])
 names(dt.subs) <- unique(dt.tile$seqnames)

 lens <- structure(sapply(dt.subs, nrow), names=unique(dt.tile$seqnames))
 curr.seqname = dt.tile$seqnames[1]
 dt.sub <- dt.subs[[dt.tile$seqnames[1]]]
 
 scores <- mclapply(seq(nrow(dt.tile)), function(x) {
   if (curr.seqname != dt.tile$seqnames[x]) {
     curr.seqname = dt.tile$seqnames
     dt.sub <- dt.subs[[dt.tile$seqnames[x]]]
   }
     print(paste(x, 'of', nrow(dt.tile)))

     #iix <- seq(from=last, to=lens[dt.tile$seqnames[x]])
     ab <- dt.sub$start >= dt.tile$start[x] & dt.sub$end <= dt.tile$end[x]
     #last <- ab[length(ab)] + last
     if (is.null(cutoff))
       return(sum(dt.sub[ab]$total.score))
      else
        return(sum(dt.sub$width[ab]))
 }, mc.cores=1)
 
 dt.tile$total.score <- unlist(scores)
 return(dt.tile)
 
}

############
#' sig.background
#'
#' @description Calculates the background probabilities from
#'     from a set of bins and event. The bins must contain a total.score, equal
#'     to the mean detectibality in that bin.
############
sig.background <- function(grl, grt) {

  ## this OVERESTIMATES the counts
  M  <- table(mcols(grl)$individual)*2
  M <- M[order(names(M))]  
  Nb <- length(grt)

  #Sb = grt$total.score / sum(grt$total.score)
  Sb = grt$total.score

  Sg = sum(grt$total.score)#Sg.100 ## Sg.100 is stored in ~/.Rprofile

  ## THIS IS BETTER
  fo <- sig.bin.counts(grl, grt)
  fo[, sample.counts := nrow(.SD), by='sample']
  setkey(fo, sample)
  M <- structure(as.numeric(fo$sample.counts[!duplicated(fo)]), names=fo$sample[!duplicated(fo)])
  M <- M[order(names(M))]

  browser()
  p <- do.call('rbind', lapply(unique(names(M)), function(i) {
    q = M[i] / length(grt) * Sb / Sg
    #q = M[i] / Sg
    #thisp = 1 - (1 - q)^Sb # sample / bin specific probabilty of no event
  }))

  rownames(p) <- names(M)

  return (p)
  
}

 ##   ## check the input
##   if (is.list(grl)) {
##     if (inherits(grl[[1]], 'GRangesList'))
##       grl <- unlist(grl)
##       num.samples.auto <- sapply(grl, !is.na)
##   }
##   if (inherits(grl, 'GRangesList'))
##     gr <- grl.unlist(grl)
##   else if (inherits(grl,'GRanges'))
##     gr <- grl
##   else   
##     stop('sig.binom: Expecting GRangesList, list of GRangesLists or GRanges')


##   ##
##   if (is.null(grw) || TRUE) {
##     grt <- gr.tile(gr.all(), w=binwidth)
##   } else {
##     bin.weight = sum(grw$score * as.numeric(width(grw))) / 3e9 * binwidth
##     csum <- cumsum(grw$score * width(grw))
##     csum.mod <- csum %% ceiling(bin.weight)
##     csum.breaks <- c(1, which(c(0,diff(csum.mod)) < 0))
    
##     sn <- as.character(seqnames(grw))[csum.breaks[seq(1, length(csum.breaks)-1)]]
##     st <- start(grw)[csum.breaks[seq(1, length(csum.breaks)-1)]]
##     ed <- end(grw)[csum.breaks[seq(2, length(csum.breaks))]]

##     ix <- ed > st
##     grt = GRanges(sn[ix], IRanges(st[ix], ed[ix]))
##   }    

##############################
#' sig.load.bins
#'
#' Load the bins filled in with sequence info
##############################
sig.load.bins <- function(id = '1e6') {

  grt.1e6 <- '/xchip/gistic/Jeremiah/tracks/bins/100map.1e6.dt.rds'
  mat.1e6 <- '/xchip/gistic/Jeremiah/tracks/bins/100map.1e6.mat.rds'
  grt.1e5 <- '/xchip/gistic/Jeremiah/tracks/bins/100map.1e5.dt.rds'
  mat.1e5 <- '/xchip/gistic/Jeremiah/tracks/bins/100map.1e5.mat.rds'
  
  if (id == '1e6')
    return(list(grt=readRDS(grt.1e6), mat=readRDS(mat.1e6)))
  else if (id == '1e5')
    return(list(grt=readRDS(grt.1e5), mat=readRDS(mat.1e5)))    
  
}

############
#' sig.2d.background.from.1d.probs
#'
#' @description Returns a matrix of dimension length(dt) x length(dt)
#' providing the background probability for each bin, factoring in
#' detectability and 1/L distribution for each bin. Inter-chromosomal
#' events are stored as negative values, to be dealt with later.
#' @param dt A data table of ranges, with 1d probs stored in \code{total.score}
#' @param no.score Flag for setting whether to NOT incorporate total.score (do 1/L only). Default FALSE
###########
sig.2d.background.from.1d.probs <- function(dt, no.score=FALSE) {

  if (!inherits(dt, "data.table"))
    stop('expecting data.table')
  
  if (!"total.score" %in% colnames(dt))
    stop('required to set total.score field on dt')

  if (no.score)
    print("Calculating 2D matrix from 1D probabilities using 1/L ONLY")
  else
    print("Calculating 2D matrix from 1D probabilities using 1/L and mappability")    
  
  qs <- dt$total.score / sum(dt$total.score)
  sn <- dt$seqnames #as.character(seqnames())

  hg19_len <- c(249250621,243199373,198022430,191154276,180915260,
                171115067,159138663,146364022,141213431,135534747,
                135006516,133851895,115169878,107349540,
                102531392,90354753,81195210,78077248,59128983,
                63025520,48129895,51304566,155270560,59373566) 
  names(hg19_len) <- c(seq(22), 'X', 'Y')
  
  glen = sum(hg19_len)
  
  mat <- do.call('rbind', lapply(seq(length(sn)), function(i) {
    if (i %% 1000 == 0)
      print(paste(i, 'of', length(sn)))
    out <- 1/(abs((i-seq(length(sn))))+1)^1 #* (hg19_len[sn[i]]/glen)
    inter <- sn != sn[i]
    out[inter] <- -1 ## set to neg, will normalize later
    if (!no.score)
      out <- out * qs * qs[i] ## factor in the qscore (mappability score)
    return(out)
  }))

  return(mat)
}

##############
#' sig.2d.normalize.background
#'
#' @description Given a background probability matrix where inter-chromosomal
#' event probabilties are marked with a (-), will return a normalized probability
#' matrix with sum(mat) == 1, and with intra-chromosomal to inter-chromosomal mass
#' ratios consistent with the fraction intrachromosomal for that sample
#'
#' @param mat Matrix of background probabilities calculated from sig.2d.map.background
#' @param frac.intra The fraction of intra-chromosomal events for this sample
#' @value Normalized matrix accounting for fraction intra-chromomsal events
##############
sig.2d.normalize.background <- function(mat, frac.intra) {

  ix <- mat > 0 ## inter flagged as neg
  mat <- abs(mat) ## reset the flag
  t.intra = sum(mat[ ix])
  t.inter = sum(mat[ !ix])
  mat[ ix] <- mat[ ix] / t.intra * frac.intra #* num.events
  mat[!ix] <- mat[!ix] / t.inter * (1-frac.intra) #* num.events
  mat[mat==0] <- min(mat[mat > 0])

  print(paste('Total probability mass: ', sum(mat), 'Frac intra', frac.intra, 'Probability mass intra: ', sum(mat[ ix])))
  return(mat)

}

#' sig.2d.permute
#'
#' Permute events, using the sample-specific rate
#' @param mats Either total empirical background rate (1/L, mapping, empirical intra/inter) or list of them (per sample)
sig.2d.permute <- function(mats, ra, num.permute=100, method='1d', cores=1) {

  # if list, combine to 3D array
  #if (is.list(mats)) {
  #  mats <- do.call('abind', c(mats, along=3))    
  #}

  if (!is.list(mats))
    nr = nrow(mats)
  else 
    nr = nrow(mats[[1]])

  #if (dim(mats)[3] != num.events)
  #  stop('sig.2d.permute: num.events must be same dim as mats')

  ## build the sample num matrix
  #sample.mats = do.call('abind', c(lapply(num.events, function(x) matrix(nrow=dim(mats)[1], ncol=dim(mats)[2], x)), along=3))

  ## build the per-bin probability
  #pmat <- mats * sample.mats

  if (method=='1d') {

    num.events = table(ra$individual)/2
    
    print('...cutting triangles (sum rows to get 1D probabilties)')
    if (is.list(mats)) {
      mats.cut <- lapply(mats, colSums) ## take cuts of the triangle to get 1D
      cutlen <- length(mats.cut[[1]])
    } else {
      mc <- colSums(mats)
      cutlen = length(mc)
    }

    #samps <- lapply(seq_along(num.events), function(x) {
    #    rando <- sample(cutlen, prob=mc, replace=TRUE, size=num.events[y]*num.permute)
    #    rand.split <- split(rando, rep(seq(num.permute), each=num.events[y]))
    #  })

    print("...doing permutations")
    g <- mclapply(seq(num.permute), function(y) {
      if (y %% 50 == 0) print(paste("Permuting", y, "of", num.permute))
      ab <- sapply(seq_along(num.events), function(x) {
        rando <- sample(length(mc), prob=mc, replace=TRUE, size=num.events[x])
        m <- rep(0, length(mc))
        m[rando] <- -log(mc[rando]*num.events[x], 10)
        return (m)
      })
      return (rowSums(ab))
    }, mc.cores=opt$cores)

    return(unlist(g))
  }
}

#################
#' sig.generate.normalized.matrices
#'
#' Wrapper function to generate the normalized background probabilitiy using the
#' empircal inter-intra chromosomal event rate
#'
#' @param dt data.table of all events for all samples
#' @param mat non-empirical background probablity matrix (1/L and mapping)
#################
sig.generate.normalized.matrices <- function(dt, mat) {

  if (!grepl("span|SPAN", paste(colnames(dt), collapse="")))
    stop("span field required of dt")
  if ("span" %in% colnames(dt))
    tspan = dt$span
  if ("SPAN" %in% colnames(dt))
    tspace = dt$SPAN
  if (!any(is.na(tspan)))
    warning("Expecting NA values for spans for inter-chromosomal. Check this")
  if (!"individual" %in% colnames(dt))
    stop("Must have \"individual\" field")
  
  #mc = mcols(grl)
  #indivs = mc$individual
  indivs = dt$individual[dt$grl.iix == 1]
  #M = table(indivs)

  ## THIS IS BETTER
  #fo <- sig.bin.counts(grl, grt)
  #mc <- mcols(grl)[grl.unlist(grl)$grl.ix[fo$query.id], ]
  #indivs = mc$individual
  #M = table(indivs)
  
#  fo[, sample.counts := nrow(.SD), by='sample']
#  setkey(fo, sample)
#  M <- structure(as.numeric(fo$sample.counts[!duplicated(fo)]), names=fo$sample[!duplicated(fo)])
#  M <- M[order(names(M))]

  num.intra <- dt[ , sum(!is.na(span)),by=individual]$V1
  #num.intra <- sapply(names(M), function(x) sum(!is.na(tspan[indivs==x])))
  #frac.intra <- num.intra / M
  frac.intra = sum(num.intra) / nrow(dt)
  #assert("frac.intra must be <= 1", frac.intra <= 1)
  return(sig.2d.normalize.background(mat, frac.intra))

  #print(paste("Generating sample-specific background matricies for", length(frac.intra), "samples"))
  #all.mats <- mclapply(frac.intra, function(x) sig.2d.normalize.background(mat, x), mc.cores=1)
  #return(all.mats)
  
}

sig.fscore <- function(ra.dt, grt, mats) {

  ## get the counts
  fo <- sig.bin.counts(ra.dt, grt)
  if (is.list(mats) && !identical(sort(names(mats)),sort(unique(fo$sample))))
    warning('sig.fscore: expecting mats to be one per sample')
  setkey(fo, sample)

  if (is.list(mats)) {
    mats.cut <- lapply(mats, colSums)
    cutlen <- length(mats.cut[[1]])
  } else {
    mc = colSums(mats)
    cutlen = length(mc)
  }

  M = table(ra.dt$individual[ra.dt$grl.iix == 1])
  name.vec <- sort(unique(fo$sample))
  
  out <- do.call('rbind',  lapply(seq_along(name.vec), function(y) {
    mm <- rep(0, cutlen)
    if (is.list(mats)) 
      mc <- mats.cut[[y]]
    hits <- fo[name.vec[y]]$tile.bin
    mm[hits] <- -log(mc[hits] * M[name.vec[y]], 10) # get F scores for this sample
    return(mm)
  }))
  out <- colSums(out)
  names(out) <- names(mats)
  
  return(out) ## sum the F scores across the samples in each bin
  
}

sig.fscore2pval <- function(f.back, f.real) {

  ecf <- ecdf(c(f.back))
  p.val <- 1 - ecf(f.real)
  p.val[p.val == 0] <- 1 / length(f.back)
  
  ## N = length(f.real)
##   p.val <- sapply(seq_along(f.real), function(x) {
##     if (x %% 100 == 0)
##       print(paste('Calculating p-val for bin', x, 'of', N))
##     sum(f.real[x] <= f.back)
##   })
##   p.val <- p.val / length(f.back)

  out <- list(pval = p.val, qval = p.adjust(p.val, 'BH'))
  
}

#' sig.bin.counts
#'
#' Place observed events into bins. Tally total number of samples that have an
#' event in the observed bin.
#' @param dt data.table containing all of the events for all of the samples
#' @param grt data.table containing the binned genome
#' @param indiv.field string that points to field storing sample identifier (default: individual)
#' @return data.table with \code{sample} and \code{sample.counts} filled in
##############
sig.bin.counts <- function(dt, grt, indiv.field = 'individual') {

  if (!inherits(dt, 'data.table'))
    stop("sig.bin.counts: Expecting data to be passed as a data.table")
  
  ## get the sample ids for each BP
  indivs <- dt$individual
  Ns = length(unique(indivs))

  ## make the Granges for the overlaps call
  gr <- dt2gr(dt)

  ## make the tile GRanges
  if (inherits(grt, "data.table"))
    gr.tiles = dt2gr(grt)
  else
    gr.tiles = grt
  
  #indivs <- rep(mcols(grl)[, indiv.field], each=2)
  #Ns = length(unique(indivs))

  ## find the event - bin overlaps
  fo <- gr2dt(gr.findoverlaps(gr, gr.tiles))
  setnames(fo, 'subject.id', 'tile.bin')
  orig.num <- nrow(fo)

  ## associate indiviuals with bins, de-duplicate
  fo[, sample := indivs[query.id]]
  setkey(fo, tile.bin, sample)
  fo <- fo[!duplicated(fo)]

  fo[, sample.counts := nrow(.SD), by='sample'] ## get total number of binned-events per sample
  fo[, event.count.per.bin := nrow(.SD), tile.bin] ## get number of events in a bin
  setkey(fo, sample) ## sort by sample count
  
  return(fo)
  
  post.num <- nrow(fo)

  ## get the probabilities from the distribution
  fo[, count := nrow(.SD), by='tile.bin']
  setkey(fo, tile.bin) ## now that counted, get rid of samples
  fo <- fo[!duplicated(fo)]
  return(fo$count)

}



################
# import.snowman.reads
################
import.snowman.reads <- function(file) {

  if (!file.exists(file))
    stop(paste('File does not exist:', file))

  dat <- read.table(file, header=TRUE, sep=',', stringsAsFactors=FALSE)
  gr.reads  <- with(dat, GRanges(contig, IRanges(start, width=nchar(seq)), strand="+", sequence=seq, rname=rname, sw_score=sw_score, cigar=paste(nchar(seq), "M", sep='')))
  return(gr.reads)
  
}


################
# import.snowman.contigs
################
import.snowman.contigs <- function(bam, bami) {

  if (!file.exists(bam))
    stop(paste('BAM does not exist:', bam))
  if (!file.exists(bami))
    stop(paste('BAM index does not exist:', bami))

  suppressWarnings(gr.contigs <- read.bam(bam=bam, bami=bami, isProperPair=FALSE, as.grl=FALSE, pairs.grl=FALSE))
  return(gr.contigs)
  
}


##################
## plot.snowman.contig
##################
plot.snowman.contig <- function(gr.this.contig, gr.this.reads) {
  
  ## make the chains
  #to.plot = 'contig_8:75093502-75101502_12'
  #gr.this.contig <- gr.contigs[gr.contigs$qname %in% to.plot]
  
  ## make the contig colored sequence grl
  cseq.set <- DNAStringSet(gr.this.contig$seq)
  names(cseq.set) <- gr.this.contig$qname
  rseq.set <- DNAStringSet(gr.this.reads$sequence)
  names(rseq.set) <- gr.this.reads$rname

  cgc <- cgChain(gr.this.contig, sn=gr.this.contig$qname)
  pac <- cgChain(gr.this.reads, sn=gr.this.reads$rname)
  #pac = paChain(cseq.set[which.max(width(cseq.set))], rseq.set)
  r2g <- gMultiply(cgc, pac)
  
  grl.contig.seq <- seq2grl(cseq.set[which.max(width(cseq.set))]) ## max because of hard clipping issue
  grs.contig <- lift(cgc, grl.contig.seq)
  td.c2g <- do.call('trackData', c(list(data=grs.contig, track.name='Contig', labels.suppress=TRUE)))
  
  ## make the reads colored sequence grl
  grs.reads <- seq2grl(rseq.set, sn=names(rseq.set))
  reads2genome <- lift(r2g, grs.reads)
  td.r2g <- do.call('trackData', c(list(data=reads2genome, track.name='Reads2Genome', labels.suppress=TRUE, height=4)))
  
  ##
  return(td.sno<-list(r2g=td.r2g, c2g=td.c2g))
  #ppdf(display(td, windows=gr.reduce(gr.this.contig)), filename=filename)
}


##
#
##
cn.compare.plot <- function(gr.cn=NULL, gr.ra1, gr.ra2) {

  ## load the copy number breaks if 
  if (is.null(gr.cn)) {
    RAVALID = c(1)
    HCC1143.CN <- read.table('~/Projects/HCC1143/HCC1143_CNBkps.csv', header=TRUE, quote="", stringsAsFactors=FALSE, sep=',')
    cnR <- with(HCC1143.CN[HCC1143.CN$RAValid %in% RAVALID,], GRanges(chr, IRanges(left, width=1), strand=ifelse(Delta < 0, '+', '-')))
    mcols(cnR) <- HCC1143.CN[HCC1143.CN$RAValid %in% c(RAVALID),]
    gr.cn <- cnR[start(cnR) != 1]
  }

  ## unlist if need
  if (inherits(gr.ra1, "GRangesList"))
    gr.ra1 <- grl.unlist(gr.ra1)
  if (inherits(gr.ra2, "GRangesList"))
    gr.ra2 <- grl.unlist(gr.ra2)
  
  ## get the distances
  dis1 <- gr.dist(gr.cn, gr.ra1)
  dis1[is.na(dis1)] <- 1e8
  dis.to1 <- apply(dis1, 1, min) + 1
  dis.to1 <- log(apply(dis1, 1, min) + 1, 10)
  
  ## get the distance (snow)
  dis2 <- gr.dist(gr.cn, gr.ra2)
  dis2[is.na(dis2)] <- 1e8
  dis.to2 <- apply(dis2, 1, min) + 1
  dis.to2 <- log(apply(dis2, 1, min) + 1, 10)
  
  ## make the plot
  df = data.frame(dist=c(dis.to1, dis.to2), Method=c(rep('1', length(dis.to1)), rep('2', length(dis.to2))))

##   p <- ggplot(data=df, aes(x=dist, color=Method)) + geom_histogram(aes(y=..density.., fill=Method), position='identity', alpha=0.3) + xlab('Log10 Distance to nearest breakpoint') +
##     geom_step(aes(y=..y..),stat="ecdf") +
##       scale_x_continuous(breaks=seq(0,8)) +
##         theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
##           scale_color_manual(values=c('blue', 'red'), labels=c('dRanger', 'Local-Assembly')) +
##             scale_fill_manual(values=c('blue', 'red'), labels=c('dRanger', 'Local-Assembly')) +
##               ggtitle('HCC1143: SegSeq breaks vs Rearrangement')

  ## display the values
  print(sum((dis.to1 < log(20e3, 10)) / length(dis.to1)))
  print(sum((dis.to2 < log(20e3, 10)) / length(dis.to2)))

  return(df)
}

## emulate GGplot colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
      hcl(h=hues, l=65, c=100)[1:n]
  }


#################
snowman.qcimport <- function(file) {

  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  ## read the table
  con  <- file(file, open = "r")
  
  df.mapq <- df.nm <- df.isize <- df.as <- df.xp <- df.len <- df.phred <- df.clip <- data.frame()
  rg <- list()
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    
    ## found a new one
    if (grepl("READGROUP", line)) {
      thisrg = gsub("READGROUP:BI:(.*)", "\\1", line)
      rg[[thisrg]] <- data.frame(readgroup = thisrg)
    } else if (grepl("total", line)) {
      rg[[thisrg]]$total = as.numeric(gsub("total,([0-9]+)", "\\1", line))
    } else if (grepl("unmap", line)) {
      rg[[thisrg]]$unmap = as.numeric(gsub("unmap,([0-9]+)", "\\1", line))
    } else if (grepl("qcfail", line)) {
      rg[[thisrg]]$qcfail = as.numeric(gsub("qcfail,([0-9]+)", "\\1", line))
    } else if (grepl("duplicate", line)) {
      rg[[thisrg]]$duplicate = as.numeric(gsub("duplicate,([0-9]+)", "\\1", line))
    } else if (grepl("supplementary", line)) {
      rg[[thisrg]]$supp = as.numeric(gsub("supplementary,([0-9]+)", "\\1", line))    
    } else if (grepl("mapq", line)) {
      df.mapq <- rbind(df.mapq, as.numeric(strsplit(gsub("mapq,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("nm", line)) {
      df.nm <- rbind(df.nm, as.numeric(strsplit(gsub("nm,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("isize", line)) {
      df.isize <- rbind(df.isize, as.numeric(strsplit(gsub("isize,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("as", line)) {
      df.as <- rbind(df.as, as.numeric(strsplit(gsub("as,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("xp", line)) {
      df.xp <- rbind(df.xp, as.numeric(strsplit(gsub("xp,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("clip", line)) {
      df.clip <- rbind(df.clip, as.numeric(strsplit(gsub("clip,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("len", line)) {
      df.len <- rbind(df.len, as.numeric(strsplit(gsub("len,([0-9]+)", "\\1", line), ",")[[1]]))
    } else if (grepl("phred", line)) {
      df.phred <- rbind(df.phred, as.numeric(strsplit(gsub("phred,([0-9]+)", "\\1", line), ",")[[1]]))
    } else {
      stop(paste("Failed to read file at line:", line))
    }
    
  }
  
  close(con)
  
  colnames(df.mapq) <- seq(from=0,to=60)
  colnames(df.nm) <- seq(from=0,to=ncol(df.nm)-1)
  colnames(df.isize) <- seq(from=0,to=ncol(df.isize)-1)
  colnames(df.xp) <- seq(from=0, to=ncol(df.xp)-1)
  colnames(df.as) <- seq(from=0, to=ncol(df.as)-1)
  colnames(df.len) <- seq(from=0, to=ncol(df.len)-1)
  colnames(df.phred) <- seq(from=0, to=ncol(df.phred)-1)
  colnames(df.clip) <- seq(from=0, to=ncol(df.clip)-1)
  readg <- sapply(rg, function(x) x$readgroup)
  
  df.mapq$readgroup  <- readg
  df.nm$readgroup    <- readg
  df.isize$readgroup <- readg
  df.xp$readgroup    <- readg
  df.as$readgroup    <- readg
  df.phred$readgroup <- readg
  df.len$readgroup   <- readg
  df.clip$readgroup  <- readg

  return(list(mapq=df.mapq, nm=df.nm, isize=df.isize, xp=df.xp, as=df.as, phred=df.phred, len=df.len, clip=df.clip))
  
}

plot.cancer.genes <- function(ra, opt, max.rank = NULL) {

  ## data table to GRanges to do overlaps
  gr <- dt2gr(ra)
  uid = mcols(gr)$uid
  mcols(gr) = NULL
  grl <- split(gr, uid)
  
  ###############################################
  ## check if anything overlaps with cancer genes
  ###############################################
  ##cgc.genes <- track.load('cgc')
  cgc = read.delim(file.path(Sys.getenv("GIT_HOME"), 'ICGC_SV', 'data', 'cancer_gene_census.tsv'), strings = F)
  genes <- readRDS(file.path(Sys.getenv("GIT_HOME"), "ICGC_SV", "data", "gr.allgenes.rds"))
  cgc.genes <- genes[genes$gene %in% cgc$Symbol]
  fo <- gr.findoverlaps(gr, cgc.genes + 5e4)
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

#' Read a sifs file
read.sifs <- function(x="/cga/fh/pancan_data/pcawg_calls/pcawg_train2_metadata/PCAWG_Data_Freeze_Train_2.0_Paired_Download_table_2014_11_18-2171_Tumour-Normal_pairs_2014_11_18.tsv")
  {
    sifs <- data.table(read.delim(x, header=TRUE, sep='\t', stringsAsFactors=FALSE))
    sifs[, short_code := gsub("([A-Z]+).*", '\\1',Project.Code) ]
    setkey(sifs, Tumour.Analyzed.Sample.Aliquot.GUUID)  
    return(sifs)
  }

#' Plot the distribution of events from breaks and sifs file
sig.type.plot <- function(ra, sifs=read.sifs("/cga/fh/pancan_data/pcawg_calls/pcawg_train2_metadata/PCAWG_Data_Freeze_Train_2.0_Paired_Download_table_2014_11_18-2171_Tumour-Normal_pairs_2014_11_18.tsv"))
  {
    if (!"analysis.id" %in% colnames(ra))
      ra$analysis.id <- gsub("(.*?).svcp_.*", "\\1", ra$uid)
    if (!"short_code" %in% colnames(ra))
      ra$short_code <- sifs$short_code[match(ra$analysis.id, sifs$Tumour.Analyzed.Sample.Aliquot.GUUID)]

    ## get counts
    ra[, per.sample.counts := nrow(.SD), by="analysis.id"]

    ## get the median-events
    ra[, median.counts := median(per.sample.counts), by=short_code]
    
    df <- data.frame(short.code = ra$short_code, analsis.id = ra$analysis.id, per.sample.counts = ra$per.sample.counts, median.counts=ra$median.counts)
    df$short.code = as.character(df$short.code)
    df <- df[!duplicated(df), ]

    df <- df[!is.na(df$short.code), ]
    df <- df[order(df$median.counts, decreasing=TRUE),]
    df$short.code <- factor(df$short.code, levels=unique(df$short.code))
    
    g <- ggplot(df, aes(short.code, per.sample.counts/2)) + geom_boxplot() + theme_bw() + xlab("") + ylab("Events per sample")
    return(g)
  }

sig.plot.gene.hits <- function(ra, short_code=NA, pad=0, exclude.fragile = FALSE, exclude.large=FALSE)
  {
    ## gene counts per type
    setkey(ra, GENE, short_code)
    ra[, gene.counts.per.type := sum(!duplicated(individual)), by=c("GENE","short_code")]
    ra[, type.count := sum(!duplicated(individual)), by="short_code"]
    UNIQ_SAMPLES = sum(!duplicated(ra$individual))
    ra[, gene.counts.per.type.p.value := min(-log(1-sum(dbinom(0:gene.counts.per.type[1], size=type.count[1], prob=gene.count.unique.sample[1] / UNIQ_SAMPLES)),10),15), by=c("GENE", "short_code")]

    ## trime down
    setkey(ra, short_code) 
    if (!is.na(short_code))
      ra <- ra[short_code]
    
    genes <- readRDS(file.path(Sys.getenv("GIT_HOME"), "ICGC_SV", "data", "gr.allgenes.rds"))
    genes <- genes[!duplicated(genes$gene), ]
    suppressWarnings(genes <- c(genes, GRanges(4, IRanges(91048684, 92523370), gene="CCSER1")))
    suppressWarnings(genes <- c(genes, GRanges(11, IRanges(78363876, 79151992), gene="TENM4")))
    suppressWarnings(genes <- c(genes, GRanges(2, IRanges(51259739, 52635055), gene="AC007682.1")))
    suppressWarnings(genes <- c(genes, GRanges(17, IRanges(31340105, 32501983), gene="ASIC2")))
    suppressWarnings(genes <- c(genes, GRanges(5, IRanges(166711804, 167691162), gene="TENM2")))
    gene.lengths <- structure(width(genes), names=genes$gene)

    ## get fragile genes
    frag <- fread(file.path(Sys.getenv("GIT_HOME"), "ICGC_SV", "data", "common_fragile_sites_fixed.bed"))
    frag <- gr.nochr(GRanges(frag$V1, IRanges(as.numeric(frag$V2), as.numeric(frag$V3))))
    frag <- sort(frag)
    genes$class = ""
    genes$class[gr.findoverlaps(genes, frag)$query.id] = "Fragile"
    genes$class[grepl("^OR", genes$gene)] <- "Olfactory"
    genes$class[grepl("^HIST", genes$gene)] <- "Histone"
    genes$class[grepl("^SNO", genes$gene)] <- "SnoRNA"
    genes$class[grepl("^HLA", genes$gene)] <- "HLA"

    ## cancer genes
    cgc = read.delim(file.path(Sys.getenv("GIT_HOME"), 'ICGC_SV', 'data', 'cancer_gene_census.tsv'), strings = F)
    genes$class[genes$gene %in% cgc$Symbol & nchar(genes$class)==0] <- "Cancer Gene Census"
    
    ## manually annotate some
    genes$class[genes$gene %in% c("CSMD1","CCSER1", "GMDS")] <- "Fragile"
    genes$class[width(genes) > 6.5e5 & genes$class == ""] <- "Large (top 1%)"
    genes$class[width(genes) > 2.45e5 & genes$class == ""] <- "Large (top 5%)"
    
    dt.genes <- data.table(gene.width=width(genes), GENE=genes$gene, gene.class=genes$class)
    setkey(dt.genes, GENE)
    setkey(ra, GENE)
    ra[, gene.count := nrow(.SD), by="GENE"]
    ra[, gene.count.unique.sample := sum(!duplicated(individual)), by="GENE"]
    setkey(ra, GENE)
    ra <- dt.genes[ra]
    
    ra[, gene.count.per.kb := nrow(.SD) / gene.width * 1000, by="GENE"]
    ra[, gene.count.unique.sample.per.kb := sum(!duplicated(individual)) / gene.width * 1000, by="GENE"]
    ra$gene.class[nchar(ra$gene.class) == 0] <- "none"
    
    ## plot params
    groups <- c("lightblue", "black", "darkgreen", "peachpuff", "grey70", "gold", "seagreen1", "seagreen4", "red")
    cols <- c("Fragile", "none", "Olfactory", "Histone", "SnoRNA", "HLA", "Large (top 1%)", "Large (top 5%)", "Cancer Gene Census")

    ## trim if need
    if (exclude.fragile)
      ra <- ra[!ra$gene.class %in% c('Fragile')]
    if (exclude.large)
      ra <- ra[!ra$gene.class %in% c('Large (top 1%)', 'Large (top 5%)')]

    
    ## plot by counts
    setkey(ra, gene.count)
    df <- data.frame(gene=ra$GENE, count=ra$gene.count.unique.sample, gene.class=ra$gene.class)
    df <- df[!duplicated(df),]
    df$gene <- as.character(df$gene)
    df$gene[is.na(df$gene)] <- df$gene[is.na(df$gene)] <- "NO GENE"
    df <- df[order(df$count, decreasing=TRUE), ]
    g1 <- ggplot(df[2:100,], aes(x=factor(gene, levels=gene), y=count, fill=gene.class), color=NA) + geom_bar(stat='identity') + theme_bw() +
      theme(axis.text.x=element_text(angle=90)) + xlab("Gene") + ylab("# Individuals") +
      scale_fill_manual(name="Gene class", values=groups, limits=cols)
    
    ## plot by counts / Kb
    setkey(ra, gene.count)
    df <- data.frame(gene=ra$GENE, count=ra$gene.count.unique.sample.per.kb, gene.class=ra$gene.class)
    df <- df[!duplicated(df),]
    df$gene <- as.character(df$gene)
    df$gene[is.na(df$gene)] <- df$gene[is.na(df$gene)] <- "NO GENE"
    df <- df[order(df$count, decreasing=TRUE), ]
    g2 <- ggplot(df[2:50,], aes(x=factor(gene, levels=df$gene), y=count, fill=gene.class), color=NA) + geom_bar(stat='identity') + theme_bw() +
      theme(axis.text.x=element_text(angle=90)) + xlab("Gene") + ylab("Breaks per Kb") +
      scale_fill_manual(name="Gene class", values=groups,limits=cols)
    
    ## plot by type pval
    df <- data.frame(gene=ra$GENE, count=ra$gene.counts.per.type.p.value, gene.class=ra$gene.class)
    df <- df[!duplicated(df),]
    df$gene <- as.character(df$gene)
    df$gene[is.na(df$gene)] <- df$gene[is.na(df$gene)] <- "NO GENE"
    df <- df[order(df$count, decreasing=TRUE), ]
    g3 <- ggplot(df[2:100,], aes(x=factor(gene, levels=gene), y=count, fill=gene.class), color=NA) + geom_bar(stat='identity') + theme_bw() +
      theme(axis.text.x=element_text(angle=90)) + xlab("Gene") + ylab("-log10 Type Enrichment P-value") +
      scale_fill_manual(name="Gene class", values=groups, limits=cols)
    
    return(list(g1, g2, g3))
  }


breaks.dt2grl <- function(ra)
  {
    ##

    hg19_len <- hg_seqlengths()
    radf <- data.frame(seqnames=ra$seqnames, start=ra$start, end=ra$end, uid=ra$uid, strand=ra$strand, col=ra$col, border=ra$border)
    ra.ul = gr.fix(dt2gr(data.table(radf)),hg19_len)
    ra.bound = split(ra.ul, ra$uid)
    
    mcols(ra.bound)$col <- ra.ul$col[seq(1,length(ra.ul), by=2)]
    mcols(ra.bound)$border <- ra.ul$border[seq(1,length(ra.ul), by=2)]

    return (ra.bound)
  }
