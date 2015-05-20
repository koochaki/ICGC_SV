#' S4 class for \code{gBreak}
#'
#' Class \code{gBreak} defines an object to store rearrangement data from
#' multiple tumors. 
#'
#' @section Slots:
#' \describe{
#'   \item{data}{length(.Object) length list containing genomic data (e.g. GRanges, GrangesLists, RleLists, path to ucsc file or ffTrack file on disk)}
#'   \item{seqinfo}{Seqinfo object}
#'   \item{colormap}{length(.Object) length named list of named vectors whose entry i specifies the colormap for the meta data field of object entry i (specified by the name) and maps unique value of that data field to colors (specified by the named vector)}
#'   \item{edges}{list of data.frames of length length(.Object) which has columns $from, $to, and optional fields $col, $lwd, and $lty to specify splined edges joining data items in the corresponding track}
#' }
#' 
#' @name gBreak-class
#' @rdname gBreak-class
#' @exportClass gBreak
setClass('gBreak', representation(data = 'data.table', sdata = 'data.table', grl.data = 'GRangesList', gr.data = 'GRanges', sifs="data.table", cache='data.table'))

gBreak <- function(...) new('gBreak', ...)

setMethod('initialize', 'gBreak', function(.Object,
                                           data = NULL,
                                           sifs.file = NULL,
                                           ...)
          {
            ## fill the data
            if (is.null(data))
              {
                .Object@data = data.table()
                .Object@grl.data = GRangesList()
                .Object@gr.data = GRanges()
                .Object@sifs = data.table()
                .Object@cache = data.table()
              }
            else
              .Object@data = data

            
            return (.Object)

          })

setMethod('show', 'gBreak', function(object)
          {
            cat(sprintf('gBreak object with %s breakpoints\n', nrow(object@data)))
          })

#' #' Import data from list of VCF files
#'
#' Reads rearrangement a list of VCF files into a data table. If a cache file is provided
#' and exists, the table is loaded directly from the cache file. If the cache file
#' does not exist, it is created and the table is written to it.
#' @param vcf.file Text file containing paths to VCF files to import, one per line
#' @param mc.cores Number of threads to use
#' @return Data table containing rearrangement information from the vcf. segnames
setGeneric('importVCF', function(x, vcf.file="", mc.cores=1) standardGeneric('importVCF'))
setMethod('importVCF', 'gBreak',
          function(x, vcf.file="", mc.cores=1) {

            ## read in the VCF list file
            print('...reading in the VCF file list')
            v <- read.delim(vcf.file, header=FALSE, stringsAsFactors=FALSE)[,1]
            v <- v[grepl("vcf", v)]
            print(paste('...importing VCF files. Total of', length(v)))
            cols = sample(colors(), length(v), replace=TRUE)
            
            ## loop through the file list and import the VCF
            b = mclapply(v, function(x) {
              print(paste('   ', basename(x), " is ", match(x, v), " of ", length(v), sep=""))
              suppressWarnings(ab <- ra_breaks(x))
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

            ## fill in the data
            x@data <- rbindlist(b,fill=TRUE)
            x@data$uid = paste(x@data$individual, x@data$grl.ix, sep="_")

            ## fill in the spans
            x@data[, span := abs(end[1]-end[2]) * (seqnames[1] == seqnames[2]), by=uid]
            x@data$span[x@data$span == 0] <- NA

            ## add the GRanges data
            x <- .fillgr(x)
            
            print(paste("SV data.table size is", object.size(x@data)/1e6, "Megabytes"))
            return (x)
})

SIFS_DEFAULT = "/cga/fh/pancan_data/pcawg_calls/pcawg_train2_metadata/PCAWG_Data_Freeze_Train_2.0_Paired_Download_table_2014_11_18-2171_Tumour-Normal_pairs_2014_11_18.tsv"
setGeneric('addSifs', function(x, sifs=SIFS_DEFAULT) standardGeneric('addSifs'))
setMethod('addSifs', 'gBreak', function(x, sifs=SIFS_DEFAULT) {

  ## check that the file exists
  if (!file.exists(sifs))
    stop(paste("sifs file does not exist:", sifs))

  s <- data.table(read.delim(sifs, header=TRUE, sep='\t', stringsAsFactors=FALSE))
  s[, short_code := gsub("([A-Z]+).*", '\\1',Project.Code) ]
  setkey(s, Tumour.Analyzed.Sample.Aliquot.GUUID)

  x@sifs = s
  if (nrow(x@data) == 0)
    return (x)
  
  ## match sifs info to the breaks
  if (!"analysis.id" %in% colnames(x@data))
    x@data$analysis.id <- gsub("(.*?).svcp_.*", "\\1", x@data$uid)
  if (!"short_code" %in% colnames(x@data))
    x@data$short_code <- s$short_code[match(x@data$analysis.id, s$Tumour.Analyzed.Sample.Aliquot.GUUID)]
  
  
  return(x)

})

## add the GRanges and GRangesList format of data
setGeneric('.fillgr', function(x) standardGeneric('.fillgr'))
setMethod('.fillgr', 'gBreak', function(x) {

  ## convert to GRanges
  x@gr.data <- dt2gr(x@data)
  mcols(x@gr.data) <- DataFrame(grl.ix=as.numeric(x@data$grl.ix), grl.iix=as.numeric(x@data$grl.iix),
                                col=as.character(x@data$col), border=as.character(x@data$border),
                                uid=as.character(x@data$uid))

  ## split into GRangesList
  x@grl.data = split(x@gr.data, x@gr.data$uid)

  ## add coloring metadata to GRangesList
  mcols(x@grl.data)$col <- x@gr.data$col[seq(1,length(x@gr.data), by=2)]
  mcols(x@grl.data)$border <- x@gr.data$border[seq(1,length(x@gr.data), by=2)]

  return (x)
})

setGeneric('filter', function(x, min.length = 2000, max.sample.count = 2000) standardGeneric('filter'))
setMethod('filter', 'gBreak', function(x, min.length = 2000, max.sample.count = 2000)
          {

            orig <- nrow(x@data)
            orig.num.samples <- length(unique(x@data$individual))
            
            ## remove breakpoints from samples with too many hits
            x@data[, sample.count := nrow(.SD), by='individual']
            x@data <- x@data[x@data$sample.count < 2000, ]
            f1 <- orig - nrow(x@data)
            f2 <- length(unique(x@data$individual))
            
            ## remove breakpoints with spans that are too small
            orig2 <- nrow(x@data)
            x@data <- x@data[is.na(x@data$span) | x@data$span >= min.length]
            f3 <- orig2 - nrow(x@data)

            ## print a summary
            cat(sprintf("Removed %s breaks for bad samples, %s breaks for low span. Reduced to %s individuals from %s\n",
                        f1, f3, f2, orig.num.samples))
            
            return(x)
            
          })

setGeneric('save', function(x, file) standardGeneric('save'))
setMethod('save', 'gBreak', function(x, file) {
  saveRDS(x, file)
})

setGeneric('queryGene', function(x, gene, short.code=NA) standardGeneric('queryGene'))
setMethod('queryGene', 'gBreak', function(x, gene, short.code=NA) {
  
  setkey(x@data, GENE)
  r <- x@data[gene]

  print(table(r$short_code))
  r[, NRDS:=NULL]
  r[, uid:=NULL]
  r[, col:=NULL]
  r[, border:=NULL]
  r[, individual:=NULL]

  ## subset by short_code
  if (!is.na(short.code)) {
    setkey(r, short_code)
    r <- r[short.code]
  }
  
  return(r)
})

setGeneric('plot'

