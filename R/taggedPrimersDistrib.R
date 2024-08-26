
taggedPrimersDistrib <- function(primersFile, t2s, fastqR1In, fastqR2In, allowedMis = 0, outputFolder = "demult", subsample = .05, splitHeader = " ", progressBar = T) {

  # primersFile is a name of a fasta file
  # t2s is a tag to sample array containing a single multiplexed library "run", "sample, "forward", "reverse"

  suppressMessages({
    library(ShortRead)
    library(lattice)
    library(Biostrings)
    library(ggplot2)
  })

  ## validate t2s
  if (!is.data.frame(t2s)) t2s <- read.table(t2s, header = T, sep=",", stringsAsFactors = F)
  libID <- t2s[1,1]
  ## validate R1 and R2
  if (identical(fastqR1In, fastqR2In))  {
    stop(libID, ": identical R1 and R2. stop.")
  }

  # If more than a lib in the t2s
  if (length(unique(t2s$run)) > 1)
  {
    stop(libID, ": more than one lib in the t2s, not allowed...")
  }
  lib_name <- unique(t2s$run)

  # uniqueness of sample name check
  if (length(t2s$sample) != length(unique(t2s$sample)))
  {
    stop(libID, ": There are duplicate samples in the tag to sample.. stop")
  }

  # read primer file
  primers <- readFasta(primersFile)

  ## Nreads in the fastqs. Since 4 lines per read in a fastq file:
  Nreads <- as.numeric(system2("gunzip", args = paste0(fastqR1In, " -c | wc -l"), stdout = TRUE)) / 4
  samplingReads <- floor(Nreads * subsample )
  message(libID, ": Sampling ", samplingReads, " reads (",subsample * 100,"%) from a total of ", Nreads)

  ## reading both fastq by chunk // see here: http://bioconductor.org/packages/devel/bioc/vignettes/ShortRead/inst/doc/Overview.pdf
  fqR1INstr <- FastqSampler(fastqR1In, n = samplingReads, ordered=TRUE)
  fqR2INstr <- FastqSampler(fastqR2In, n = samplingReads, ordered=TRUE)
  seed <- sample(10000,1)
  set.seed(seed);fqR1IN <- yield(fqR1INstr) ## seed important to have paired reads from both files.
  set.seed(seed);fqR2IN <- yield(fqR2INstr)

  # IDs are supposed to be matching all along !
  if (!identical(sub(" .*", "", id(fqR1IN)), sub(" .*", "", id(fqR2IN)))) {
    stop(libID, ": Not the same order in the fastq files. Maybe check header of the fastq to parse correctly?")
  }

  # make a named vector of it
  prim <- as.character(sread(primers))
  names(prim) <- as.character(id(primers))

  # split fwd and rev
  fr <- sapply(strsplit(basename(names(prim)), "-"), `[`,1)
  fr_ <- unique(fr)

  fwds <- prim[grep(fr_[1], names(prim))]
  revs <- prim[grep(fr_[2], names(prim))]

  ## prepare the array
  tagdist <- array(NA, c(length(fwds)+1, length(revs)+1))
  rownames(tagdist) <- c(names(fwds), "rev_only")
  colnames(tagdist) <- c(names(revs), "fwd_only")

  # counter
  cpt <- 1
  # total combis
  combs <- length(fwds) * length(revs)

  # progress bar
  n_iter <- combs # Number of iterations of the loop

  # Initializes the progress bar
  if (progressBar) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = n_iter, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "=")

  for (fwd_ in names(fwds))
  {
    # ge the fwd primer
    fwd <- prim[as.character(fwd_)]

    ## re oriente the reads that match the fwd primer in the R2
    # we need to remove reads that are below the width of primer
    fqR1IN_ <- fqR1IN[width(fqR1IN) > as.numeric(nchar(fwd)) & width(fqR2IN) > as.numeric(nchar(fwd))]
    # fetch the same reads in R2
    fqR2IN_ <- fqR2IN[width(fqR1IN) > as.numeric(nchar(fwd)) & width(fqR2IN) > as.numeric(nchar(fwd))]

    idy <- vmatchPattern(fwd, sread(narrow(fqR2IN_, start = 1, end =as.numeric(nchar(fwd)))), fixed=FALSE, max.mismatch=allowedMis)
    # get the indexing
    yy <- elementNROWS(idy)
    ## if fwd found in R2, we swap the read from R2 file to R1 file and vice versa
    tmpR2 <- fqR2IN_[which(yy==1)]
    tmpR1 <- fqR1IN_[which(yy==1)]
    ##
    fqR1IN_[which(yy==1)] <- tmpR2
    fqR2IN_[which(yy==1)] <- tmpR1
    ####
    # redo the indexing in the R1
    idx <- vmatchPattern(fwd, sread(narrow(fqR1IN_, start = 1, end =as.numeric(nchar(fwd)))), fixed=FALSE, max.mismatch=allowedMis)
    # get the indexing
    xx <- elementNROWS(idx)
    tmpR1 <- fqR1IN_[which(xx==1)]
    tmpR2 <- fqR2IN_[which(xx==1)]

    # get the ones without the fwd also
    tmpR1no <- fqR1IN_[which(xx==0)]
    tmpR2no <- fqR2IN_[which(xx==0)]

    # paste in tagdist
    tagdist[fwd_, "fwd_only"] <- sum(xx)
    # counter for checking reverse
    cptx <- sum(xx)

    # now for each reverse primer associated with the fwd
    for (rev_ in names(revs))
    {
      #revComb <- t2s_fwdComb$reverse[1]
      # get the rev seq and reverse complement it
      rev <- prim[as.character(rev_)]
      # we need to remove reads that are below the width of primer
      tmpR1_ <- tmpR1[width(tmpR1) > as.numeric(nchar(rev)) & width(tmpR2) > as.numeric(nchar(rev))]
      # fetch the same reads in R2
      tmpR2_ <- tmpR2[width(tmpR1) > as.numeric(nchar(rev)) & width(tmpR2) > as.numeric(nchar(rev))]
      # also for the fwd that did not match the tagged fwd
      tmpR1no_ <- tmpR1no[width(tmpR1no) > as.numeric(nchar(rev)) & width(tmpR2no) > as.numeric(nchar(rev))]
      # fetch the same reads in R2
      tmpR2no_ <- tmpR2no[width(tmpR1no) > as.numeric(nchar(rev)) & width(tmpR2no) > as.numeric(nchar(rev))]

      ## now count the rev
      # check rev primer in all R2 reads
      idy <- vmatchPattern(rev, sread(narrow(tmpR2_, start = 1, end = nchar(rev))), fixed=FALSE, max.mismatch=allowedMis)
      # check fwd primer in rev comp reads
      yy <- elementNROWS(idy)
      ## get the reads
      tmpR1_ <- tmpR1_[which(yy==1)]
      tmpR2_ <- tmpR2_[which(yy==1)]
      # paste in tagdist
      tagdist[fwd_, rev_] <- sum(yy)
      # and check within the reads without the fwd to count the reads with inly the rev
      idyno <- vmatchPattern(rev, sread(narrow(tmpR2no_, start = 1, end = nchar(rev))), fixed=FALSE, max.mismatch=allowedMis)
      # check fwd primer in rev comp reads
      yyno <- elementNROWS(idyno)
      ## get the reads
      tagdist["rev_only", rev_] <- sum(yyno)
      # follow up
      #message(paste0(cpt ,"/", combs, " - primer combination ",fwd_, " x ", rev_))
      if (progressBar) setTxtProgressBar(pb, cpt)
      cpt <- cpt + 1
    }
  }

  df <- data.frame(stack(tagdist)[,, drop = F])
  ## if blank in the t2s
  t2s$Combination <- "Expected sample"
  t2s$Combination[grep("blank", t2s$sample, ignore.case = T)] <- "Blank control"
  t2s$Combination[grep("neg", t2s$sample, ignore.case = T)] <- "PCR neg. control"

  t2s$Combination <- factor(t2s$Combination, c("Expected sample", "PCR neg. control",  "Blank control"))

  p_dist <- ggplot(df, aes(x = col, y = row, fill = log10(value+1))) +
    geom_tile(color = "grey") + ## linewidth = .3 scale_fill_gradient(low = "white", high = "lightgreen") +
    geom_point(inherit.aes = F, data = t2s, mapping = aes(x = reverse, y = forward, shape = Combination), size = 3) +
    labs(x = "reverse primer", y = "forward primer", fill = "#Reads (log10)", shape = "") +
    scale_shape_manual(values = c(`Expected sample` = 21, `Blank control` = 23, `PCR neg. control` = 22)) +
    scale_fill_distiller(palette = "GnBu", direction = 1) +  #(palette = "virids", labels = c(0, 1, 2, 3, 4, 5), values = c(0, 1, 2, 3, 4, 5)) +
    theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0("Distribution of ", samplingReads, " sampled reads in ", lib_name))

  ggsave(filename = paste0(outputFolder, "TaggedPrimersDistrib_",lib_name,".pdf"), p_dist, width = 7, height = 6)
  write.table(tagdist, paste0(outputFolder, "TaggedPrimersDistrib_",lib_name,".tsv"), col.names = NA, quote = F, sep = "\t")

  ## free memory
  close(fqR1INstr)
  close(fqR2INstr)
  invisible(gc())
  if (progressBar) close(pb)
}
