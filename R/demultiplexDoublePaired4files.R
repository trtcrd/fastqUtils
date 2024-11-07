
demultiplexDoublePaired4files <- function(primersFile, t2s, fastqR1In, fastqR2In, allowedMis = 0, 
                                          outputFolder = "demult", overwrite = F, splitHeader = " ", withIndels = F,
                                          chunkSize = 1000000, progressBar = T) {

  ## to do, trim hoverhang, i.e. trimming the reverse complement rev primer in fwd reads ?
  ## option NOT to trim the primers
  ## compared to DTD in SLIM, with the first function we did not filter NNNNNNN reads, that passed and were probably filter in dada2
  ## now we do remove these
  ## Taking those NNN reads away make the fastqs identical to DTD (just not same read order)

  # primersFile is a name of a fasta file
  # t2s is a tag to sample array containing a single multiplexed library "run", "sample, "forward", "reverse"

  
  ## vmatchPattern with indels allowed: the base function does not work: "vmatchPattern() does not support indels yet" ?? 
  # found this here: https://support.bioconductor.org/p/58350/
  vmatchPattern2 <- function(pattern, subject,
                             max.mismatch=0, min.mismatch=0,
                             with.indels=FALSE, fixed=TRUE,
                             algorithm="auto")
  {
    if (!is(subject, "XStringSet")) subject <- Biostrings:::XStringSet(NULL, subject)
    algo <- Biostrings:::normargAlgorithm(algorithm)
    if (Biostrings:::isCharacterAlgo(algo)) stop("'subject' must be a single (non-empty) string ", "for this algorithm")
    pattern <- Biostrings:::normargPattern(pattern, subject)
    max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
    min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
    with.indels <- Biostrings:::normargWithIndels(with.indels)
    fixed <- Biostrings:::normargFixed(fixed, subject)
    algo <- Biostrings:::selectAlgo(algo, pattern, max.mismatch, min.mismatch, with.indels, fixed)
    C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject, max.mismatch, min.mismatch, 
                    with.indels, fixed, algo, "MATCHES_AS_RANGES", PACKAGE="Biostrings")
    unlisted_ans <- IRanges(start=unlist(C_ans[[1L]], use.names=FALSE), width=unlist(C_ans[[2L]], use.names=FALSE))
    relist(unlisted_ans, C_ans[[1L]])
  }
  
  
  suppressMessages({
    library(ShortRead)
    library(Biostrings)
  })
  ## time tracking
  start.time <- Sys.time()

  ## validate t2s (if file name or passed data.frame)
  if (!is.data.frame(t2s)) t2s <- read.table(t2s, header = T, sep=",", stringsAsFactors = F)
  libID <- t2s[1,1]

  ## validate R1 and R2
  if (identical(fastqR1In, fastqR2In))  {
    stop(libID, ": identical R1 and R2. stop.")
  }

  ## check whether there is an output folder
  path.exp <- file.path(getwd(), outputFolder)
  if(!dir.exists(file.path(path.exp))) {
    dir.create(path.exp)
  } else if (overwrite) {
    system2("rm", args = paste("-r", path.exp))
    message(paste0(outputFolder, " folder will be overwritten."))
    dir.create(path.exp)
  } else {
    stop(paste0(outputFolder, " folder for the run already exists, stop."))
  }

  ## read primers
  primers <- readFasta(primersFile)


  # uniqueness of sample names check
  if (length(t2s$sample) != length(unique(t2s$sample))) {
    stop(libID, ": there are duplicate samples in the tag to sample, please check and fix ", names(unique(t2s$sample)[unique(t2s$sample)>1]))
  }
  # uniqueness of primer combinations
  if (length(table(paste0(t2s$forward, "_", t2s$reverse))) < nrow(t2s)) {
    stop(libID, ": some primer combinations are duplicated, please check and fix ", names(table(paste0(t2s$forward, "_", t2s$reverse))[table(paste0(t2s$forward, "_", t2s$reverse))>1]))
  }

  ## reading both fastq by chunk // see here: http://bioconductor.org/packages/devel/bioc/vignettes/ShortRead/inst/doc/Overview.pdf
  fqR1INstr <- FastqStreamer(fastqR1In, n = chunkSize)
  fqR2INstr <- FastqStreamer(fastqR2In, n = chunkSize)

  ##### progress and stats
  # how many chunks will we have?
  ## Nreads in the fastqs. Since 4 lines per read in a fastq file:
  Nreads <- as.numeric(system2("gunzip", args = paste0(fastqR1In, " -c | wc -l"), stdout = TRUE)) / 4
  Nbchunks <- ceiling(Nreads / chunkSize)
  message(libID, ": ", Nreads, " reads to demultiplex. Will do ", Nbchunks, " chunk(s) of ", chunkSize," reads.")

  n_iter <- Nbchunks
  # Initializes the progress bar
  if (progressBar) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = n_iter, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "=")

  ## preparing a summary stat file
  stats <- array(0, c(length(unique(t2s$sample))+1, 5))
  rownames(stats) <- c("total", t2s$sample)
  colnames(stats) <- c("R1.R2.in", "R1.R2.shorter.than.primer", "NNNN", "R1.R2.out", "unique.reads")
  stats["total","R1.R2.in"] <- Nreads * 2
  stats["total","R1.R2.shorter.than.primer"] <- 0 ## for both R1 and R2

  ## width filtering by taking the longest of the two primers
  widthCutOff <- max(width(head(sread(primers), 1)), width(tail(sread(primers), 1)))

  #### collect headers of each sample, to check whether we have reads ascribed to multiple samples (not possible, only happen if we are too permissive with allowedMis)
  list_headers_sample <- c()
  list_samples        <- c()
  
  # dummy counter
  cpt <- 1
  repeat {

    fqR1IN <- yield(fqR1INstr)
    fqR2IN <- yield(fqR2INstr)
    if (length(fqR1IN) == 0) {
      break
    }

    if (!identical(sub(" .*", "", id(fqR1IN)), sub(" .*", "", id(fqR2IN)))) {
      message(libID, ": ", "Not the same order in the fastq files. Maybe check header of the fastq to parse correctly?")
      break
    }
    # sample names
    samples <- t2s$sample
    rownames(t2s) <- t2s$sample
    fwdCombs <- unique(t2s$forward)

    # make a named vector of it
    prim <- as.character(sread(primers))
    names(prim) <- as.character(id(primers))

    # we need to remove reads that are below the width of primer
    fqR1IN_ <- fqR1IN[width(fqR1IN) > widthCutOff & width(fqR2IN) > widthCutOff]
    # fetch the same reads in R2
    fqR2IN_ <- fqR2IN[width(fqR1IN) > widthCutOff & width(fqR2IN) > widthCutOff]
    ### at each chunk we remove the bad reads
    stats["total","R1.R2.shorter.than.primer"] <- stats["total","R1.R2.shorter.than.primer"] + (length(fqR1IN) - length(fqR1IN_)) + (length(fqR2IN) - length(fqR2IN_))

    for (fwdComb in fwdCombs)
    {
      #fwdComb <- fwdCombs[1]
      t2s_fwdComb <- subset(t2s, t2s$forward == fwdComb)
      # ge the fwd primer
      fwd <- prim[as.character(fwdComb)]
      if (is.na(fwd))
      {
        message(paste0(libID, ": The forward primer:", fwdComb, " does not exist in the fasta, skipping ", sample))
        next
      }

      ### searching the fwd primer in R1
      id.fwd.R1 <- vmatchPattern2(fwd, sread(narrow(fqR1IN_, start = 1, end = as.numeric(nchar(fwd)))), fixed=FALSE, 
                                  max.mismatch=allowedMis, with.indels = withIndels)
      # get the indexing
      id.fwd.R1.x <- elementNROWS(id.fwd.R1)
      tmpR1fwd         <- fqR1IN_[which(id.fwd.R1.x==1)]
      tmpR2fwd.for.rev <- fqR2IN_[which(id.fwd.R1.x==1)]

      ### searching the fwd primer in R2
      id.fwd.R2 <- vmatchPattern2(fwd, sread(narrow(fqR2IN_, start = 1, end = as.numeric(nchar(fwd)))), fixed=FALSE, 
                                  max.mismatch=allowedMis, with.indels = withIndels)
      # get the indexing
      id.fwd.R2.x <- elementNROWS(id.fwd.R2)
      tmpR2fwd         <- fqR2IN_[which(id.fwd.R2.x==1)]
      tmpR1fwd.for.rev <- fqR1IN_[which(id.fwd.R2.x==1)]

      # now for each reverse primer associated with the fwd
      for (revComb in t2s_fwdComb$reverse)
      {
        #revComb <- t2s_fwdComb$reverse[1]
        # get the rev seq and reverse complement it
        rev <- prim[as.character(revComb)]
        if (is.na(rev))
        {
          message(paste0(libID, ": The reverse primer:", revComb, " does not exist in the fasta, skipping ", sample))
          next
        }
        # sample name
        sample <- subset(t2s_fwdComb, t2s_fwdComb$reverse == revComb)[,"sample"]

        ### searching the rev primer in the reads of R2 that matched the R1fwd
        id.rev.R2 <- vmatchPattern2(rev, sread(narrow(tmpR2fwd.for.rev, start = 1, end = as.numeric(nchar(rev)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R2.x     <- elementNROWS(id.rev.R2)
        ## These are the reads that match the fwd in R1 and rev in R2 (but now switched orientation)
        tmpR2fwd.rev.ok <- tmpR2fwd.for.rev[which(id.rev.R2.x==1)]

        ### searching the rev primer in the reads of R1 that matched the R2fwd
        id.rev.R1 <- vmatchPattern2(rev, sread(narrow(tmpR1fwd.for.rev, start = 1, end = as.numeric(nchar(rev)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R1.x <- elementNROWS(id.rev.R1)
        ## These are the reads that match the rev in R1 and fwd in R2 (but now switched orientation)
        tmpR1fwd.rev.ok <- tmpR1fwd.for.rev[which(id.rev.R1.x==1)]

        ##### NOW reverse primer
        ## searching reverse in R1
        id.rev.R1 <- vmatchPattern2(rev, sread(narrow(fqR1IN_, start = 1, end = as.numeric(nchar(rev)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R1.x <- elementNROWS(id.rev.R1)
        tmpR1rev         <- fqR1IN_[which(id.rev.R1.x==1)]
        tmpR2rev.for.fwd <- fqR2IN_[which(id.rev.R1.x==1)]
        # searching forward in the reads of R2 that matched the R2rev
        id.rev.R2 <- vmatchPattern2(fwd, sread(narrow(tmpR2rev.for.fwd, start = 1, end = as.numeric(nchar(fwd)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R2.x     <- elementNROWS(id.rev.R2)
        ## These are the reads that start with the rev in R1 and fwd in R2 (but now switched orientation)
        ## --> expected 5' to 3' orientation are these two files
        tmpR2rev.fwd.ok <- tmpR2rev.for.fwd[which(id.rev.R2.x==1)]

        ## searching reverse in R2
        id.rev.R2 <- vmatchPattern2(rev, sread(narrow(fqR2IN_, start = 1, end = as.numeric(nchar(rev)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R2.x <- elementNROWS(id.rev.R2)
        tmpR2rev         <- fqR2IN_[which(id.rev.R2.x==1)]
        tmpR1rev.for.fwd <- fqR1IN_[which(id.rev.R2.x==1)]
        # searching forward in the reads of R2 that matched the R2rev
        id.rev.R2 <- vmatchPattern2(fwd, sread(narrow(tmpR1rev.for.fwd, start = 1, end = as.numeric(nchar(fwd)))), fixed=FALSE, 
                                    max.mismatch=allowedMis, with.indels = withIndels)
        # get the indexing
        id.rev.R2.x     <- elementNROWS(id.rev.R2)
        ## These are the reads that start with the rev in R2 and fwd in R1 (but now switched orientation)
        ## --> expected 5' to 3' orientation are these two files
        tmpR1rev.fwd.ok <- tmpR1rev.for.fwd[which(id.rev.R2.x==1)]

        ## Any NN? // fixing naming of files // to fix above at some point
        tmpR1fwd.rev.exp <- tmpR1rev.fwd.ok[nFilter(0)(tmpR1rev.fwd.ok) & nFilter(0)(tmpR2fwd.rev.ok)]
        tmpR2rev.fwd.exp <- tmpR2fwd.rev.ok[nFilter(0)(tmpR2fwd.rev.ok) & nFilter(0)(tmpR1rev.fwd.ok)]
        #
        tmpR1rev.fwd.exp <- tmpR1fwd.rev.ok[nFilter(0)(tmpR1fwd.rev.ok) & nFilter(0)(tmpR2rev.fwd.ok)]
        tmpR2fwd.rev.exp <- tmpR2rev.fwd.ok[nFilter(0)(tmpR2rev.fwd.ok) & nFilter(0)(tmpR1fwd.rev.ok)]

        ## stats
        ### at each chunk we add the reads
        stats[sample,"NNNN"] <- stats[sample,"NNNN"] + ((length(tmpR1fwd.rev.ok) - length(tmpR1fwd.rev.exp))) + ((length(tmpR1rev.fwd.ok) - length(tmpR1rev.fwd.exp))) + ((length(tmpR2fwd.rev.ok) - length(tmpR2fwd.rev.exp))) + ((length(tmpR2rev.fwd.ok) - length(tmpR2rev.fwd.exp)))
        stats[sample,"R1.R2.out"] <- stats[sample,"R1.R2.out"] + length(tmpR1fwd.rev.exp) + length(tmpR1rev.fwd.exp) + length(tmpR2fwd.rev.exp) + length(tmpR2rev.fwd.exp)

        ### collect header of R1fwd.rev
        list_headers_sample <- c(list_headers_sample, as.character(id(tmpR1fwd.rev.ok)))
        list_samples        <- c(list_samples, rep(sample, length(tmpR1fwd.rev.ok)))
        
        # and write the fastq (+ 1 to start just after the primer, "a" because we use the streamer)
        writeFastq(narrow(tmpR1fwd.rev.exp, start = nchar(fwd)+1, end = width(tmpR1fwd.rev.exp)), paste0(outputFolder,"/", sample, "_R1fwd.rev.fastq.gz"), "a")
        writeFastq(narrow(tmpR1rev.fwd.exp, start = nchar(rev)+1, end = width(tmpR1rev.fwd.exp)), paste0(outputFolder,"/", sample, "_R1rev.fwd.fastq.gz"), "a")
        #
        writeFastq(narrow(tmpR2fwd.rev.exp, start = nchar(fwd)+1, end = width(tmpR2fwd.rev.exp)), paste0(outputFolder,"/", sample, "_R2fwd.rev.fastq.gz"), "a")
        writeFastq(narrow(tmpR2rev.fwd.exp, start = nchar(rev)+1, end = width(tmpR2rev.fwd.exp)), paste0(outputFolder,"/", sample, "_R2rev.fwd.fastq.gz"), "a")
      }
    }
    if (progressBar) setTxtProgressBar(pb, cpt)
    cpt <- cpt + 1
  }

  close(fqR1INstr)
  close(fqR2INstr)
  ## export the stats
  stats[,"R1.R2.in"][2:nrow(stats)] <- NA
  stats[,"R1.R2.shorter.than.primer"][2:nrow(stats)] <- NA
  stats["total", "NNNN"] <- sum(stats[,"NNNN"][2:nrow(stats)])
  stats["total", "R1.R2.out"] <- sum(stats[,"R1.R2.out"][2:nrow(stats)])
  stats[,"unique.reads"] <- stats[,"R1.R2.out"] / 2
  stats <- data.frame(SampleID = rownames(stats), stats)
  write.table(stats, paste0(outputFolder, "/00_Stats_demultiplexing_", unique(t2s$run), ".tsv"), sep = "\t", quote = F, row.names = F)
  
  ## make a dataframe of headers to check
  header_df <- data.frame(samples = list_samples, headers = list_headers_sample)
  if(nrow(header_df) != length(unique(header_df$headers))) {
    ## export it 
    prob <- names(table(header_df$headers)[table(header_df$headers)>1])
    header_df$problematic <- ifelse(header_df$headers %in% prob, TRUE, FALSE)
    header_df <- header_df[order(header_df$headers),]
    header_df <- subset(header_df, header_df$problematic == T)
    write.table(header_df, paste0(outputFolder, "/01_Problematic_reads_", unique(t2s$run), ".tsv"), sep = "\t", quote = F, row.names = F)
    warning("There are ", length(prob), "/", stats["total","unique.reads"], " reads that have been ascribed to multiple samples, see ", 
            paste0("01_problematic_reads_", unique(t2s$run), ".tsv"), " file for details.\n--> You should reduce the 'allowedMis' parameter")
    
  }
  
  ## free memory
  rm(tmpR1fwd.rev.exp, tmpR1rev.fwd.exp, tmpR2fwd.rev.exp, tmpR2rev.fwd.exp)
  rm(tmpR1fwd.rev.ok, tmpR1rev.fwd.ok, tmpR2fwd.rev.ok, tmpR2rev.fwd.ok)
  rm( fqR1IN_, fqR2IN_, fqR1IN, fqR2IN)
  invisible(gc())
  if (progressBar) close(pb)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if (units(time.taken) == "secs") {
    time.taken_num <- round(as.numeric(end.time - start.time), 2)  ## in seconds
    message(libID, ": demultiplexed ", stats["total","unique.reads"], " reads (", round(stats["total","unique.reads"] *100 / (stats["total","R1.R2.in"]/2), 1)  ,"%). it took ", round(time.taken_num, 2), " seconds")
  } else {
    time.taken_num <- round(as.numeric(end.time - start.time), 2) / 60 ## in hours (after a minute, minutes is the default unit)
    message(libID, ": demultiplexed ", stats["total","unique.reads"], " reads (", round(stats["total","unique.reads"] *100 / (stats["total","R1.R2.in"]/2), 1)  ,"%). it took ", round(time.taken_num, 2), " hours")
  }
}
