


demultiplexNanopore <- function(primersFile, t2s, fastqIn, allowedMis = 0, trim = T, withIndels = F, outputFolder = "demult", 
                                overwrite = F, chunkSize = 100000, progressBar = T) {
  
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
    C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject, max.mismatch, 
                    min.mismatch, with.indels, fixed, algo, "MATCHES_AS_RANGES", PACKAGE="Biostrings")
    unlisted_ans <- IRanges(start=unlist(C_ans[[1L]], use.names=FALSE), width=unlist(C_ans[[2L]], use.names=FALSE))
    relist(unlisted_ans, C_ans[[1L]])
  }

  suppressMessages({
    library(ShortRead)
    library(Biostrings)
    library(dada2)
  })
  
  # tracking time
  start.time <- Sys.time()
  
  ## validate t2s
  if (!is.data.frame(t2s)) t2s <- read.table(t2s, header = T, sep=",", stringsAsFactors = F)
  libID <- t2s[1,1]
  
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
    stop(libID, ": some primer combinations are duplicated, please check and fix ", 
         names(table(paste0(t2s$forward, "_", t2s$reverse))[table(paste0(t2s$forward, "_", t2s$reverse))>1]))
  }
  
  
  ## reading both fastq by chunk // see here: http://bioconductor.org/packages/devel/bioc/vignettes/ShortRead/inst/doc/Overview.pdf
  fqINstr <- FastqStreamer(fastqIn, n = chunkSize)
  
  ##### progress and stats
  # how many chunks will we have?
  ## Nreads in the fastqs. Since 4 lines per read in a fastq file:
  Nreads <- as.numeric(system2("gunzip", args = paste0(fastqIn, " -c | wc -l"), stdout = TRUE)) / 4
  Nbchunks <- ceiling(Nreads / chunkSize)
  message(libID, ": ", Nreads, " reads to demultiplex. Will do ", Nbchunks, " chunk(s) of ", chunkSize," reads.")
  
  # progress bar
  n_iter <- Nbchunks # Number of iterations of the loop
  
  # Initializes the progress bar
  if (progressBar) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = n_iter, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "=")
  
  ## preparing a summary stat file
  stats <- array(0, c(length(unique(t2s$sample))+1, 5))
  rownames(stats) <- c("total", t2s$sample)
  colnames(stats) <- c("Reads.in", "Reads.shorter.than.primer", "Neg.width.artefact", "NNNN", "Reads.out")
  stats["total","Reads.in"] <- Nreads
  stats["total","Reads.shorter.than.primer"] <- 0 
  
  ## width filtering by taking the longest of the two primers
  widthCutOff <- max(width(head(sread(primers), 1)), width(tail(sread(primers), 1)))
  
  
  #### collect headers of each sample, to check whether we have reads ascribed to multiple samples 
  ## (not possible, only happen if we are too permissive with allowedMis)
  list_headers_sample <- c()
  list_samples        <- c()
  
  # dummy counter
  cpt <- 1
  repeat {
    
    fqIN <- yield(fqINstr)
    
    if (length(fqIN) == 0) {
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
    fqIN.clean <- fqIN[width(fqIN) > widthCutOff]
    ### at each chunk we remove the bad reads
    stats["total","Reads.shorter.than.primer"] <- stats["total","Reads.shorter.than.primer"] + (length(fqIN) - length(fqIN.clean))
    
    for (fwdComb in fwdCombs) {
      #fwdComb <- fwdCombs[1]
      t2s_fwdComb <- subset(t2s, t2s$forward == fwdComb)
      # ge the fwd primer
      fwd <- prim[as.character(fwdComb)]
      if (is.na(fwd)) {
        message(paste0(libID, ": The forward primer:", fwdComb, " does not exist in the fasta, skipping "))
        next
      }
      
      ### searching the fwd primer
      # if primer expected right at the begining: 
      #if (!nano) idx <- vmatchPattern2(fwd, sread(narrow(fqIN.clean, start = 1, end = as.numeric(nchar(fwd)))), fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
      # if primer can be somewhere after (e.g. as in nanopore data), search in the first half of the read
      idx <- vmatchPattern2(fwd, sread(narrow(fqIN.clean, start = 1, end = floor(as.numeric(width(fqIN.clean)/2)))), 
                            fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
      
      # get the indexing
      xx <- elementNROWS(idx)
      ## fastq containing the fwd primer only
      tmpFq.fwd1 <- fqIN.clean[which(xx==1)]
      
      ### check for fwd RC
      fwdRC <- dada2::rc(fwd)
      
      #if (!nano) idxRC <- vmatchPattern2(fwdRC, sread(narrow(fqIN.clean, start = 1, end = as.numeric(nchar(fwdRC)))), fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
      # if primer can be somewhere after (e.g. as in nanopore data), search in the first half of the read
      idxRC <- vmatchPattern2(fwdRC, sread(narrow(fqIN.clean, start = 1, end = floor(as.numeric(width(fqIN.clean)/2)))), 
                              fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
      
      # get the indexing
      xx <- elementNROWS(idxRC)
      ## fastq containing the fwd primer only
      tmpFq.fwd2 <- fqIN.clean[which(xx==1)]
      ## RC the fastq
      tmpFq.fwd2 <- Biostrings::reverseComplement(tmpFq.fwd2)
      
      ## concatenate the fastqs
      tmpFastqFWD <- append(tmpFq.fwd1, tmpFq.fwd2)
      ### reads are all oriented correctly now
      

      # now for each reverse primer associated with the fwd
      for (revComb in t2s_fwdComb$reverse) {
        #revComb <- t2s_fwdComb$reverse[1]
        # get the rev seq and reverse complement it
        rev <- prim[as.character(revComb)]
        ## RC of the reverse primers 
        revRC <- dada2::rc(rev)
        
        if (is.na(revRC)) {
          message(paste0(libID, ": The reverse primer:", revComb, " does not exist in the fasta, skipping "))
          next
        }
        # sample name
        sample <- subset(t2s_fwdComb, t2s_fwdComb$reverse == revComb)[,"sample"]

        # check for the reverse at the end
        #if (!nano) idy <- vmatchPattern2(revRC, sread(narrow(tmpFastqFWD, start = (as.numeric(width(tmpFastqFWD)) - as.numeric(nchar(rev))), end = as.numeric(width(tmpFastqFWD)))), fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
        # check in the second half of the read
        idy <- vmatchPattern2(revRC, sread(narrow(tmpFastqFWD, start = floor(as.numeric(width(tmpFastqFWD)/2))), 
                                           end = as.numeric(width(tmpFastqFWD))), fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
        
        # get the indexing
        yy <- elementNROWS(idy)
        ## fastq containing the fwd and rev primers in the expected orientation
        tmpFastqFWD.REV <- tmpFastqFWD[which(yy==1)]
        
        ## we filter reads with any N
        tmpFastqFWD.REV.clean <- tmpFastqFWD.REV[nFilter(0)(tmpFastqFWD.REV)]
        
        ## stats
        ### at each chunk we add the reads
        stats[sample,"NNNN"] <- stats[sample,"NNNN"] + ((length(tmpFastqFWD.REV) - length(tmpFastqFWD.REV.clean)))
        
        
        ## get the final indexing of fwd and rev primers to export. 
        idx <- vmatchPattern2(fwd, sread(narrow(tmpFastqFWD.REV.clean, start = 1, end = floor(as.numeric(width(tmpFastqFWD.REV.clean)/2)))),
                              fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
        idy <- vmatchPattern2(revRC, sread(narrow(tmpFastqFWD.REV.clean, start = floor(as.numeric(width(tmpFastqFWD.REV.clean)/2))), end = as.numeric(width(tmpFastqFWD.REV.clean))), 
                              fixed=FALSE, max.mismatch=allowedMis, with.indels = withIndels)
        
        ## in rare case, start and end would give a negative width (especially loose matches), whe check and remove the reads, as well as throwing a warning
        ## concat the indexes to make sure they match, and remove the reads with negative width
        idx.fwd <- as.data.frame(idx)
        idy.rev <- as.data.frame(idy)
        ##
        idx.fwd$inREV <- ifelse(idx.fwd$group %in% idy.rev$group, T, F)
        idy.rev$inFWD <- ifelse(idy.rev$group %in% idx.fwd$group, T, F)
        ## uniques IDs 
        unique.ids <- unique(idx.fwd$group[idx.fwd$inREV == T], idy.rev$group[idy.rev$inFWD == T])
        idx.clean <- subset(idx.fwd, idx.fwd$group %in% unique.ids)
        idy.clean <- subset(idy.rev, idy.rev$group %in% unique.ids)
        ## now non neg width
        non.neg <- ifelse(idy.clean[,"end"] - idx.clean[,"start"] >= 0, T, F)
        ## throw an error if any 
        if (!is.na(table(non.neg)["FALSE"])) {
          stats[sample,"Neg.width.artefact"] <- as.numeric(table(non.neg)["FALSE"])
          neg.artefact <- tmpFastqFWD.REV.clean[unique.ids[which(non.neg == F)]]
          writeFasta(neg.artefact, file = paste0(outputFolder,"/", sample, "_neg.width.artefact.fasta"))
          message(as.numeric(table(non.neg)["FALSE"]), " read(s) in sample ", sample, " has primers matching at positions yielding a negative width, exporting a fasta and discarding from output.")
        }
        ## export fastq in good order 
        exp <- unique.ids[which(non.neg == T)]
        idx.clean <- idx.clean[idx.clean$group %in% exp,]
        idy.clean <- idy.clean[idy.clean$group %in% exp,]
        
        stats[sample,"Reads.out"] <- stats[sample,"Reads.out"] + length(tmpFastqFWD.REV.clean[exp])
        
        ## some negative start if no indel allowed (not full primers in the reads I guess)
        if (!is.na(table(idx.clean[,"start"] <= 0)["TRUE"])) {
          message(sample, " has some negative values for forward primer match, not full primer is detected in the read I guess. Exporting from 1st nucleotide" ) 
          idx.clean[,"start"][idx.clean[,"start"] <= 0] <- 1
        }
        if (!is.na(table(idy.clean[,"end"] > width(tmpFastqFWD.REV.clean[exp]))["TRUE"])) {
          message(sample, " has some negative values for forward primer match, not full primer is detected in the read I guess. Exporting until end of read" ) 
          idy.clean[,"end"][idy.clean[,"end"] > width(tmpFastqFWD.REV.clean[exp])] <- width(tmpFastqFWD.REV.clean[exp])
        }
        
        
        
        ### collect header
        list_headers_sample <- c(list_headers_sample, as.character(id(tmpFastqFWD.REV.clean)))
        list_samples        <- c(list_samples, rep(sample, length(tmpFastqFWD.REV.clean)))

        if (!trim) writeFastq(narrow(tmpFastqFWD.REV.clean[exp], start = idx.clean[,"start"], end = idy.clean[,"end"]), paste0(outputFolder,"/", sample, "_untrimmed.fastq.gz"), "a")
        if (trim)  writeFastq(narrow(tmpFastqFWD.REV.clean[exp], start = idx.clean[,"end"]+1, end = idy.clean[,"start"] -1), paste0(outputFolder,"/", sample, ".fastq.gz"), "a")
        
        
      }
    }
    if (progressBar) setTxtProgressBar(pb, cpt)
    cpt <- cpt + 1
  }

  close(fqINstr)
  
  ## export the stats
  stats[,"Reads.in"][2:nrow(stats)] <- NA
  stats[,"Reads.shorter.than.primer"][2:nrow(stats)] <- NA
  stats["total", "NNNN"] <- sum(stats[,"NNNN"][2:nrow(stats)])
  stats["total", "Reads.out"] <- sum(stats[,"Reads.out"][2:nrow(stats)])

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
  rm(fqIN, fqIN.clean, tmpFastqFWD, tmpFastqFWD.REV, tmpFastqFWD.REV.clean, tmpFq.fwd1, tmpFq.fwd2)
  invisible(gc())
  if (progressBar) close(pb)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if (units(time.taken) == "secs") {
    time.taken_num <- round(as.numeric(end.time - start.time), 2)  ## in seconds
    message(libID, ": demultiplexed ", stats["total","Reads.in"], " reads (", round(stats["total","Reads.out"] *100 / stats["total","Reads.in"], 1)  ,"% reads matched primers combinations). it took ", round(time.taken_num, 2), " seconds")
  } else {
    time.taken_num <- round(as.numeric(end.time - start.time), 2) / 60 ## in hours (after a minute, minutes is the default unit)
    message(libID, ": demultiplexed ", stats["total","Reads.in"], " reads (", round(stats["total","Reads.out"] *100 / stats["total","Reads.in"], 1)  ,"% reads matched primers combinations). it took ", round(time.taken_num, 2), " hours")
  }
}
