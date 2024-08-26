
dada2workflow4files <- function(params, t2s, inputFolder = NULL, outputFolder = NULL, overwrite = F, novaseq = T, minLen = 60, cores_per_job = 4) {

  ##
  suppressMessages({
    library(doMC)
    library(dada2)
    library(seqinr)
    library(ShortRead)
    library(Biostrings)
  })

  ## Choose how many cores to use
  registerDoMC(cores = floor(cores_per_job / 2)) ## make possibly (cores_per_job / 2 * cores_per_job) cores used at the same time (cores_per_job = 4 ---> 2 parallel dada workflows with 4 each = 8 cores)

  ## validate params file (primers, R1R2, output folder) and t2s
  params <- read.table(params, header = T, sep=",", stringsAsFactors = F)
  ## check lib names
  if (length(unique(params$library_name)) != length(params$library_name))  {
    OK <- F
    stop("some library names are identical in params file. please fix ", paste(names(table(params$library_name)[table(params$library_name) > 1]), collapse = " "))
  }
  rownames(params) <- params$library_name

  ## check R1 R2
  if (length(unique(params$R1)) != length(unique(params$R2)))  {
    if (length(unique(params$R1)) < length(unique(params$R2))) n <- paste(params$R1[table(params$R1)[table(params$R1) > 1]], collapse = " ")
    if (length(unique(params$R1)) > length(unique(params$R2))) n <- paste(table(params$R1)[table(params$R2) > 1], collapse = " ")
    OK <- F
    stop("some identical files in R1 or R2. please fix ", n)
  }

  for (i in 1:nrow(params)) {
    if (identical(params$R1[i], params$R2[i]))  {
      OK <- F
      stop(paste0(params$R1[i], " is the set as R1 and R2. please fix ", params$library_name[i]))
    }
  }

  ### t2s and primers file
  t2s_name <- gsub(".csv", "", t2s)
  t2s <- read.table(t2s, header = T, sep=",", stringsAsFactors = F)
  ## check colnames in t2s
  if (!all(colnames(t2s) %in% c("run","sample","forward","reverse"))) {
    OK <- F
    stop("please fix columns name of t2s file to get: 'run', 'sample', 'forward', 'reverse'. stop")
  }
  ## check for duplicate
  for (i in unique(params$library_name)) {
    t2s_tmp <- subset(t2s, t2s$run == i)
    if (length(t2s_tmp$sample) != length(unique(t2s_tmp$sample))) {
      spls <- paste(names(table(t2s_tmp$sample)[table(t2s_tmp$sample) > 1]), collapse = " ")
      OK <- F
      stop("there are duplicate samples in the tag to sample of ", i,". please fix ", spls)
    }
    if (length(table(paste0(t2s_tmp$forward, "_", t2s_tmp$reverse))) < nrow(t2s_tmp)) {
      comb <- names(which(table(paste0(t2s_tmp$forward, ",", t2s_tmp$reverse))>1))
      s1 <- sapply(strsplit(comb, ","), `[`, 1)
      s2 <- sapply(strsplit(comb, ","), `[`, 2)
      spls <- t2s_tmp[t2s_tmp$forward %in% s1 & t2s_tmp$reverse %in% s2, "sample"]
      OK <- F
      stop("some primer combinations are duplicated in ",i,". please fix ", paste(spls, collapse=" "))
    }
    ## check primer file
    primers <- readFasta(params[i, "primers"])
  }


  ## check if lib names in params file are all in t2s // warning only if other way around
  if (FALSE %in% names(table(unique(params$library_name) %in% unique(t2s$run)))) {
    spls <- paste(unique(params$library_name)[!unique(params$library_name) %in% unique(t2s$run)], " ")
    OK <- F
    stop("some libraries in params file are not in the t2s ", i,". please fix ", spls)
  }

  if (FALSE %in% names(table(unique(t2s$run) %in% unique(params$library_name)))) {
    spls <- paste(unique(t2s$run)[!unique(t2s$run) %in% unique(params$library_name)], " ")
    message("be aware that some libraries in t2s file are not in the params file. please check ", spls)
  }

  ## by default outputFolder is the "dada2_" t2s name
  if (is.null(outputFolder)) outputFolder <- paste0("run_", t2s_name, "_dada2")
  ## check whether the run output folder exist (by default the names of the t2s, can be specified)
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


  ### Error model inference with binned NovaSeq basecalling
  # see: https://github.com/benjjneb/dada2/issues/938
  # https://github.com/benjjneb/dada2/issues/1307
  loessErrfun_mod <- function (trans) {
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow = 0, ncol = length(qq))
    for (nti in c("A", "C", "G", "T")) {
      for (ntj in c("A", "C", "G", "T")) {
        if (nti != ntj) {
          errs <- trans[paste0(nti, "2", ntj), ]
          tot <- colSums(trans[paste0(nti, "2", c("A","C","G","T")),])
          rlogp <- log10((errs + 1)/tot)
          rlogp[is.infinite(rlogp)] <- NA
          df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
          mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
          pred <- predict(mod.lo, qq)
          maxrli <- max(which(!is.na(pred)))
          minrli <- min(which(!is.na(pred)))
          pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
          pred[seq_along(pred) < minrli] <- pred[[minrli]]
          est <- rbind(est, 10^pred)
        }
      }
    }
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-07
    est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
    err <- rbind(1 - colSums(est[1:3,]), est[1:3,], est[4,], 1 - colSums(est[4:6,]), est[5:6,], est[7:8,],
                 1 - colSums(est[7:9,]), est[9,], est[10:12,], 1 - colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),"2", c("A", "C", "G", "T"))
    colnames(err) <- colnames(trans)
    return(err)
  }

  ## for each library in the params file
  out <- foreach(i = 1:nrow(params)) %dopar% {
    # get path folder and t2s
    libID <- params$library_name[i]
    if (is.null(inputFolder)) path.in <- paste0(getwd(), "/run_", t2s_name, "_demultiplex/", libID, "/")
    if (!is.null(inputFolder)) path.in <- paste0(getwd(), "/", inputFolder, "/", libID, "/")
    # path exp
    path.exp <- file.path(getwd(), outputFolder, libID)
    if(!dir.exists(path.exp)) dir.create(path.exp)
    t2s_tmp <- subset(t2s, t2s$run == libID)

    ## path input of the library to process

    ### first path for the fastq files
    # R1fwd.rev and R2rev.fwd
    R1fwd.rev <- sort(list.files(path.in, pattern="_R1fwd.rev", full.names = TRUE))
    R2rev.fwd <- sort(list.files(path.in, pattern="_R2rev.fwd", full.names = TRUE))
    # R1rev.fwd and R2fwd.rev
    R1rev.fwd <- sort(list.files(path.in, pattern="_R1rev.fwd", full.names = TRUE))
    R2fwd.rev <- sort(list.files(path.in, pattern="_R2fwd.rev", full.names = TRUE))
    ##

    # Extract sample names // Be carefull to check before going through the next step, but should be fine.
    sample.names <- sapply(strsplit(basename(R1fwd.rev), "_R1fwd.rev"), `[`,1)

    ### creating directories
    path.tmp1 <- file.path(path.exp, "tmp1")
    if(!dir.exists(path.tmp1)) dir.create(path.tmp1)
    path.tmp2 <- file.path(path.exp, "tmp2")
    if(!dir.exists(path.tmp2)) dir.create(path.tmp2)
    path.rds <- file.path(path.exp, "output_rds")
    if(!dir.exists(path.rds)) dir.create(path.rds)

    ### starting the workflow FWD
    filtR1fwd.rev <- file.path(path.tmp1, basename(R1fwd.rev))
    filtR2rev.fwd <- file.path(path.tmp1, basename(R2rev.fwd))
    ## fix path.in to path.out
    filtR1fwd.rev <- gsub(path.in, path.exp, filtR1fwd.rev)
    filtR2rev.fwd <- gsub(path.in, path.exp, filtR2rev.fwd)
    # basic filtering with default DADA2 settings
    filteringFWD <- filterAndTrim(fwd=R1fwd.rev, rev=R2rev.fwd, filt=filtR1fwd.rev, filt.rev=filtR2rev.fwd,
                                  multithread=cores_per_job, verbose=F, minLen = minLen) ###, orient.fwd=FWD)
    ### starting the workflow REV
    filtR1rev.fwd <- file.path(path.tmp1, basename(R1rev.fwd))
    filtR2fwd.rev <- file.path(path.tmp1, basename(R2fwd.rev))
    # basic filtering with default DADA2 settings
    filteringREV <- filterAndTrim(fwd=R1rev.fwd, rev=R2fwd.rev, filt=filtR1rev.fwd, filt.rev=filtR2fwd.rev,
                                  multithread=cores_per_job, verbose=F, minLen = minLen) ###, orient.fwd=FWD)

    ## summing stats
    filtering <- filteringFWD + filteringREV
    rownames(filtering) <- gsub("_R1fwd.rev.fastq.gz", "", rownames(filtering))

    ## gather stats
    readsCounts <- data.frame(array(NA, c(nrow(filtering), 7)))
    colnames(readsCounts) <- c("reads.in", "reads.out", "orient", "cutadapt", "primer.free", "length", "dada2")
    rownames(readsCounts) <- sample.names
    readsCounts[,1] <- filtering[,1]
    readsCounts[,2] <- filtering[,2]

    ##keep only samples with reads
    keep <- filtering[,"reads.out"] > 100

    ### path update
    filtR1fwd.rev_kept <- filtR1fwd.rev[keep]
    filtR2rev.fwd_kept <- filtR2rev.fwd[keep]
    #
    filtR1rev.fwd_kept <- filtR1rev.fwd[keep]
    filtR2fwd.rev_kept <- filtR2fwd.rev[keep]

    # gather infos
    ## gather last error estimates after the 10 rounds
    lastE <- array(NA, c(length(filtR1fwd.rev_kept), 6))
    colnames(lastE) <- c("sample_id", "lastE_R1fwd.rev", "lastE_R2rev.fwd", "lastE_R1rev.fwd", "lastE_R2fwd.rev", "duration")

    sampleNames <- sapply(strsplit(basename(filtR1fwd.rev_kept), "_R1fwd.rev"), `[`, 1)

    ############################################################
    # error model to be adjusted for novaseq vs miseq (or binned vs non binned)
    if (novaseq)  errFun <- loessErrfun_mod
    if (!novaseq) errFun <- loessErrfun
    ############################################################

    errM.R1fwd.rev <- learnErrors(filtR1fwd.rev_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    errM.R2rev.fwd <- learnErrors(filtR2rev.fwd_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    ## save the error model for both fwd and rev
    saveRDS(errM.R1fwd.rev, file.path(path.rds, paste0("Mod_", libID, "_errM.R1fwd.rev.rds")))
    saveRDS(errM.R2rev.fwd, file.path(path.rds, paste0("Mod_", libID, "_errM.R2rev.fwd.rds")))
    ##
    errM.R1rev.fwd <- learnErrors(filtR1rev.fwd_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    errM.R2fwd.rev <- learnErrors(filtR2fwd.rev_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    ## save the error model for both fwd and rev
    saveRDS(errM.R1rev.fwd, file.path(path.rds, paste0("Mod_", libID, "_errM.R1rev.fwd.rds")))
    saveRDS(errM.R2fwd.rev, file.path(path.rds, paste0("Mod_", libID, "_errM.R2fwd.rev.rds")))

    ## and the plot
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.R1fwd.rev.pdf")))
    print(plotErrors(errM.R1fwd.rev, nominalQ = T))
    dev.off()
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.R2rev.fwd.pdf")))
    print(plotErrors(errM.R2rev.fwd, nominalQ = T))
    dev.off()
    #
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.R1rev.fwd.pdf")))
    print(plotErrors(errM.R1rev.fwd, nominalQ = T))
    dev.off()
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.R2fwd.rev.pdf")))
    print(plotErrors(errM.R2fwd.rev, nominalQ = T))
    dev.off()

    # derep now
    drp.R1fwd.rev <- derepFastq(filtR1fwd.rev_kept)
    drp.R2rev.fwd <- derepFastq(filtR2rev.fwd_kept)
    #
    drp.R1rev.fwd <- derepFastq(filtR1rev.fwd_kept)
    drp.R2fwd.rev <- derepFastq(filtR2fwd.rev_kept)

    ## denoising
    dd.R1fwd.rev <- dada(drp.R1fwd.rev, err=errM.R1fwd.rev, selfConsist=F, multithread=cores_per_job, verbose = F)
    dd.R2rev.fwd <- dada(drp.R2rev.fwd, err=errM.R2rev.fwd, selfConsist=F, multithread=cores_per_job, verbose = F)
    #
    dd.R1rev.fwd <- dada(drp.R1rev.fwd, err=errM.R1rev.fwd, selfConsist=F, multithread=cores_per_job, verbose = F)
    dd.R2fwd.rev <- dada(drp.R2fwd.rev, err=errM.R2fwd.rev, selfConsist=F, multithread=cores_per_job, verbose = F)
    #

    # trimOverhang is important here as we have fully overlapping reads
    merger.FWD <- mergePairs(dd.R1fwd.rev, drp.R1fwd.rev, dd.R2rev.fwd, drp.R2rev.fwd, trimOverhang=TRUE, verbose = F)
    sta.FWD <- makeSequenceTable(merger.FWD)
    rownames(sta.FWD) <- gsub("_R1fwd.rev.fastq.gz", "", rownames(sta.FWD))
    #
    merger.REV <- mergePairs(dd.R1rev.fwd, drp.R1rev.fwd, dd.R2fwd.rev, drp.R2fwd.rev, trimOverhang=TRUE, verbose= F)
    sta.REV <- makeSequenceTable(merger.REV)
    rownames(sta.REV) <- gsub("_R1rev.fwd.fastq.gz", "", rownames(sta.REV))

    ## finally merging the table in same orientation and summing reads per samples
    ## tryRC to make sure they are merged !
    sta <- mergeSequenceTables(sta.FWD, sta.REV, tryRC = T, repeats = "sum")
    #colnames(sta) <- rc(colnames(sta))

    # save it
    saveRDS(sta, file.path(path.rds, paste0("ASVtable_", libID, "_merger.rds")))
    ## get the last errors estimates if not converged
    lastE[,"sample_id"] <- sampleNames
    lastE[,"lastE_R1fwd.rev"] <- tail(dada2:::checkConvergence(errM.R1fwd.rev),1)
    lastE[,"lastE_R2rev.fwd"] <- tail(dada2:::checkConvergence(errM.R2rev.fwd),1)
    lastE[,"lastE_R1rev.fwd"] <- tail(dada2:::checkConvergence(errM.R1rev.fwd),1)
    lastE[,"lastE_R2fwd.rev"] <- tail(dada2:::checkConvergence(errM.R2fwd.rev),1)

    ## mapping the number of reads
    readsCounts[names(keep[keep==T]), "dada2"] <- rowSums(sta[names(keep[keep==T]),])

    ## export the cleaning stats
    write.table(readsCounts, file = paste0(path.rds, "/00_readsCounts_stats_", libID, ".tsv"), quote = F, sep="\t", row.names = T, fileEncoding = "UTF-8", col.names = NA)

    ### export the stats
    write.table(lastE, file = paste0(path.rds, "/01_dada2_stats_",libID,".tsv"), quote = F, sep="\t", row.names = T, col.names = NA, fileEncoding = "UTF-8")


  }

}
