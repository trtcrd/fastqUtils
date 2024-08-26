
dada2workflow <- function(params, t2s, inputFolder = NULL, outputFolder = NULL, overwrite = F, novaseq = T, minLen = 60, cores_per_job = 4) {

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
    fwdFiles <- sort(list.files(path.in, pattern="_fwd", full.names = TRUE))
    revFiles <- sort(list.files(path.in, pattern="_rev", full.names = TRUE))
    ##

    # Extract sample names // Be carefull to check before going through the next step, but should be fine.
    sample.names <- sapply(strsplit(basename(fwdFiles), "_fwd"), `[`,1)

    ### creating directories
    path.tmp1 <- file.path(path.exp, "tmp1")
    if(!dir.exists(path.tmp1)) dir.create(path.tmp1)
    path.tmp2 <- file.path(path.exp, "tmp2")
    if(!dir.exists(path.tmp2)) dir.create(path.tmp2)
    path.rds <- file.path(path.exp, "output_rds")
    if(!dir.exists(path.rds)) dir.create(path.rds)

    ### starting the workflow
    filtFwdFiles <- file.path(path.tmp1, basename(fwdFiles))
    filtRevFiles <- file.path(path.tmp1, basename(revFiles))
    ## fix path.in to path.out
    filtFwdFiles <- gsub(path.in, path.exp, filtFwdFiles)
    filtRevFiles <- gsub(path.in, path.exp, filtRevFiles)
    # basic filtering with default DADA2 settings
    filtering <- filterAndTrim(fwd=fwdFiles, rev=revFiles, filt=filtFwdFiles, filt.rev=filtRevFiles,
                               multithread=cores_per_job, verbose=F, minLen = minLen) ###, orient.fwd=FWD)

    ## stats
    rownames(filtering) <- gsub("_fwd.fastq.gz", "", rownames(filtering))

    ## gather stats
    readsCounts <- data.frame(array(NA, c(nrow(filtering), 7)))
    colnames(readsCounts) <- c("reads.in", "reads.out", "orient", "cutadapt", "primer.free", "length", "dada2")
    rownames(readsCounts) <- sample.names
    readsCounts[,1] <- filtering[,1]
    readsCounts[,2] <- filtering[,2]

    ##keep only samples with reads
    keep <- filtering[,"reads.out"] > 100

    ### path update
    filtFwdFiles_kept <- filtFwdFiles[keep]
    filtRevFiles_kept <- filtRevFiles[keep]
    #

    # gather infos
    ## gather last error estimates after the 10 rounds
    lastE <- array(NA, c(length(filtFwdFiles_kept), 4))
    colnames(lastE) <- c("sample_id", "lastE_fwd", "lastE_rev", "duration")

    sampleNames <- sapply(strsplit(basename(filtFwdFiles_kept), "_fwd"), `[`, 1)

    ############################################################
    # error model to be adjusted for novaseq vs miseq (or binned vs non binned)
    if (novaseq)  errFun <- loessErrfun_mod
    if (!novaseq) errFun <- loessErrfun
    ############################################################

    errM.FWD <- learnErrors(filtFwdFiles_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    errM.REV <- learnErrors(filtRevFiles_kept, multithread=cores_per_job, randomize=TRUE, errorEstimationFunction = errFun, verbose = F)
    ## save the error model for both fwd and rev
    saveRDS(errM.FWD, file.path(path.rds, paste0("Mod_", libID, "_errM.FWD.rds")))
    saveRDS(errM.REV, file.path(path.rds, paste0("Mod_", libID, "_errM.REV.rds")))
    ##
    ## and the plot
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.FWD.pdf")))
    print(plotErrors(errM.FWD, nominalQ = T))
    dev.off()
    pdf(file.path(path.rds, paste0("Mod_plot_", libID, "_errM.REV.pdf")))
    print(plotErrors(errM.REV, nominalQ = T))
    dev.off()

    # derep now
    drp.FWD <- derepFastq(filtFwdFiles_kept)
    drp.REV <- derepFastq(filtRevFiles_kept)
    #

    ## denoising
    dd.FWD <- dada(drp.FWD, err=errM.FWD, selfConsist=F, multithread=cores_per_job, verbose = F)
    dd.REV <- dada(drp.REV, err=errM.REV, selfConsist=F, multithread=cores_per_job, verbose = F)
    #
    # trimOverhang is important here as we have fully overlapping reads
    merger <- mergePairs(dd.FWD, drp.FWD, dd.REV, drp.REV, trimOverhang=TRUE, verbose = F)
    sta <- makeSequenceTable(merger)
    rownames(sta) <- gsub("_fwd.fastq.gz", "", rownames(sta))
    #
    # save it
    saveRDS(sta, file.path(path.rds, paste0("ASVtable_", libID, "_merger.rds")))
    ## get the last errors estimates if not converged
    lastE[,"sample_id"] <- sampleNames
    lastE[,"lastE_fwd"] <- tail(dada2:::checkConvergence(errM.FWD),1)
    lastE[,"lastE_rev"] <- tail(dada2:::checkConvergence(errM.REV),1)

    ## mapping the number of reads
    readsCounts[names(keep[keep==T]), "dada2"] <- rowSums(sta[names(keep[keep==T]),])

    ## export the cleaning stats
    write.table(readsCounts, file = paste0(path.rds, "/00_readsCounts_stats_", libID, ".tsv"), quote = F, sep="\t", row.names = T, fileEncoding = "UTF-8", col.names = NA)

    ### export the stats
    write.table(lastE, file = paste0(path.rds, "/01_dada2_stats_",libID,".tsv"), quote = F, sep="\t", row.names = T, col.names = NA, fileEncoding = "UTF-8")


  }

}
