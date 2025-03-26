
#### Loop over libraries for demultiplexing
demultRun <- function(params, t2s, accountMixedOrientation = F, allowedMis = 0, outputFolder = NULL, overwrite = F, cores = 4, splitHeader = " ", chunkSize = 1000000) {

  ##
  suppressMessages({
    library(doMC)
    library(seqinr)
    library(ShortRead)
    library(Biostrings)
  })
  ## Choose how many cores to use
  registerDoMC(cores = cores)
  OK <- T ## ok to run the %dopar% loop?

  ## by default outputFolder is the "demultiplex_" t2s name
  if (is.null(outputFolder)) outputFolder <- paste0("run_", gsub(".csv", "", t2s), "_demultiplex")

  ## validate params file (primers, R1R2, output folder) and t2s
  params <- read.table(params, header = T, sep=",", stringsAsFactors = F)
  ## check lib names
  if (length(unique(params$library_name)) != length(params$library_name))  {
    OK <- F
    stop("some library names are identical in params file. please fix ", paste(names(table(params$library_name)[table(params$library_name) > 1]), collapse = " "))
  }
  rownames(params) <- params$library_name

  ## check if both R1 R2 exists in the path given
  missing_files <- c()
  for (i in params$library_name) {
    if (!file.exists(file = params[i,"R1"])) missing_files <- c(missing_files, params[i,"R1"])
    if (!file.exists(file = params[i,"R2"])) missing_files <- c(missing_files, params[i,"R2"])
  }
  if (length(missing_files)>0) {
    OK <- F
    stop(" These files are not found in the path given, please check:\n", paste(missing_files, collapse = "\n"))
  }

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
      stop(paste0(params$R1[i], " is set for both R1 and R2. please fix ", params$library_name[i]))
    }
  }

  ### t2s and primers file
  t2s <- read.table(t2s, header = T, sep=",", stringsAsFactors = F)
  ## check colnames in params file and t2s
  if (!all(colnames(t2s) %in% c("run","sample","forward","reverse"))) {
    OK <- F
    stop("please fix columns name of t2s file to get: 'run', 'sample', 'forward', 'reverse'. stop")
  }
  if (!all(colnames(params) %in% c("library_name","primers","R1","R2"))) {
    OK <- F
    stop("please fix columns name of params file file to get: 'library_name', 'primers', 'R1', 'R2'. stop")
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
    stop("some libraries in params file are not in the t2s ", i,". please fix ", spls)
    OK <- F
  }

  if (FALSE %in% names(table(unique(t2s$run) %in% unique(params$library_name)))) {
    spls <- paste(unique(t2s$run)[!unique(t2s$run) %in% unique(params$library_name)], " ")
    message("be aware that some libraries in t2s file are not in the params file. please check ", spls)
  }

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

  ## for each lib, this is done within the functions

  #################
  if(OK) {
    out <- foreach(i = 1:nrow(params)) %dopar% {
      # get exp folder and t2s
      path.exp.lib <- paste0(outputFolder, "/", params$library_name[i])
      t2s_tmp <- subset(t2s, t2s$run == params$library_name[i])
      lib <- params$library_name[i]
      # run the demultiplexer
      if (accountMixedOrientation) {
        demultiplexDoublePaired4files(primersFile= params[lib, "primers"],
                                      t2s = t2s_tmp,
                                      fastqR1In = params[lib, "R1"],
                                      fastqR2In = params[lib, "R2"],
                                      outputFolder = path.exp.lib,
                                      progressBar =F)
      } else if (!accountMixedOrientation) {
        demultiplexDoublePaired(primersFile= params[lib, "primers"],
                                t2s = t2s_tmp,
                                fastqR1In = params[lib, "R1"],
                                fastqR2In = params[lib, "R2"],
                                outputFolder = path.exp.lib,
                                progressBar =F)
      } else { stop("stop because accountMixedOrientation is not set. do you want to account for mixed read orientation (4 files) or not (2 files, with swapped reads in R1-R2)?") }
      message(params$library_name[i], ": done with demultiplexing. now getting tag distrib plot")
      # run the tag distrib
      taggedPrimersDistrib(primersFile= params[lib, "primers"],
                           t2s = t2s_tmp,
                           fastqR1In = params[lib, "R1"],
                           fastqR2In = params[lib, "R2"],
                           outputFolder = paste0(outputFolder, "/"),
                           progressBar = F)
    }
  }
}
