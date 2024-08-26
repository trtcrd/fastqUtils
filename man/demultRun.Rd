\encoding{UTF-8}
\name{demultRun}
\alias{demultRun}
\title{Demultiplexing a list of illumina amplicon libraries}

\description{
Demultiplexing a list of illumina amplicon libraries in parallel.
}


\usage{
  demultRun(params, t2s, accountMixedOrientation = F, allowedMis = 0, outputFolder = NULL, overwrite = F, cores = 4, splitHeader = " ", chunkSize = 1000000)
}


\arguments{
  \item{params}{A .csv file containing the instructions to demultiplex a list of libraries. Must contains column names "library_name", "primers" (fasta containing primers names and sequences), "R1", "R2" where primers, R1 and R2 are paths to the files.}
  \item{t2s}{A "tag-to-sample" .csv file containing all informations for demultiplexing samples. Must contains columns names "run" (for library name, to.do: harmonize this field across functions), "sample" (sample name, cannot be duplicated in a given library), "forward" and "reverse" (for the tagged primers combination).}
  \item{accountMixedOrientation}{Should we account for mixed orientation in the R1 and R2 (default = F)? This happens with pcr-free library preparation protocols, i.e. orientation of reads is not controlled. If TRUE, 4 files will be produced from the R1 and R2 files: R1fwd, R21rev, R2fwd, R2rev.}
  \item{allowedMis}{Number of mismatches allowed in searching primers sequences in the reads (default = 0).}
  \item{outputFolder}{Output folder. Can be specified. Default is "run_", t2s_name, "_demultiplex"}
  \item{overwrite}{If outputFolder already exists, will throw an error. Shall we overwrite the folder (default = FALSE)?}
  \item{cores}{CPU to be used. If more than one, multiple libraries will be demultiplexed in parallel (Default = 4).}
  \item{splitHeader}{Character for parsing the fastq header names and checking consistency between R1 and R2 files (default = " "). Most likely not to be edited, unless a specific error about that is thrown.}
  \item{chunkSize}{How many reads to be processed per batch (default = 1000000). The higher the number, the more RAM is needed. Default involves about 1 Go per library}
}

\details{
The R1 and R2 fastqs files are processed by chunks, and multiplexed libraries can be processed in parallel (one library per thread). Can account for mixed orientation of reads. This function is a wrapper of the \code{demultiplexDoublePaired} (if not accounting fior read orientation) and the \code{demultiplexDoublePaired4files} (if accounting for read orientation) functions. See \href{https://github.com/trtcrd/SLIM?tab=readme-ov-file#prepare-and-upload-your-data/}{here} for an example of \code{t2s} and \code{primers} files format). The \code{taggedPrimersDistrib} function is also run on each library.
}

\examples{
demultRun(params = "my_run_params.csv", t2s = "t2s_my_run.csv", cores = 12)
}
