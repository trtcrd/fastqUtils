\encoding{UTF-8}
\name{dada2workflow}
\alias{dada2workflow}
\title{dada2 workflow for demultiplexed libraries}

\description{
Runs the default dada2 workflow demultiplexed libraries (where possible mixed orientation of reads was not accounted for).
}

\usage{
  dada2workflow(params, t2s, inputFolder = NULL, outputFolder = NULL,
    overwrite = F, novaseq = T, minLen = 60, cores_per_job = 4)
}


\arguments{
  \item{params}{A .csv file containing the instructions to demultiplex a list of libraries. Must contains column names "library_name", "primers" (fasta containing primers names and sequences), "R1", "R2" where primers, R1 and R2 are paths to the files.}
  \item{t2s}{A "tag-to-sample" .csv file containing all informations for demultiplexing samples. Must contains columns names "run" (for library name, to.do: harmonize this field across functions), "sample" (sample name, cannot be duplicated in a given library), "forward" and "reverse" (for the tagged primers combination). }
  \item{inputFolder}{By default, will look into the default \code{outputFolder} of the \code{demultRun} function. Can be manually set}
  \item{outputFolder}{Output folder. Can be specified. Default is "run_", t2s_name, "_dada2"}
  \item{overwrite}{If outputFolder already exists, will throw an error. Shall we overwrite the folder (default = FALSE)?}
  \item{novaseq}{If fastqs contains binned Q-scores (see \href{https://github.com/benjjneb/dada2/issues/938}{here} and \href{https://github.com/benjjneb/dada2/issues/1307}{here} for more information), we use a modified version of the \code{loessErrfun()} function in dada2 (shared by Guillem Salazar \href{https://github.com/benjjneb/dada2/issues/938#issuecomment-774051061}{here}).}
  \item{minLen}{\code{minLen} parameter of the \code{filterAndTrim} function of the dada2 workflow.}
  \item{cores_per_job}{CPU dedicated to each dada2 workflow (Default = 4). This is scaled with the number of demultiplexed libraries to be processed (by default: cores_per_job/2 libraries can be processed in parallel, each with cores_per_job cores).}
}

\details{
Runs the default dada2 workflow on each demultiplexed library. The two fastqs files fwd and rev are processed and merged. Check whether you fastqs contains binned scores (most likely with NovaSeq data, and most likely not if MiSeq data). The workflow exports an ASV table for each library, but does not filter chimeras, nor taxonomically annotate the ASVs
}

\examples{
demultRun(params = "my_run_params.csv", t2s = "t2s_my_run.csv", cores = 12)
dada2workflow(params = "my_run_params.csv", t2s = "t2s_my_run.csv", cores_per_job = 4)
}
