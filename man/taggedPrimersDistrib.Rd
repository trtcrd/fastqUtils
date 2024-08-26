\encoding{UTF-8}
\name{taggedPrimersDistrib}
\alias{taggedPrimersDistrib}
\title{Sample a multiplexed library and count the reads obtained for each possible combination of tagged primer}

\description{
Randomly sample the fastq files of a multiplexed library and count the number of reads containing each possible combination of tagged primers (i.e. not only the expected combination for the list pf samples). The function exports the reads counts and produces a heatmap.
}

\usage{
  taggedPrimersDistrib(primersFile, t2s, fastqR1In, fastqR2In, allowedMis = 0, outputFolder = "demult", subsample = .05, splitHeader = " ", progressBar = T)
}

\arguments{
  \item{primersFile}{Path to a fasta file containing tagged primers names and sequences).}
  \item{t2s}{A "tag-to-sample" .csv file containing all informations for demultiplexing samples of a single library. Must contains columns names "run" (for library name, to.do: harmonize this field across functions), "sample" (sample name, cannot be duplicated in a given library), "forward" and "reverse" (for the tagged primers combination).}
  \item{fastqR1In}{Path to the R1 fastq file.}
  \item{fastqR2In}{Path to the R2 fastq file.}
  \item{allowedMis}{Number of mismatches allowed in searching primers sequences in the reads (default = 0).}
  \item{outputFolder}{Output folder. Can be specified. Default is "demult"}
  \item{subsample}{Proportion of reads to be sampled from the library (default = 0.05)}
  \item{splitHeader}{Character for parsing the fastq header names and checking consistency between R1 and R2 files (default = " "). Most likely not to be touched, unless a specific error about that is thrown.}
  \item{progressBar}{Showing a progress bar. Default = T}
}

\details{
This can be seen as a quality control step. Expected combinations of primers (i.e. our samples) should have a large number of reads. Some pipetting mistakes during PCR plate preparation (e.g. adding two primers combinations in a single well) can be seen on the heatmap.
}

\examples{
taggedPrimersDistrib(primersFile = "my_primers.fasta", t2s = "t2s_my_library.csv", fastqR1In = "Lib1_R1.fastq.gz", fastqR2In = "Lib1_R2.fastq.gz")
}
