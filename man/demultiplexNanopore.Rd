\encoding{UTF-8}
\name{demultiplexNanopore}
\alias{demultiplexNanopore}
\title{Demultiplexing a nanopore amplicon library}

\description{
Demultiplexing a nanopore library that is double tagged, i.e. forward and reverse primers. Reads are re-oriented using the forward primer.
}

\usage{
  demultiplexNanopore(primersFile, t2s, fastqIn, allowedMis = 0, trim = T, withIndels = F, outputFolder = "demult", overwrite = F, chunkSize = 100000, progressBar = T)
}

\arguments{
  \item{primersFile}{Path to a fasta file containing tagged primers names and sequences).}
  \item{t2s}{A "tag-to-sample" .csv file containing all informations for demultiplexing samples of a single library. Must contains columns names "run" (for library name, to.do: harmonize this field across functions), "sample" (sample name, cannot be duplicated in a given library), "forward" and "reverse" (for the tagged primers combination).}
  \item{fastqIn}{Path to the R1 fastq file.}
  \item{allowedMis}{Number of mismatches allowed in searching primers sequences in the reads (default = 0). If you set it too high (i.e. 2 or more), you may have reads ascribed to multiple samples. In that case, a warning will be thrown with a file with more details.}
  \item{trim}{Trimming the primers from the reads? Default = T}
  \item{withIndels}{Allow indels when searching for the primers? Default = F}
  \item{outputFolder}{Output folder. Can be specified. Default is "demult"}
  \item{overwrite}{If outputFolder already exists, will throw an error. Shall we overwrite the folder (default = FALSE)?}
  \item{chunkSize}{How many reads to be processed per batch (default = 1000000). The higher the number, the more RAM is needed. Default involves about 1 Go per library.}
  \item{progressBar}{Showing a progress bar. Default = T}
  }

\details{
Demultiplex a single nanopore amplicon library and re-orient reads of mixed orientation.
}

\examples{
demultiplexNanopore(primersFile = "my_primers.fasta", t2s = "t2s_my_library.csv", fastqIn = "Lib_nano.fastq.gz")
}
