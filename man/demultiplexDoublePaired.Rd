\encoding{UTF-8}
\name{demultiplexDoublePaired}
\alias{demultiplexDoublePaired}
\title{Demultiplexing illumina amplicon libraries}

\description{
Demultiplexing an illumina amplicon library (double tags, or dual indexing) without accounting for possible mixed orientation of reads in R1 and R2.
}

\usage{
  demultiplexDoublePaired(primersFile, t2s, fastqR1In, fastqR2In,
    allowedMis = 0, outputFolder = "demult", overwrite = F,
    splitHeader = " ", chunkSize = 1000000, progressBar = T)
}

\arguments{
  \item{primersFile}{Path to a fasta file containing tagged primers names and sequences).}
  \item{t2s}{A "tag-to-sample" .csv file containing all informations for demultiplexing samples of a single library. Must contains columns names "run" (for library name, to.do: harmonize this field across functions), "sample" (sample name, cannot be duplicated in a given library), "forward" and "reverse" (for the tagged primers combination).}
  \item{fastqR1In}{Path to the R1 fastq file.}
  \item{fastqR2In}{Path to the R2 fastq file.}
  \item{allowedMis}{Number of mismatches allowed in searching primers sequences in the reads (default = 0). If you set it too high (i.e. 2 or more), you may have reads ascribed to multiple samples. In that case, a warning will be thrown with a file with more details.}
  \item{outputFolder}{Output folder. Can be specified. Default is "demult".}
  \item{overwrite}{If outputFolder already exists, will throw an error. Shall we overwrite the folder (default = FALSE)?}
  \item{splitHeader}{Character for parsing the fastq header names and checking consistency between R1 and R2 files (default = " "). Most likely not to be touched, unless a specific error about that is thrown.}
  \item{chunkSize}{How many reads to be processed per batch (default = 1000000). The higher the number, the more RAM is needed. Default involves about 1 Go per library.}
  \item{progressBar}{Showing a progress bar. Default = T.}
  }

\details{
Demultiplex a single amplicon library without accounting for possible mixed reads orientation in the R1 and R2 files. Two fastqs files will be produced for each sample: fwd and rev.
}

\examples{
demultiplexDoublePaired(primersFile = "my_primers.fasta", t2s = "t2s_my_library.csv", fastqR1In = "Lib1_R1.fastq.gz", fastqR2In = "Lib1_R2.fastq.gz")
}
