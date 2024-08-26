# fastqUtils

Set of functions to demultiplex illumina amplicon libraries and running the dada2 workflow. It can account for mixed orientations of reads in R1 and R2 files. Contains convenient wrappers to process a list of libraries in parallel.

## Installation

The fastqUtils package can be installed via:
```
install.packages("remotes")
remotes::install_github("trtcrd/fastqUtils")
```


## Content
The current version contains the following functions:

+ **demultRun()**
Process a list of illumina libraries in parallel (double tagged or dual-indexed). This function is a wrapper of the demultiplexDoublePaired() function (if not accounting fior read orientation) or of the demultiplexDoublePaired4files() function (if accounting for reads mixed orientation). It also run the taggedPrimersDistrib() function on the libraries.
+  **dada2workflow()**
Runs the dada2 default workflow on the demultiplexed libraries (where mixed orientation was not accounted for)
+ **dada2workflow4files()**
Runs the dada2 default workflow on the demultiplexed libraries (where mixed orientation was accounted for)
+ **demultiplexDoublePaired()**
Demultiplex an illumina libraries (not accounting for possible reads mixed orientation)
+ **demultiplexDoublePaired4files()**
Demultiplex an illumina libraries (accounting for possible reads mixed orientation)
+ **taggedPrimersDistrib()**
Sample a multiplexed library and count the reads obtained for each possible combination of tagged primer

## Input file format

Three files are mandatory:
The "params" file, that contains all the instructions to demultiplex the libraries.

library_name|primers|R1|R2|
--- | --- | --- | --- |
lib1|primers_V9.fasta|lib1_R1_S1_sub.fastq.gz|lib1_R2_S1_sub.fastq.gz|
lib2|primers_V9.fasta|lib2_R1_S2_sub.fastq.gz|lib2_R2_S2_sub.fastq.gz|
lib3|primers_V9.fasta|lib3_R1_S3_sub.fastq.gz|lib3_R2_S3_sub.fastq.gz|
...|...|...|...|...

The primers file, a fasta file that contains the tagged primers names and sequences.

```
>forwardPrimer-A
ACCTGCCTAGCGTYG
>forwardPrimer-B
GAATGCCTAGCGTYG
>reversePrimer-B
GAATCTYCAAATCGG
>reversePrimer-C
ACTACTYCAAATCGG
```

The tag-to-sample file, a .csv file

```
run,sample,forward,reverse
library_1,sample_1,forwardPrimer-A,reversePrimer-B
library_1,sample_2,forwardPrimer-B,reversePrimer-C
library_2,sample_3,forwardPrimer-A,reversePrimer-B
library_2,sample_4,forwardPrimer-B,reversePrimer-C
```

## Usage example

There is two wrappers that can be called sequentially:
**demultRun()** demultiplex a list of illumina amplicon libraries (dual-indexed) in parallel.
**dada2workflow()** runs the default dada2 workflow on a list of libraries in parallel.

```
library(fastqUtils)
demultRun(params = "run_params.csv", t2s = "t2s_run_test.csv", accountMixedOrientation = T, cores = 4)
dada2workflow4files(params = "run_params.csv", t2s = "t2s_run_test.csv", cores_per_job = 2)
```

### version 0.1.0 ###

First version
