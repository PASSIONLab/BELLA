# BELLA - Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper

BELLA is a computationally efficient and highly accurate long-read to long-read aligner and overlapper. BELLA uses a k-mer seed-based approach to detect overlaps between noisy, long-read data. BELLA provides a novel algorithm for pruning k-mers that are unlikely to be useful in overlap detection and whose presence would only incur unnecessary computational costs. This reliable k-mers detection algorithm explicitly maximizes the probability of retaining k-mers that belong to unique regions of the genome.
To achieve fast overlapping without sketching, BELLA uses sparse matrix-matrix multiplication and utilizes high-performance software and libraries developed for this sparse matrix subroutine. BELLA’s overlap detection has been coupled with a state-of-the-art [seed-and-extend banded-alignment](https://github.com/seqan/seqan) method. BELLA’s alignment step implements a new method to separate true alignments from false positives depending on the alignment score.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

**COMPILER:** the software **requires gcc-6 or higher** with OpenMP to be compiled.  
To run the evaluation test python3 and simplesam package are required. It can be installed via pip: 
```
pip install simplesam
```

### Installing

Clone the repository and enter it:

```
cd bella
```
Build using makefile:

```
ln -s makefile-nersc Makefile && make bella
```

## Running BELLA

To run with default setting:
```
./bella -i <text-file-listing-all-input-fastq-files> -o <out-filename> -d <depth>
```
BELLA requires a text file containing the path to the input fastq file(s) as the argument for the -i option.
Example: [input-example.txt](https://github.com/giuliaguidi/bella/files/2620924/input-example.txt)

To show the usage:
```
./bella -h
```

Optional flag description: 
```
-i : list of fastq(s) (required)
-o : output filename (required)
-d : depth (required)
-k : k-mer length [17]
-a : fixed alignment threshold [50]
-x : alignment x-drop factor [7]
-e : error rate [auto estimated from fastq]
-m : total RAM of the system in MB [auto estimated if possible or 8,000 if not]
-z : skip the pairwise alignment [false]
-w : relaxMargin parameter for alignment on edges [300]
-c : alignment score deviation from the mean [0.1]
-n : filter out alignment on edge [false]
-r : kmerRift: bases separating two k-mers used as seeds for a read [1,000]
-K : all (non-overlapping and separated by <kmerRift> bases) k-mers as alignment seeds [false]
-f : k-mer list from Jellyfish (required if #DEFINE JELLYFISH enabled)
-p : output in PAF format [false]
```
**NOTE**: to use [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) k-mer counting is necessary to enable **#DEFINE JELLYFISH.**

The parallelism depends on the available number of threads and on the available RAM [Default: 8000MB]. Use -DLINUX for Linux or -DOSX for macOS at compile time to estimate available RAM from your machine.

## Output Format

BELLA outputs alignments in a format similar to [BLASR's M4 format](https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format). Example output (tab-delimited):

```HTML
[A ID] [B ID] [# shared k-mers] [alignment score] [overlap length] [n=B fwd, c=B rc] [A start] [A end] [A length] [B start] [B end] [B length]
```
The positions are zero-based and are based on the forward strand, whatever which strand the sequence is mapped.
If **-p** option is used, BELLA outputs alignments in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md). Example output (tab-delimited):

```HTML
[A ID] [A length] [A start] [A end] ["+" = B fwd, "-" = B rc] [B ID] [B length] [B start] [B end] [alignment score] [overlap length] [mapping quality]
```

## Performance Evaluation

The repository contains also the code to get the recall/precision of BELLA and other long-read aligners (Minimap, Minimap2, DALIGNER, MHAP and BLASR).

**Ground truth generation for real data set**: SAMparser.py allows to transform the BWA-MEM/Minimap2 .sam output file in a simpler format usable as input to the evaluation code when using real data set. 

```
python3 SAMparser.py <bwamem/minimap2-output>
```

**Ground truth generation for synthetic data set**: mafconvert.py allows to transform the .MAF file from PBSIM (Pacbio read simulator) in a simpler format usable as input to the evaluation code when using synthetic data set.

```
python mafconvert.py axt <maf-file> > <ground-truth.txt>
```

To run the evaluation program:
```
cd bench
```
```
make result
```
```
./result -G <grouth-truth-file> [-B <bella-output>] [-m <minimap/minimap2-output>] [-D <daligner-output>] [-L <blasr-output>] [-H <mhap-output>] [-M <mecat-output>] [-i <mecat-idx2read-file>]
```
If the output of BELLA is in PAF format, you should run it using minimap2 **-m** flag.

To show the usage:
```
./result -h
```
**NOTE**: add -z flag if simulated data is used.

## Demo 

You can download an _E. coli_ 30X dataset [here](https://bit.ly/2EEq3JM) to test BELLA. For this dataset, you can use the following single mapped ground truth to run the evaluation code: [ecsample_singlemapped_q10.txt](https://github.com/giuliaguidi/bella/files/3143607/ecsample_singlemapped_q10.txt). A detailed description of the procedure we use to generate the ground truth for real data can be found in our [preprint](https://doi.org/10.1101/464420).

You can run the evaluation code located in /bench folder as: 

```./result -G ecsample_singlemapped_q10.txt -B <bella-output>```

## Citation

To cite our work or to know more about our methods, please refer to:

> BELLA: Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper. Giulia Guidi, Marquita Ellis, Daniel Rokhsar, Katherine Yelick, Aydın Buluç. bioRxiv 464420; doi: https://doi.org/10.1101/464420.

## Authors

* [**Giulia Guidi**](https://sites.google.com/berkeley.edu/gguidi/)
* [**Aydın Buluç**](https://people.eecs.berkeley.edu/~aydin/)
* [**Marquita Ellis**](https://sites.google.com/view/about-mme)

## Contributors

* [**Daniel Rokhsar**](https://mcb.berkeley.edu/labs/rokhsar/)
* [**Katherine Yelick**](https://people.eecs.berkeley.edu/~yelick/)

## Copyright Notice
 
Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper (BELLA), Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) Giulia Guidi and Marco Santambrogio. All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
 
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 

## Acknowledgments

Funding provided in part by DOE ASCR through the [Exascale Computing Project](https://www.exascaleproject.org/), and computing provided by [NERSC](https://www.nersc.gov/). Thanks to Rob Egan and Steven Hofmeyr for valuable discussions. Thanks to [Politecnico di Milano](https://www.polimi.it/en/) for key collaborations.
