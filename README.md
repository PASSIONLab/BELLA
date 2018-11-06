# BELLA - Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper

BELLA is a computationally-efficient and highly-accurate long-read to long-read aligner and overlapper. BELLA implements a k-mer seed based approach for finding overlaps between pairs of reads. The feasibility of this approach has been demonstrated through a mathematical model based on Markov chains. To achieve fast overlapping without sketching, BELLA exploits sparse matrix-matrix multiplication and utilizes high-performance software and libraries developed for this sparse matrix subroutine.
BELLA applies a simple yet novel procedure for pruning k-mers. We demonstrated that this reliable k-mer selection procedure retains nearly all valuable information with high probability. Our overlap detection has been coupled with state-of-the-art [seed-and-extend banded-alignment methods](https://github.com/seqan/seqan). BELLA implements two different modes for the alignment step. The sensitive mode is fast and achieves high recall, while the precise mode is slightly slower but achieves significantly higher precision. The precise mode uses a new method to separate true alignments from false positives depending on the alignment score.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The software requires gcc-6 or higher with OpenMP to be compiled.  
To run the evaluation test python3 and simplesam package are required. It can be installed via pip: 
```
pip install simplesam
```

### Installing

Clone the repository and enter it:

```
cd longreads
```
Build using makefile:

```
make bella
```

## Running BELLA

To run with default setting:
```
./bella -i <listoffastq> -o <out-filename> -d <depth>
```

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
-a : alignment score threshold [50]"
-p : alignment x-drop factor [3]
-e : error rate [0.15]
-z : skip the pairwise alignment [false]
-w : relaxMargin parameter for alignment on edges [300]
-c : alignment score deviation from the mean [0.1]
-v : PRECISE MODE, use adaptive alignment threshold [false]
-x : ALIGNMENT CONSTRAINT, filter out alignment on edge [false]
-f : k-mer list from Jellyfish (required if #DEFINE JELLYFISH enabled)
```
**NOTE**: to use [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) k-mer counting is necessary to enable **#DEFINE JELLYFISH.**

The multi-threading can be set either depending (a) on the maximum number of thread or (b) on the available RAM. Option (b) should be preferred when medium to large genomes are used. It requires to enable **#DEFINE RAM** in mtspgemm2017/overlapping.h as well as the kind of Operating System used, macOS or Linux.

## Output Format

BELLA outputs alignments in a format similar to [BLASR's M4 format](https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format). Example output (tab-delimited):

```HTML
[A ID] [B ID] [# shared k-mers] [alignment score] [n=B fwd, c=B rc] [A start] [A end] [A length] [B start] [B end] [B length]
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
cd analysis
```
```
make check
```
```
./evaluation -g <grouth-truth-file> [-b <bella-output>] [-m <minimap/minimap2-output>] [-d <daligner-output>] [-l <blasr-output>] [-p <mhap-output>]
```

To show the usage:
```
./evaluation -h
```
**NOTE**: add -z flag if simulated data is used.

To know about the evaluation procedure design and the ground truth generation please refer to:

> Preprint available on biorxiv.

## Authors

* [**Giulia Guidi**](https://linkedin.com/in/giulia-guidi/)
* [**Aydın Buluç**](https://people.eecs.berkeley.edu/~aydin/)
* [**Marquita Ellis**](http://crd.lbl.gov/departments/computer-science/CLaSS/class-staff/marquita-ellis/)

Aydın Buluç and Giulia Guidi are part of the [PASSION Lab](https://passion.lbl.gov/), a research team at the [Computational Research Division](http://crd.lbl.gov) of Lawrence Berkeley National Laboratory and at the [Department of Electrical Engineering and Computer Sciences](https://eecs.berkeley.edu/) of UC Berkeley. 

## Contributors

* [**Daniel Rokhsar**](https://mcb.berkeley.edu/labs/rokhsar/)
* [**Katherine Yelick**](https://people.eecs.berkeley.edu/~yelick/?_ga=2.137275831.646808918.1523950603-1375276454.1515506755)

## Copyright Notice
 
Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper (BELLA), Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) Giulia Guidi and Marco Santambrogio. All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
 
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 

## Acknowledgments

Funding provided in part by DOE ASCR through the Exascale Computing Project, and computing provided by NERSC. Thanks to Rob Egan and Steven Hofmeyr for discussion.
