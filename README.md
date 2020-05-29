# BELLA - Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper

BELLA is a computationally efficient and highly accurate long-read to long-read aligner and overlapper. BELLA uses a k-mer based approach to detect overlaps between noisy, long reads. We demonstrated the feasibility of the k-mer based approach through a mathematical model based on Markov chains. BELLA provides a novel algorithm for pruning k-mers that are unlikely to be useful in overlap detection and whose presence would only incur unnecessary computational costs. Our reliable k-mers detection algorithm explicitly maximizes the probability of retaining k-mers that belong to unique regions of the genome.
BELLA achieves fast overlapping without sketching using sparse matrix-matrix multiplication (SpGEMM), implemented utilizing high-performance software and libraries developed for this sparse matrix subroutine. Any novel sparse matrix format and multiplication algorithm would be applicable to overlap detection and enable continued performance improvements. We coupled BELLA's overlap detection with [our newly developed vectorized seed-and-extend banded-alignment algorithm](https://github.com/giuliaguidi/xavier).
The choice of the optimal k-mer seed occurs through our binning mechanism, where k-mer positions within a read pair are used to estimate the length of the overlap and to "bin" k-mers to form a consensus.We developed and implemented a new method to separate true alignments from false positives depending on the alignment score. This method demonstrates that the probability of false positives decreases exponentially as the overlap length between sequences increases.

## Content

*	[Getting Started](#getting-started)
	*	[Dependencies](#dependencies)
	*	[Compile](#compile)
	*	[Run](#run)
	*	[Error Rate](#error-rate)
	*	[Memory Usage](#memory-usage)
*	[Output Format](#output-format)
*	[Performance Evaluation](#performance-evaluation)
*	[Demo](#demo)
*	[I get 0 outputs, what is likely going wrong?](#i-get-0-outputs-what-is-likely-going-wrong)
*	[Citation](#citation)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Dependencies

* **COMPILER:** the software **requires gcc-6 or higher** with OpenMP to be compiled.
* [**CUDA**](https://docs.nvidia.com/cuda/) to compile and use GPU-accelerated pairwise alignment. You **do not** need CUDA to use CPU-based pairwise alignment. Our stand-alone GPU-based pairwise alignment, named **LOGAN**, can be found [**here**](https://github.com/albertozeni/LOGAN).

* **Python3** and **simplesam** are required to generare the ground truth data. You can install simplesam via [pip](https://pip.pypa.io/en/stable/installing/): 
```
pip install simplesam
```

### Compile

Clone the repository and enter it:

```
git clone https://github.com/giuliaguidi/bella
cd bella
```
Build using makefile:

```
ln -s makefile-nersc Makefile
make bella (CPU-only) OR make bella-gpu (CPU/GPU)
```

### Run

To run with default setting:
```
./bella -f <text-file-listing-all-input-fastq-files> -o <out-filename> -c <coverage> -q
```
BELLA requires a text file containing the path to the input fastq file(s) as the argument for the -f option.
Example: [input-example.txt](https://github.com/giuliaguidi/bella/files/2620924/input-example.txt)

To show the usage:
```
./bella -h
```

Optional flag description: 
```
-f : List of fastq(s)	(required)
-o : Output filename	(required)
-c : Dataset coverage	(required)
-k : KmerSize [17]			
-a : User-defined alignment threshold [FALSE, -1]   		
-x : SeqAn xDrop [7]   									 
-e : Error rate [0.15]   				 
-q : Estimare error rate from the dataset [FALSE]   	 
-u : Use default error rate setting [FALSE]  		   	 
-m : Total RAM of the system in MB [auto estimated if possible or 8,000 if not]  	 
-z : Do not run pairwise alignment [FALSE]   			 
-d : Deviation from the mean alignment score [0.10]  	 
-w : Bin size binning algorithm [500]   	 
-p : Output in PAF format [FALSE]   		 
-r : Probability threshold for reliable range [0.002]     
-g : GPUs available [1, only works when BELLA is compiled for GPU] 
-s : K-mer counting split count can be increased for large dataset [1]	 
```
### Error Rate

The error rate is an important parameter in BELLA as it is used to choose which k-mers contribute to the overlap detection.

The user should either:

* **-e** = suggest an error rate
* **-q** = confirm that the data has quality values and we can estimate the error rate from the data set
* **-u** = confirm that we can use a default error rate (0.15)

### Memory Usage

The parallelism during the overlap detection phase depends on the available number of threads and on the available RAM [Default: 8000MB].

Use **-DOSX** or **-DLINUX** at compile time to estimate available RAM from your machine. 
If your machine has more RAM than the default one, using **-DOSX** or **-DLINUX** would **make the ovelap detection phase faster**. 

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

* **Ground truth generation for real data set**: SAMparser.py allows to transform the Minimap2 .sam output file in a simpler format usable as input to the evaluation code when using real data set. 

```
minimap2 -ax map-pb  ref.fa pacbio-reads.fq > aln.sam   # for PacBio subreads
samtools view -h -Sq 10 -F 4 aln.sam > mapped_q10.sam	# remove reads with quality values smaller than 10
samtools view -h mapped_q10.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -S -h > unique_mapped_q10.sam	# remove reads mapped to multiple locations
python3 SAMparser.py <bwamem/minimap2-output>
```

* **Ground truth generation for synthetic data set**: mafconvert.py allows to transform the .MAF file from PBSIM (Pacbio read simulator) in a simpler format usable as input to the evaluation code when using synthetic data set.

```
python mafconvert.py axt <maf-file> > <ground-truth.txt>
```

To run the evaluation program:
```
cd bench
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

## I get 0 outputs, what is likely going wrong?

<p align="justify">
Error rate estimation might have gone wrong. If the error estimated is greater than 1, the adaptive alignment threshold would be so high that no alignments would pass the threshold. Please check if your fastq file has proper quality values. If not, please define an error rate using command line options.
</p>

## Citation

To cite our work or to know more about our methods, please refer to:

> BELLA: Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper. Giulia Guidi, Marquita Ellis, Daniel Rokhsar, Katherine Yelick, Aydın Buluç. bioRxiv 464420; doi: https://doi.org/10.1101/464420.

## Authors

* [**Giulia Guidi**](https://sites.google.com/berkeley.edu/gguidi/)
* [**Aydın Buluç**](https://people.eecs.berkeley.edu/~aydin/)
* [**Marquita Ellis**](https://sites.google.com/view/about-mme)

## Contributors

* [**Dan Rokhsar**](https://mcb.berkeley.edu/labs/rokhsar/)
* [**Kathy Yelick**](https://people.eecs.berkeley.edu/~yelick/)
* [**Alberto Zeni**](https://it.linkedin.com/in/alberto-zeni-b61077158)
* [**Elizabeth Koning**](https://www.linkedin.com/in/elizabeth-koning)
* [**Qi Zhou**](https://www.linkedin.com/in/qizhou1512/?originalSubdomain=it)

## Copyright Notice

<p align="justify">
Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper (BELLA), Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) Giulia Guidi and Marco Santambrogio. All rights reserved.
</p>

<p align="justify">
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
</p>
<p align="justify">
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 
</p>

## Acknowledgments

Funding provided in part by DOE ASCR through the [Exascale Computing Project](https://www.exascaleproject.org/), and computing provided by [NERSC](https://www.nersc.gov/). Thanks to Rob Egan and [Steven Hofmeyr](https://crd.lbl.gov/departments/computer-science/CLaSS/members/class-staff/steven-hofmeyr/) for valuable discussions. Thanks to [NECST Laboratory](https://necst.it/) and Ed Younis for key collaborations.

