# BELLA - Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper

BELLA is a computationally efficient and highly accurate long-read to long-read aligner and overlapper. BELLA uses a k-mer seed-based approach to detect overlaps between noisy, long-read data. BELLA provides a novel algorithm for pruning k-mers that are unlikely to be useful in overlap detection and whose presence would only incur unnecessary computational costs. This reliable k-mers detection algorithm explicitly maximizes the probability of retaining k-mers that belong to unique regions of the genome.
To achieve fast overlapping without sketching, BELLA uses sparse matrix-matrix multiplication and utilizes high-performance software and libraries developed for this sparse matrix subroutine. BELLA’s overlap detection has been coupled with a state-of-the-art [seed-and-extend banded-alignment](https://github.com/seqan/seqan) method. BELLA’s alignment step implements a new method to separate true alignments from false positives depending on the alignment score.

## Content

*	[Getting Started on Linux](#getting-started-on-linux)
	*	[Dependencies](#dependencies)
	*	[Compile](#compile)
	*	[Run](#run)
	*	[Error Rate](#error-rate)
	*	[K-mer Counting](#k-mer-counting)
	*	[Memory Usage](#memory-usage)
*	[Output Format](#output-format)
*	[Performance Evaluation](#performance-evaluation)
*	[Demo](#demo)
*	[I get 0 outputs, what is likely going wrong?](#i-get-0-outputs-what-is-likely-going-wrong)
*	[Citation](#citation)

## Getting Started on Linux

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. This version of BELLA **only works on Linux-based machines**. To run BELLA on macOS-based machines, please switch to **mac branch**.

### Dependencies

* **COMPILER:** the software **requires gcc-6** with OpenMP to be compiled.
* **BOOST/1.67.0** to use Gerbil kmerCounting.
You can install BOOST/1.67.0 using [conda](https://anaconda.org/anaconda/boost):
```
conda install -c anaconda boost
```
* [**CUDA**](https://docs.nvidia.com/cuda/) to compile and use GPU-accelerated Gerbil. You **do not** need CUDA to use CPU-based Gerbil.

* **Python3** and **simplesam** are required to generare the ground truth data. You can install simplesam via [pip](https://pip.pypa.io/en/stable/installing/): 
```
pip install simplesam
```

### Compile

Clone the repository, its submodule, and enter it:

```
git clone https://github.com/giuliaguidi/bella
cd bella
```
OR:
```
git clone https://github.com/giuliaguidi/bella
cd bella
git submodule init
git submodule update
```
Build using makefile:

```
ln -s makefile-mac Makefile && make bella
```

### Run

To run with default setting:
```
./bella -i <text-file-listing-all-input-fastq-files> -o <out-filename> -d <coverage>
```
BELLA requires a text file containing the path to the input fastq file(s) as the argument for the -i option.
Example: [input-example.txt](https://github.com/giuliaguidi/bella/files/2620924/input-example.txt)

To show the usage:
```
./bella -h
```

Optional flag description: 
```
-f : List from Jellyfish (required if Jellyfish kmerCounting is used)
-i : List of fastq(s)	(required)
-o : Output filename	(required)
-d : Dataset coverage	(required)
-k : KmerSize [17]
-a : User-defined alignment threshold [FALSE, 0]
-x : SeqAn xDrop [7]
-e : Error rate [0.15]
-q : Estimare error rate from the dataset [FALSE]
-u : Use default error rate setting [FALSE]
-g : Use Gerbil as kmerCounter [FALSE]
-y : Enable GPU [FALSE]
-m : Total RAM of the system in MB [auto estimated if possible or 8,000 if not]
-z : Do not run pairwise alignment [FALSE]
-c : Deviation from the mean alignment score [0.10]
-r : KmerRift: bases separating two k-mers [kmerSize]
-s : Common k-mers threshold to compute alignment [auto estimated if possible]
-b : Bin size binning algorithm [500]
-p : Output in PAF format [FALSE]
-w : Probability threshold for reliable range [0.002]
```
### Error Rate

The error rate is an important parameter in BELLA as it is used to choose which k-mers contribute to the overlap detection.

The user should either:

* **-e** = suggest an error rate
* **-q** = confirm that the data has quality values and we can estimate the error rate from the data set
* **-u** = confirm that we can use a default error rate (0.15)

### K-mer Counting

BELLA can run with three different k-mer counting options:

* **Default**: BELLA uses its own fast k-mer counter based on a [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter) data structure. This is the fastest CPU-based option but it is limited by the available RAM. If BELLA goes **out-of-memory during the k-mer counting stage**, you should use BELLA version on the **master branch** and use Gerbil k-mer counter. **Gerbil k-mer counter only works on Linux-based machines**.
* **Jellyfish**: BELLA uses [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) k-mer counter. It is necessary to install Jellyfish, add **-DJELLYFISH** when compiling BELLA, and give Jellyfish output file to BELLA as input parameter. Differently from Gerbil, the k-mer counting does not happen within BELLA.

### Memory Usage

The parallelism during the overlap detection phase depends on the available number of threads and on the available RAM [Default: 8000MB].

Use **-DOSX** at compile time to estimate available RAM from your machine. 
If your machine has more RAM than the default one, using **-DOSX** would **make the ovelap detection phase faster**. 

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

Error rate estimation might have gone wrong. If the error estimated is greater than 1, the adaptive alignment threshold would be so high that no alignments would pass the threshold. Please check if your fastq file has proper quality values. If not, please define an error rate using command line options.

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
* [**Qi Zhou**](https://it.linkedin.com/in/qizhou1512)

## Copyright Notice
 
Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper (BELLA), Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) Giulia Guidi and Marco Santambrogio. All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
 
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 

## Acknowledgments

Funding provided in part by DOE ASCR through the [Exascale Computing Project](https://www.exascaleproject.org/), and computing provided by [NERSC](https://www.nersc.gov/). Thanks to Rob Egan and [Steven Hofmeyr](https://crd.lbl.gov/departments/computer-science/CLaSS/members/class-staff/steven-hofmeyr/) for valuable discussions. Thanks to [Politecnico di Milano](https://www.polimi.it/en/) for key collaborations.
