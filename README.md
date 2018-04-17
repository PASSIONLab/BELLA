# BELLA - Berkeley Efficient Long-read to Long-read Overlapper and Aligner

BELLA is a computationally-efficient and highly-accurate long-read to long-read aligner and overlapper. BELLA implements a k- mer seed based approach for finding overlaps between pairs of reads. The feasibility of this approach has been demonstrated through a mathematical model based on Markov chains. To achieve fast overlapping without sketching, BELLA exploits sparse matrix-matrix multiplication and utilizes high-performance software and libraries developed for this sparse matrix subroutine.
BELLA applies a simple yet novel procedure for pruning k-mers. We demonstrated that this reliable k-mer selection procedure retains nearly all valuable information with high probability. Our overlap detection has been coupled with state-of-the-art seed-and-extend banded-alignment methods. BELLA attains high recall.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The software require g++-6 and OpenMP to be compiled.
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

Describe optional flags.

### Performance evaluation

The repository contains also the code to get the recall/precision of BELLA and other long-read aligners (Minimap, Minimap2, DALIGNER, MHAP and BLASR).

**Ground truth generation for real data set**: SAMparser.py allows to transform the BWA-MEM .sam outfile in a simpler format usable as input to the evaluation code when using real data set. 

```
python3 SAMparser.py <bwamem-output>
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
**NOTE**: add -z flag if synthetic data is used.

### Overlapping feasibility via Markov Chain Model

Explain what these tests test and why

```
Give an example
```

## Built With

* [GNU Make](https://www.gnu.org/software/make/)

## Authors

* **Giulia Guidi**
* **Aydin Buluc** [Personal website](https://people.eecs.berkeley.edu/~aydin/)

## Contributors

* **Marquita Ellis**
* **Daniel Rokhsar**
* **Katherine Yelick** [Personal website](https://people.eecs.berkeley.edu/~yelick/?_ga=2.137275831.646808918.1523950603-1375276454.1515506755)

## License

Add license.

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
