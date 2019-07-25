# Gerbil: A fast and memory-efficient k-mer counter with GPU-support

A basic task in bioinformatics is the counting of k-mers in
genome strings. The k-mer counting problem is to build a histogram of
all substrings of length k in a given genome sequence. 
Gerbil is a k-mer counter that is specialized for high effiency when counting k-mers for large k. 
The software is decribed [here](https://almob.biomedcentral.com/articles/10.1186/s13015-017-0097-9).

To cite Gerbil in publication, please use
> Marius Erbert, Steffen Rechner, and Matthias Müller-Hannemann,
> Gerbil: A fast and memory-efficient k-mer counter with GPU-support,
> Algorithms for Molecular Biology (2017) 12:9, open access.

## Changelog

### Version 1.11
  * Minor Bugfixes and enhanced tolerance when reading malformed fasta files

### Version 1.1
  * Improved performance while reading compressed input files
  * Improved GPU performance
  * Added option for creating `fasta` output
  * Minor bugfixes

### Version 1.0
  * Initial upload
  

## Install

Gerbil is developed and tested at Linux operating systems. Migrating it to other OS like Windows is a current issue. It follows a description of the installation process at Ubuntu 16.04.

1. Install 3rd-party libraries and neccessary software:

        sudo apt-get install git cmake g++ libboost-all-dev libz3-dev libbz2-dev

2. Download the Source Files. 

        git clone https://github.com/uni-halle/gerbil.git
        
3. Compile the Sources. Gerbil comes with a CMake Script that shold work for various operating systems. CMake will automatically detect whether all mandatory and optional libraries are available at your system.

        cd gerbil
        mkdir build
        cd build
        cmake ..
        make

The `build` directory should now contain a binary `gerbil`.

## Usage

        gerbil [option|flag]* <input-file> <temp-directory> <output-file>

Gerbil can be controlled by several command line options and flags.

| Option               | Description   | Default |
|:---------------------|:--------------| -------:|
| `‑k <int>`   | Set the length of k-mers. Supported k currently range from 8 to 136. Support for values larger than 136 can easily be activated if needed. | 28 |
| `‑m <int>`          | Set the length m of minimizers.      |   auto |
| `‑e <int>MB`  | Restrict the maximal size of main memory Gerbil is allowed to use to x MB.      |    auto |
| `‑e <int>GB`  | Restrict the maximal size of main memory Gerbil is allowed to use to x GB.      |    auto |
| `‑o <opt>`    | Change the format of the output. Valid options for `<opt>` are `gerbil`, `fasta` and `none`.      |    `gerbil` |
| `‑t <int>`          | Set the maximal number of parallel threads to use.      |    auto |
| `‑l <int>`               | Set the minimal occurrence of a k-mer to be outputted.      |    3 |
| `‑i`                   | Enable additional debug output.      |    |
| `‑g`                   | Enable GPU mode. Gerbil will automatically detect CUDA-capable devices and will use them for counting in the second phase.      |     |
| `‑v`                   | Show version number.      |     |
| `‑d`                   | Disable normalization of k-mers. If normalization is disabled, a k-mer and its reverse complement are considered as different k-mers. If normalization is enabled, we map both k-mer and its reverse complement to the same k-mer.       |     |
| `‑s`                   | Perform a system check and display information about your system.     |     |
| `‑x 1`                 | Stop execution after Phase One. Do not remove temporary files and `binStatFile` (with statistical information). When using this option, no `output` parameter is allowed. |     |
| `‑x 2`            | Execute only Phase Two. Requires temporary files and `binStatFile`. No `input` parameter is allowed. |     |
| `‑x b`            | Do not remove `binStatFile`. |     |
| `‑x h`            | Create a histogram of k-mers in a human readable format in output directory. |     |

## Input Formats

Gerbil supports the following input formats of genome read data in raw and compressed format: 
 * `fastq`, `fastq.gz`, `fastq.bz2`
 * `fasta`, `fasta.gz`, `fastq.bz2`
 * `staden`
 * `txt`: A plain text file with one path per line. This way, multiple input files can be processed at once.

## Output Format

Gerbil uses an output format that is easy to parse and requires little space. The counter of each occuring k-mer is stored in binary form, followed by the corresponding byte-encoded k-mer. Each four bases of a k-mer are encoded in one single byte. We encode `A` with `00`, `C` with `01`, `G` with `10` and `T` with `11`. Most counters of k-meres are slightly smaller than the coverage of the genome data. We exploit this property by using only one byte for counters less than 255. A counter greater than or equal to 255 is encoded in five bytes. In the latter case, all bits of the first byte are set to 1. The remaining four bytes contain the counter in a conventional 32-bit unsigned integer.

Examples (`X` means undefined):

| Counter | k-mer   | Encoding                      |
|:--------|:--------|:------------------------------|
| 67      | AACGTG  | `01000011 00000110 1110XXXX` |
| 345     | TGGATC  | `11111111 00000000 00000000 00000001 01011001 11101000 1101XXXX` |

The output file can be converted into `fasta` format by running the command

        toFasta <gerbil-output> <k> [<fasta-output>]

