## Getting Started

These instructions will guide you throught the utilization of the scripts we developed related to the main project, [BELLA](https://bitbucket.org/aydozz/longreads).

### Prerequisites

```
Text it
```

## Overlapping Feasibility via Markov Chain Model

In the paper, we present a theory on how k-mers (subsequences) can be used to anchor alignments between two long reads, allowing an accurate overlap discovery among all the reads in a dataset.  
To run the script implementing the Markov Chain Model we developed to show the overlapping feasibility, do:
```
python exactMarkov.py
```
The default function setting is:
```python
func(overlap=500, p=0.85, k=17)
```

The program with default setting computes the probability to get 17 consecutive successes on both the sequences given an overlap length of 500bp and a sequencing success probability equal to 0.85.

## Scripts for Consolidating dAligner alignment results in BELLA format

```
Run translateDalignerOut.sh with a script specifiying the desired local alignments (.las file) and corresponding (DAZZ_DB) database instance,
along with the necessary executables. (The full detailed requirements list is specified in translateDalignerOut.sh).
An example setup script (for NERSC's Cori) is example-environment-script.sh.   
```

## Further script?

```
Text here
```

## Authors

* **Giulia Guidi**
* [**Aydın Buluç**](https://people.eecs.berkeley.edu/~aydin/)
* **Marquita Ellis**

## Contributors

* **Daniel Rokhsar**
* [**Katherine Yelick**](https://people.eecs.berkeley.edu/~yelick/?_ga=2.137275831.646808918.1523950603-1375276454.1515506755)

## Copyright Notice
 
Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper (BELLA), Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) Giulia Guidi and Marco Santambrogio. All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
 
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
