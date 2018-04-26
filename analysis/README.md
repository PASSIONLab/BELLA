## Getting Started

These instructions will guide you throught the utilization of the scripts we developed related to the main project, [BELLA](https://bitbucket.org/aydozz/longreads).

### Prerequisites

```
Text it
```

## Overlapping feasibility via Markov Chain Model

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

## dAligner Parser

```
Text here
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

## License

Add license.

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
