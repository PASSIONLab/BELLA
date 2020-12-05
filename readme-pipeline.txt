# Compile BELLA and its evaluation program:

```
make bella
cd bench
make result
cd ..
```

# Run pipeline script

Pipeline script is dataset-specific right now, to run with a different dataset, make a copy and modify the path in the script accordingly.

```
./run-bella-pipeline.sh <kmer-len> <window-size> <xdrop> <lower-bound> <upper-bound>
```

```window-size``` has to be specified and = 0 when using regular k-mer.
