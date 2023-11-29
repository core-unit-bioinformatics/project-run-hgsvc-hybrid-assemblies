# Notes on samples

## Verkko (env) updates

Last commit before updating MBG and GraphAligner to get bug fixes for last set of samples:

```
commit: #dbd3f7c88d9e28b052164c650e6ed56b7ba837de
  - graphaligner=1.0.17
  - mbg=1.0.15
```

Updated to

```
commit: #2fe632d38ea615fb93d1af63c462078490e285f9
  - graphaligner=1.0.18
  - mbg=1.0.16
```

for samples:

- YRI trio: NA19238, NA19239, NA19240
- CHS trio: HG00512, HG00513, HG00514
- HG00096
- HG00732 / PUR mother


## Removed

NA19320 - cell line does not grow, insufficient ONT

## Confirmed contamination

NA18939 HiFi - resequencing

## Potential contamination

HG04036 HiFi - v1.4.1+dirty assembly completed
NA21487 HiFi - v1.4.1+dirty assembly completed

# Notes on Verkko

Production version currently is v1.4+dirty [added commits #3119b39 and #4f6a54e]

# Notes on nomenclature

## Trio kmer DBs

HXT - Illumina HiSeq X Ten
NVS - Illumina NovaSeq (6000)
