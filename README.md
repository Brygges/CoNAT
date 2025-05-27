# CoNAT
A Toolkit for the Identification and Clustering of Novel Animal Toxins.

> [!IMPORTANT]
> Files necessary for CoNAT to run are not accessible in this Github Repo. This repository serves to highlight the code itself.

### Identification
The identification directory contains the three different scripts employed for identification of novel animal toxins.
For ease of use, a shell-script has been supplied (`run.sh`) that accepts an input abundance file and will yield an output directory with results from analysis.
```
$ bash run.sh <input> <output-dir> 
```

### Steps in-between
Structural predictions are not currently part of the pipeline - as such, users must first go to external ressources for this purpose.
For CoNAT to function, a tab-seperated file must be supplied containing query, targets and e-values. These files are made using Foldseek, by making an all-against-all search against the database itself. We refer to Foldseek's own documentation on their [GitHub Page](https://github.com/steineggerlab/foldseek).
```
$ foldseek easy-search <database> <database> result.tsv tmp --format-output "query,target,evalue" -e 1.0E-3
```

### Clustering
CoNAT can be run using either default settings or bare-bones mode (`--bare`), depending on whether superpositions should be carried out.
```
$ python CoNAT.py
$ python CoNAT.py --bare
```
The Python script will yield a `index.html` file, which can be opened using the `flaskrun.py` script.
