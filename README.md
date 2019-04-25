# seq2Fun
__seq2Fun__ provides function to annotate protein or coding sequences with KEGG pathways
## Description
Annotate the protein and coding sequences with KEGG pathways could be done with [_KAAS_](https://www.genome.jp/kegg/kaas/) website. Here the _seq2Fun_ try to do same thing with _blast_ function.
## Requirement
The blast+ should be installed before installing _seq2Fun_ package
## Installation
```
library(devtools)
install_github("guokai8/seq2Fun")
``` 

## Getting started

```
library(seq2Fun)
## check if blast tools have been installed
blast_help()
db <- preparedb(species = "Arabidopsis thaliana", seqtype = "AA", savedb = TRUE) 
## Take several minutes
str(db, 1)
###savedb will write out the sequences file in the work directory
makeblastdb(db$db, dbtype = "prot")
###make blast db  
seqs <- db$db[sample(500, 10)] ## random choose 10 sequences
ann <- seq2Fun(query = seqs, dbfile = db$db, anno = db$anno, evalue = 1e-10, num_threads = 2)
head(ann)
```

## Contact information

For any questions please contact guokai8@gmail.com
