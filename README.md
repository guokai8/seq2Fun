# seq2Fun
__seq2Fun__ provides function to annotate protein or coding sequences with KEGG pathways
## Description
Annotate the protein and coding sequences with KEGG pathways could be done with [_KAAS_](https://www.genome.jp/kegg/kaas/) website. Here the _seq2Fun_ try to do same thing with _blast+_ and _diamond_ software.
## Requirement
The _blast+_ or _diamond_ should be installed before installing _seq2Fun_ package
## Installation
```
library(devtools)
install_github("guokai8/seq2Fun")
``` 

## Software Usage

```
library(seq2Fun)
## check if blast+ or diamond tools have been installed
blast_help()
head(listspecies)
###Find the correct species name
db <- preparedb(species = "Arabidopsis thaliana", seqtype = "AA", savedb = TRUE) 
## Take several minutes
str(db, 2)
###savedb will write out the sequences file in the work directory
makeblastdb(db, dbtype = "prot")
###make blast db, set runblast = FALSE if you prefer diamond
seqs <- db$db[sample(500, 10)] ## random choose 10 sequences
ann <- seq2fun(query = seqs, db = db, type = "blastp", evalue = 1e-10, num_threads = 2)
## set bidirectional = TRUE if you prefer bidirectional blast
## set runblast = FALSE if you prefer diamond
head(ann)
```
## Note
The seq2Fun downloads and uses KEGG data. Non-academic uses may require a KEGG
license agreement (details at http://www.kegg.jp/kegg/legal.html).

## Contact information

For any questions please contact guokai8@gmail.com
