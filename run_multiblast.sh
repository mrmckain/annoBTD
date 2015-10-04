#!/bin/bash

~/Downloads/ncbi-blast-2.2.30+/bin/makeblastdb -in $1 -dbtype nucl
~/Downloads/ncbi-blast-2.2.30+/bin/makeblastdb -in $2 -dbtype nucl
~/Downloads/ncbi-blast-2.2.30+/bin/makeblastdb -in $3 -dbtype nucl

~/Downloads/ncbi-blast-2.2.30+/bin/blastn -query $5 -db $2 -outfmt 6 -evalue 1e-1 > $4_trnas.blastn
~/Downloads/ncbi-blast-2.2.30+/bin/blastn -query $5 -db $3 -outfmt 6 -evalue 1e-1 > $4_rrnas.blastn
~/Downloads/ncbi-blast-2.2.30+/bin/tblastx -query $4 -db $1 -outfmt 6 -evalue 1e-1 > $4_genes.tblastx

