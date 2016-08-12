#!/bin/bash

/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/makeblastdb -in $1 -dbtype nucl
/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/makeblastdb -in $2 -dbtype nucl
/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/makeblastdb -in $3 -dbtype nucl

/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/blastn -query $5 -db $2 -outfmt 6 -evalue 1e-1 > $4_trnas.blastn
/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/blastn -query $5 -db $3 -outfmt 6 -evalue 1e-1 > $4_rrnas.blastn
/var/www/plastidDB/annoBTD_multiref/ncbi-blast-2.2.30+/bin/tblastx -query $4 -db $1 -outfmt 6 -evalue 1e-1 > $4_genes.tblastx

