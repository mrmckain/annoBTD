#!/bin/bash
perl annoBTD/annotate_plastome.pl $1 $4
perl annoBTD/get_annotated_regions_fromverdant.pl $2 $3 $4
./annoBTD/run_multiblast.sh $4_annotated_regions_fromverdant_genes.fsa $4_annotated_regions_fromverdant_trnas.fsa $4_annotated_regions_fromverdant_rrnas.fsa $4_test_seqs_orffinder.fsa $1
perl Scripts/Plastome_Annotation/match_orfs_to_blast_NEW.pl $4_annotated_regions_fromverdant_genes.fsa $4_test_seqs_orffinder.fsa $4_test_seqs_orffinder.fsa_genes.tblastx $4_test_orffinder.txt $1 $4_annotated_regions_fromverdant_trnas.fsa $4_annotated_regions_fromverdant_rrnas.fsa $4_test_seqs_orffinder.fsa_trnas.blastn $4_test_seqs_orffinder.fsa_rrnas.blastn $4


