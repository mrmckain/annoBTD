#!/bin/bash
perl annoBTD_multiref/annotate_plastome.pl $1 $3
COUNT=0;
for file in $(echo $2 | tr "," " "); 
do
    COUNT=$(($COUNT + 1))
    perl annoBTD_multiref/get_annotated_regions_fromverdant.pl ./sequenceFiles/$file ./files/$file $3 $COUNT
done
./annoBTD_multiref/run_multiblast.sh $3_annotated_regions_fromverdant_genes.fsa $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_orffinder_seqs.fsa $1
perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $3_orffinder_seqs.fsa_trnas.blastn $3_annotated_regions_fromverdant_trnas.fsa $1 tRNA $3 $3_orffinder_coordinates.txt

perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $3_orffinder_seqs.fsa_rrnas.blastn $3_annotated_regions_fromverdant_rrnas.fsa Schizachyrium.fsa rRNA $3 $3_orffinder_coordinates.txt

perl annoBTD_multiref/identify_best_ref_for_orf_withscoring_v0.7.pl $3_orffinder_seqs.fsa_genes.tblastx $3_annotated_regions_fromverdant_genes.fsa $3_orffinder_seqs.fsa protein $3 $3_orffinder_coordinates.txt

perl annoBTD_multiref/match_orfs_to_blast_v2.2.pl $3_annotated_regions_fromverdant_genes.fsa $3_orffinder_seqs.fsa $3_best_orfs_for_refs_SCORE.txt $3_orffinder_coordinates.txt $1 $3_annotated_regions_fromverdant_trnas.fsa $3_annotated_regions_fromverdant_rrnas.fsa $3_orffinder_seqs.fsa_trnas.blastn $3_orffinder_seqs.fsa_rrnas.blastn $3 $3_best_tRNA_ref_SCORE.txt $3_best_rRNA_ref_SCORE.txt

