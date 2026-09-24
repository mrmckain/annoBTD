ARAGORN 1.2.41 (Laslett & Canback 2004, Nucleic Acids Res 32:11-16), GPL v3.
Source from the Debian archive (aragorn_1.2.41.orig.tar.gz). Build:
    cc -O2 -w -o ../../bin/aragorn aragorn1.2.41.c -lm
annoBTD uses it (bin/aragorn) to place tRNA gene ends: post_filter_annotation.pl
snaps each tRNA call to the overlapping ARAGORN gene, then applies the per-identity
window convention in trna_window_convention.tsv. Note: ARAGORN misnames the
intron-containing plastid tRNAs (trnI-GAU as "Thr (cgt)", trnG-UCC as "Ser (cga)");
only its gene ends and intron status are used, never its identity.
