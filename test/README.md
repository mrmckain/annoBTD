# annoBTD regression harness

Measures annotation accuracy against GenBank ground truth, so changes to the
boundary logic can be scored instead of eyeballed.

## Quick start

```bash
cd annoBTD/test
./fetch_test_set.sh        # download 7 plastomes + build truth tables (once)
./run_regression.sh        # leave-one-out annotate + score all 7
```

Needs BLAST+ on `PATH`, or `BLAST_BIN=/path/to/ncbi-blast/bin`.

Score a single genome, or point at a different output directory:

```bash
./run_regression.sh --only Nicotiana_tabacum --outdir results_myfix
```

## What it does

For each genome in the test set, every *other* genome is used as a guide
reference and the pipeline is run leave-one-out. The result is compared to a
truth table built from that genome's own GenBank record.

The test set spans land plants so a change that helps one lineage and hurts
another is visible:

| genome | note |
|---|---|
| *Nicotiana tabacum* | the classic reference plastome |
| *Arabidopsis thaliana* | rosid |
| *Spinacia oleracea* | caryophyllid |
| *Zea mays*, *Oryza sativa* | monocots |
| *Pinus thunbergii* | gymnosperm, drastically reduced IR, no *ndh* genes |
| *Marchantia polymorpha* | liverwort, most distant guide set |

## Files

| file | role |
|---|---|
| `fetch_test_set.sh` | downloads GenBank + FASTA from NCBI, builds truth tables |
| `gb_to_truth.pl` | GenBank flatfile → truth table |
| `compare_annotations.pl` | boundary-accuracy diff, prediction vs truth |
| `summarize.pl` | pulls headline numbers into `results/summary.tsv` |
| `run_regression.sh` | drives the whole leave-one-out loop |
| `refresh_references.sh` | rebuilds the derived reference tables below from the GenBank cache |
| `unit_boundary.pl` | unit tests for the boundary-finding subs |

### Derived reference tables

Everything in this directory that is not a script or a truth table is derived from
two sources: `../../taxonomy_all.tsv` (every complete plastome in GenBank with
family and order) and `gb_cache/` (its flat files, one per record). The harness
reads the derived tables by these plain names; `./refresh_references.sh` rebuilds
them into `refresh.new/` and `--install` swaps them in, moving the previous set to
`archive/refs_<date>/`.

| table | built by | what it holds |
|---|---|---|
| `record_meta.tsv` | `extract_record_meta.pl` | submitter group (last author + institution) and RefSeq source per record |
| `all_lengths.tsv` | `extract_guide_lengths.pl` | gene and exon lengths of every cached record (names normalised by `gene_synonyms.tsv`) |
| `../reference_sketches_full.tsv` | `build_reference_sketches.pl --gb-dir` | MinHash sketch of every cached genome: the full guide database |
| `taxonomy_consistency.tsv` | `check_taxonomy_consistency.pl` | sketch-versus-declared-lineage verdicts; MISLABELED records never guide or vote |
| `lineage_assignments.tsv` | `assign_lineages.pl` | lineages assigned to records that declared none (audit trail) |
| `species_units.tsv`, `species_lengths.tsv`, `species_consensus.tsv` | `build_species_consensus.pl` | one vote per submitter group, one form per species; the representative record of each species |
| `gene_expect.tsv`, `gene_profile.tsv` | `build_expect_from_lengths.pl` | species-voted length expectations by genus / family / order / all |
| `lineage_species_counts.tsv` | `refresh_references.sh` | species per family and order, the presence prior's denominator |
| `trna_window_residuals.tsv`, `trna_exon_lengths.tsv` | curated by hand | tRNA boundary convention used by the ARAGORN refinement and `../recurate_trna.pl` |

Superseded versions of these tables, the record-voted profile pipeline's
intermediates (`cds_rows*.tsv`, `guide_quality*.tsv`), truth-table backups from the
tRNA recuration and every result directory not listed under "Current reference
numbers" live in `archive/`. Nothing reads from there; it can be deleted.

## Reading the output

`results/summary.tsv` has one row per genome. The columns that matter for
boundary work:

- `both_exact%` — both ends land exactly on the GenBank coordinates
- `5p_exact%` / `3p_exact%` — each end scored separately
- `5p_median` / `3p_median` — **signed** median offset in nt

Offsets are signed relative to the direction of transcription, on both strands:

- **negative** = predicted boundary lies *upstream* of truth (gene called too long)
- **positive** = predicted boundary lies *downstream* of truth (gene called too short)

A consistently positive `5p_median` means starts are being called too far into
the gene — the signature of picking the wrong start codon.

Each genome also gets `<genome>.detail.tsv`, sorted worst-offset-first, with the
truth and predicted coordinates side by side. That is the file to open when a
particular gene keeps needing hand correction.

## Why not `genbank2verdant.pl`

That script is for importing annotations, not for benchmarking. It hardcodes
`"+"` for every feature, so it cannot represent strand; it flattens `join()`
exons without numbering them; and its `/gene/ || /tRNA/ || /rRNA/` test matches
those substrings anywhere on the line. `gb_to_truth.pl` reads strand from
`complement()`, handles both `complement(join(...))` and
`join(complement(...),...)` (needed for trans-spliced *rps12*), recovers rRNA
names from `/product` when `/gene` is absent, and normalises tRNA naming
(`tRNA-Ala (UGC)`, `trnA`, `trnA-UGC` all become `trnA-UGC`).

Features carrying only a `locus_tag` — hypothetical ORFs such as
`/note="ORF105"` — are skipped. No reference set can transfer a name to them, so
counting them as misses would understate recall, and leaving them in a guide
file injects junk reference names.

## Measured results

Leave-one-out over the 9-genome set (1331 truth features), current vs the code as found:

| | recall | 5' exact | 3' exact | both exact | wall clock |
|---|---|---|---|---|---|
| baseline (as found) | 79.1% | 86.7% | 92.3% | 81.5% | 1280s |
| current | 79.2% | **91.8%** | **94.2%** | **87.3%** | **177s** |

Same recall, boundary accuracy up 5.8 points, 7.2x faster.

**Recall is limited by reference availability, not by the algorithm.** Genomes with
a close relative in the guide set do well; genomes without one do not, and no
boundary work fixes that because the gene is never proposed:

| target | closest guide (k-mer) | recall |
|---|---|---|
| Zea mays | Sorghum, 0.806 | 98.7% |
| Oryza sativa | Sorghum, 0.397 | 95.3% |
| Sorghum bicolor | Zea, 0.803 | 88.5% |
| Nicotiana tabacum | Daucus, 0.214 | 87.7% |
| Pinus thunbergii | Nicotiana, 0.042 | 44.8% |
| Marchantia paleacea | Nicotiana, 0.025 | 40.7% |

Zea scored 59.7% recall when the test set held no other grass and 98.7% once
Sorghum and Oryza were added. Expanding the reference database is worth more than
any remaining algorithmic change.

## Current reference numbers

The arms below are the ones every later change is measured against. All nine
truth tables (and the 12 Asparagaceae and 30 grass tables in `asp/` and
`curated/`) sit on one tRNA convention (`../recurate_trna.pl`, ARAGORN 1.2.41).
"both exact" counts features with both boundaries correct, out of 1,331 truth
features.

| arm | command | recall / precision | both exact |
|---|---|---|---|
| leave-one-out, guides = the other 8 test genomes | `./run_regression.sh --outdir results_loo_gen` | | 1150 |
| curated DB (`../reference_sketches.tsv`), 8 guides | `./run_regression.sh --db ../reference_sketches.tsv --nguides 8 --outdir results_db8_oldnew_cl` | | 1168 |
| full DB (`../reference_sketches_full.tsv`), 8 guides | `./run_regression.sh --db ../reference_sketches_full.tsv --nguides 8 --outdir results_db8_full_cl` | 97.2 / 93.6 | 1208 |
| hybrid: curated DB when it has 3 same-family guides, else full DB | `./run_regression.sh --db ../reference_sketches.tsv --db2 ../reference_sketches_full.tsv --nguides 8 --outdir results_db8_hybrid_gen` | 98.1 / 94.8 | 1246 |
| full DB, no guide from the target's genus | add `--exclude-genus` (`results_db8_full_nogenus`) | 97.1 / 94.0 | 1198 |
| hybrid, no guide from the target's genus | add `--exclude-genus` (`results_db8_hybrid_nogenus`) | 98.0 / 95.3 | 1226 |

The genus-excluded arms are the numbers to quote for a genome from a lineage the
database has never seen annotated. Two variants were measured and rejected, and
their result directories are kept for the record: ranking guides by within-species
disagreement (`results_db8_hybrid_w0.005`, `_w0.02`: 1235, 1231) and the presence
prior taken from the guides' nearest species instead of the family
(`results_db8_{full,hybrid}_local`: no gain, one more spurious ORF).

## Self-annotation: the upper bound

`./run_regression.sh --self` annotates each genome using **only its own
annotation** as the guide. Reference divergence is removed entirely - the
pipeline is handed the exact answer - so anything short of a perfect result is
the algorithm's own defect. It is the sharpest diagnostic in this directory and
it should be run before chasing anything in the leave-one-out numbers.

| | recall | 5' exact | 3' exact | both exact |
|---|---|---|---|---|
| as found | 94.8% | 96.5% | 97.5% | 94.5% |
| current | **99.2%** | **100.0%** | **99.8%** | **99.8%** |

**No algorithmic boundary errors remain.** All three inexact features are flagged
`inconsistent_copies` - the source record annotates a gene's two IR copies at
different lengths, which no single reference can satisfy. Seven of nine genomes
are exact on every matched feature.

It immediately exposed a bug no leave-one-out run could have isolated:
`get_annotated_regions_fromverdant.pl` keyed its extracted sequences by gene
name, so same-named loci overwrote each other. Arabidopsis has 45 tRNA loci but
only 31 references were being extracted - the three distinct *trnS*
isoacceptors collapsed to one, and 65 tRNAs went missing across the set.

The remaining self-annotation error is almost entirely CDS:

| feature type | error rate with a perfect reference |
|---|---|
| tRNA | 1% (4/391) |
| rRNA | 5% (3/66) |
| CDS | **8% (69/863)** |

Three defects came out of it, all invisible to leave-one-out:

1. **Splice-site refinement was moving correct boundaries.** Handed the true
   intron, it shifted it away whenever that intron lacked a GT donor - which the
   survey above says is a third of them. It is now gated on the closest guide's
   k-mer containment (`ANNOBTD_SPLICE_THRESHOLD`, default 0.5): off when a close
   relative is available and the alignment can be trusted, on when it cannot.
   Worth 3.4 points of self-annotation accuracy, and still worth 0.4 on
   leave-one-out where it does apply.
2. **`nearstart` was one codon late.** It is found by probing gene residue 1, so
   a hit puts residue 0 one residue earlier. The fallback path subtracted its own
   offset but the initial probe subtracted nothing, putting every affected start
   +3 nt downstream - *ndhE*, *ycf1*, *ndhH*, *psaA*.
3. **The start codon was chosen by "whichever is further along".** `match_start`
   (an ATG) and `alt_start` (the reference's own start codon) were compared with
   `>=`. Genes whose real start is not ATG - *ycf1* - took a downstream ATG
   instead. Both candidates are now scored by distance from the reference-implied
   position.

A fourth defect was in the harness rather than annoBTD, and had been distorting
every number in this file: `compare_annotations.pl` walked predictions in
coordinate order and let the first one claim a truth row. Where a gene is IR
duplicated and the prediction set holds an extra copy, the spurious low-coordinate
call took the match and the *correct* call was scored tens of thousands of
nucleotides out. Pairing is now assigned nearest-pair-first, which is
order independent. Eight of the nine apparent wrong-IR-copy errors were this.

The ninth was real, and closed the last item from the original diagnostic: the
reverse-strand proximity search in the short-exon path never updated `$distance`,
so its comparison always succeeded and the LAST candidate won instead of the
nearest. Arabidopsis *rpl16_exon1* is a 9 nt exon whose sequence occurs twice in
the genome, and it was being placed 70 kb away from its own exon 2.

A fifth defect, found by pulling on *psaI* and *psbJ*: the ORF overlap filter was
discarding the ORF that matched the gene. Overlapping ORFs on one strand are in
different reading frames, so a covered ORF can be a different real gene rather
than a redundant fragment - Arabidopsis psbJ's ORF matched its truth coordinates
exactly and was dropped for a slightly longer neighbour in another frame. The
filter now also requires a candidate to be substantially shorter than whatever
covers it (`ANNOBTD_MIN_LEN_RATIO`, default 0.6). That is worth ~750 extra ORFs
out of 11,500 - a few seconds - and it recovered psbJ in both genomes and psaI in
two of three.

A sixth defect, found by pulling on *atpF*: the reference N-terminus was located
in the ORF protein with a **three amino acid probe**. Three residues are not
specific - in a 150-residue ORF a given 3-mer turns up by chance - and a spurious
early hit pushed `$nearstart` tens of codons upstream. Everything downstream is
bounded by it, since `match_start` and `alt_start` only search at or before
`$nearstart*3`, so the real start became unreachable: *atpF* exon 2 lost 57 nt in
Marchantia and 22 nt in Pinus. `locate_ref_start()` now takes the longest probe
that matches and prefers one that matches exactly once.

Subtracting `$mod_shift` from the forward 3' key was also tried, since it fixes
Marchantia *atpF_exon1*. It is wrong in general - it introduced a systematic -2 nt
on *petN*, *psbZ*, *psaJ* and *rpl32*, costing 0.8 points - and was reverted.

A seventh: the short-exon search chose among candidates by proximity alone.
*petB* exon 1 is 6 nt, so its sequence occurs 42 times in the Spinacia genome, and
the wrong candidate sat 47 nt *closer* to exon 2 than the right one. The intron
following a real short exon opens with GT in all 27 short-exon introns of the
reference set (*petB*, *petD*, *rpl16* across nine genomes), so `pick_short_exon()`
now prefers candidates carrying that donor, with proximity as the tie-break and a
fall back to plain proximity when no candidate has one.

Note the contrast with the splice refinement above, which had to be gated off:
there, the alignment had already determined the boundary and the splice signal
was second-guessing it. Here the alignment gives no answer at all - a 6-mer
matches 42 places - so the same signal is the only information available.

An eighth, and the last of the 5' errors: the forward start coordinate was
`Start + $start_match - $for_frame + 1 + $mod_shift`. A position in the
frame-shifted ORF sequence actually sits at `Start + $for_frame + $p + 1`. The two
agree for frames 0 and 1 - stop-to-stop ORFs have length divisible by 3, so
`$mod_shift` is `(3-$for_frame)%3` and the terms cancel - but at frame 2 the old
form lands 3 nt early. Daucus *psaI* was the only frame-2 case in the set, which
is why it survived every earlier round.

Of the 4 features still inexact under a perfect reference, **3 are not annotator
errors at all**. Oryza *trnH*, Spinacia *trnA_exon1* and Spinacia *rrn23* each have
two IR copies carrying **identical sequence** but annotated at different lengths in
the source record - 75 vs 78 nt, 38 vs 36, 2810 vs 2811, one always a prefix of
the other. A reference set holds one sequence per gene, so whichever extent it
takes, the other copy scores as an error. annoBTD gives both copies the same
extent, which is self-consistent and arguably better than the record.

`gb_to_truth.pl` now marks these `inconsistent_copies` in the flags column and
`compare_annotations.pl` reports them separately, so the ceiling stays visible
instead of looking like a defect. Six features of 1331 (0.5%) are affected.

The last algorithmic error, Marchantia *atpF_exon1*, was a units mismatch.
`$end_match` reaches the coordinate formulas from two paths on **different
scales**: when the reference ends in a stop codon the callers measure on
`length($orf_seq)` - nucleotides - but otherwise `match_end()` measures on
`length($orf_prot)*3`, which drops any trailing partial codon. The callers treated
both identically. The reverse branches happened to carry a `+$mod_shift` that
cancelled it; the forward branches did not.

Both paths are now normalised onto the nucleotide scale at the call site, and the
reverse branches' compensating `+$mod_shift` was removed with it. Correcting the
units inside `match_end()` instead was tried first and was wrong - it double
counts on the reverse branches, taking the algorithmic error count from 1 to 10.

## The reference database

Two sketch databases exist. `../reference_sketches.tsv` is the curated pool:
`taxonomy_cache.tsv` (7,912 structurally verified GenBank plastomes with tribe /
subfamily / family / order). `../reference_sketches_full.tsv` is every complete
plastome in GenBank (67,461 records at the last refresh), built from the flat-file
cache by `refresh_references.sh`; guides from it pass the lineage-consistency,
species-collapse and taxonomy-consistency filters described in `run_regression.sh`.
The hybrid arm above uses the curated pool when it holds three same-family guides
and the full one otherwise. Three scripts use a pool:

```bash
# once: sketch the pool (resumable, ~1 accession/sec against NCBI's rate limit)
perl build_reference_sketches.pl taxonomy_cache.tsv reference_sketches.tsv

# per genome: rank the pool and fetch the closest guides
perl select_guides.pl target.fsa reference_sketches.tsv -n 3 --one-per-genus > ranking.tsv
guides=$(./fetch_guides.sh ranking.tsv ./sequenceFiles ./files)
./run_annotation_pipeline.sh target.fsa "$guides" MySpecies
```

Holding exact k-mer sets for ~8,000 genomes would mean ~1.2 billion k-mers, so
`build_reference_sketches.pl` stores a bottom-N MinHash sketch instead - 2,000
hashes per genome, a few megabytes for the whole pool. Measured against exact
Jaccard on the nine test genomes the sketch is within 0.021, and it picks the same
top guide for eight of nine (the exception is Marchantia, whose neighbours are all
at noise level anyway).

Selection recovers phylogeny without consulting the taxonomy columns - they are
carried through only so the choice can be sanity checked:

| target | closest guide chosen | jaccard |
|---|---|---|
| *Nicotiana tabacum* | *Nicotiana sylvestris* | 1.000 |
| *Zea mays* | *Miscanthus sacchariflorus* | 0.701 |
| *Daucus carota* | *Elaeoselinum gummiferum* | 0.580 |
| *Arabidopsis thaliana* | *Catolobus pendulus* | 0.567 |
| *Pinus thunbergii* | (nothing above noise) | 0.022 |

*N. tabacum* pulling *N. sylvestris* at 1.000 is correct, not a bug: tobacco's
plastome is maternally inherited from that parent.

### Measured effect

Leave-one-out over the 9-genome set, guides from the other eight test genomes vs
guides chosen from the database:

| | recall | exact features | spurious |
|---|---|---|---|
| 8 other test genomes | 79.2% | 922 / 1331 (69.3%) | 70 |
| 3 database guides | 83.2% | 975 / 1331 (73.3%) | 101 |
| 3 guides + reference screening | **84.1%** | **988 / 1331 (74.2%)** | 87 |

Restricted to angiosperms - the database currently holds no gymnosperm or
bryophyte, so *Pinus* and *Marchantia* cannot benefit:

| | recall | exact features |
|---|---|---|
| 8 other test genomes | 88.0% | 821 / 1066 (77.0%) |
| 3 database guides | **94.7%** | **890 / 1066 (83.5%)** |

*Arabidopsis* goes from 105 exact features to 141, *Nicotiana* 121 to 143.

**Three guides, not more.** Sweeping the count showed each guide past the third
adds spurious calls faster than it adds information - *Arabidopsis* scores 132
exact at n=3 and 106 at n=6. Spurious calls still rise against the old set
(70 to 101) and that is the open problem with this stage.

Two of the losses are artifacts of the comparison rather than the method:
*Sorghum* and *Oryza* had each other and *Zea* as guides in the leave-one-out set,
which is about as close as a guide can get, while the sketched subset here is only
486 of the 7,912 available. Building the full database should recover those.

## Quality checks

Two scripts diagnose an annotation without needing a truth set, so they work on
real submissions.

### check_annotation.pl - does the called protein make sense?

```bash
perl check_annotation.pl plastome.fsa MySpecies_VERDANT_cleaned_annotation.txt \
     --refs MySpecies_annotated_regions_fromverdant_genes.fsa
```

Joins each gene's exons in transcription order, translates, and reports genes whose
protein is malformed: length not a multiple of three, no recognised start codon,
no terminal stop, internal stops, or a length far from the guide consensus. A
boundary can be wrong in ways coordinates alone never show; translating what was
actually called catches it.

It discriminates cleanly. On GenBank's own annotation of the test genomes it
reports 0-4 problems; on annoBTD's output for the same genomes, 1-8. *rps12* is
skipped - it is the one trans-spliced plastid gene, so joining its exons by
proximity is meaningless.

### check_references.pl - is a guide's own annotation sound?

```bash
perl check_references.pl MySpecies_annotated_regions_fromverdant_genes.fsa \
     --tsv problems.tsv --drop cleaned.fsa
```

GenBank holds real annotation errors, and annoBTD transfers whatever the winning
reference says, so a reference whose own CDS does not translate is a direct route
to a bad call. Two independent signals: **intrinsic** (the reference's own protein
is malformed) and **comparative** (it disagrees with the other guides for that
gene, in length or in 5-mer identity). The comparative test is judged against the
typical similarity within that gene rather than a fixed threshold, so one bad copy
does not drag its innocent peers down with it.

Real examples from the sampled database:

```
ndhH   198 nt vs peer median 1182 (0.17x)
psbM   207 nt vs peer median 105 (1.97x); shares 0% of 5-mers with its peers, against 47% typical
psbN   294 nt vs peer median 132 (2.23x); shares 0% of 5-mers with its peers, against 50% typical
```

The pipeline runs this automatically before BLAST and drops what it flags; set
`ANNOBTD_SCREEN_REFS=0` to disable. Measured on Oryza, dropping 7 bad references
moved recall 93.9% to 98.0%, spurious 11 to 7, exact features 124 to 130.

## Experiment: solving a locus de novo

`denovo_exons.pl` takes a window, a reference **protein**, and solves the internal
exon structure - using nothing about the guide's own exon boundaries.

```bash
perl denovo_exons.pl plastome.fsa 44213 46192 - ycf3_protein.fa --exons 3
```

The structure is constrained rather than copied: the joined exons must translate
close to the reference, introns need a GT donor and a plausible length, and the
frame must run continuously with a start, a terminal stop and no internal stops.
Among structures the protein cannot separate, a canonical A-ending acceptor breaks
the tie.

`test/denovo_benchmark.pl` runs it over every multi-exon protein-coding gene in
the test set, taking the reference protein from the closest other genome.

| | loci | share |
|---|---|---|
| exact match to GenBank | 51 | 61% |
| different coordinates, **identical protein** | 9 | 11% |
| differs | 13 | 16% |
| no solution | 10 | 12% |
| **correct protein** | **60** | **72%** |

Excluding *Pinus* and *Marchantia*, whose nearest available relative is 0.02
jaccard away, it is 58 of 70 (83%). *Zea* is 9/9.

### What the experiment established

**A different structure is not automatically a worse one.** For *ycf3* in
Arabidopsis, Daucus and Sorghum the solver puts the junction 2 nt from GenBank and
the two encode an **identical protein** - the intron is flanked by a short repeat,
so it can slide. What differs is the splice signal: the solver's placement reads
canonical GT...AC where the record reads GC...GG, and five of the eight genomes
annotate GT...AC for that same intron. The benchmark scores those as failures.

**The binding constraint is the GT donor.** *clpP* scores 0 of 4, and the reason is
not the search or the scoring - the true structure scores *higher* than what the
solver returns (0.622 against 0.605). Its second intron carries a **TG** donor in
four of five genomes and TT in the fifth, never GT, so the correct answer is
unreachable by construction. Making the donor a scored preference rather than a
requirement would fix clpP at the cost of a much larger search.

**Two of my own assumptions had to go.** A minimum exon length of 20 nt made
*petB*, *petD* and *rpl16* - whose first exons are 6, 8 and 9 nt - return nothing
at all, 0 of 27; allowing 3 nt took them to 67-100%. And positional identity as
the score rejected the correct Spinacia *ycf3* at 0.05, because that protein
carries three fewer N-terminal residues than tobacco's; an indel-tolerant k-mer
measure fixed it.

### Why this is worth pursuing

The production path matches each exon against its own per-exon reference, so a
guide that decomposes a gene differently transfers the wrong decomposition. That
is the *trnI-GAU* / *rps19* spurious-versus-missing bookkeeping above, and it is
the lineage-specific splicing problem in its original form. Solving the locus as a
unit sidesteps it, because the guide supplies only a protein.

It is not a replacement: it needs the locus anchored first, which the current
pipeline provides, and it needs a reasonably close protein - the two genomes
without one score 1/6 and 1/7.

## Reference selection

Exon structure is lineage specific - grasses differ from tobacco in the exon
counts of several genes - and annoBTD transfers whatever the winning guide says.
`rank_references_by_kmer.pl` ranks guides by canonical k-mer similarity to the
target and the pipeline writes the ranking to `<sp>_reference_ranking.txt`;
`closest_reference()` in the matcher then prefers a relative's copy of a gene.

The ranking recovers phylogeny cleanly, with an order of magnitude between
in-lineage and out-of-lineage:

```
Zea_mays  ->  Sorghum 0.806  Oryza 0.372 | eudicots ~0.07 | Pinus 0.031  Marchantia 0.020
```

An earlier attempt took the *majority* exon length across all guides. That is
exactly wrong when the target's own lineage is a minority in the guide set, and it
was removed. Set `ANNOBTD_DISABLE_RANKING=1` to fall back to score-only selection.

## Determinism

The pipeline used to produce a **different annotation on every run**: three
identical runs on tobacco gave three different output files and both-exact scores
spanning 85.7-86.5%. The cause was unsorted `keys %hash` iteration, which Perl
randomises per process; wherever behaviour depended on that order - overlap
resolution, which ORF wins a reference - the output moved with it.

All 36 such loops now iterate in sorted order, and repeated runs are byte
identical. Anything that compares two runs depends on this: before the fix, any
difference under about a point was indistinguishable from run-to-run noise.

## Debugging a specific coordinate

`match_orfs_to_blast_v2.2.pl` will dump the intermediate values behind every CDS
coordinate it emits:

```bash
ANNOBTD_DEBUG_COORDS=/tmp/coords.tsv ./run_annotation_pipeline.sh ...
```

Each row gives the sub and branch that fired, the ORF bounds, the reading frame,
and the `start_match` / `start_alt` / `end_match` values behind the call. Joining
that against a `.detail.tsv` is how the intron-junction slips were traced to
reference choice rather than to the coordinate arithmetic. The hook costs nothing
when the variable is unset.

## Splice signals in plastids

Plastid group II introns are **GT...AY**, not the spliceosomal GT...AG. Measured
over the 135 cis-spliced introns in the truth set: the donor is GT in 90, and the
acceptor is AC (51) or AT (46) against a single AG. Six are GT...AA - tobacco
*rpl2* among them - so a hard AY test rejects real introns; `splice_score` grades
the acceptor instead. The tRNA introns of *trnI*, *trnA* and *trnL* are a
different structural class and are excluded.

## What is left

- **trans-spliced *rps12*** (18 features) regressed from 44.4% to 33.3%. Its
  exons are not adjacent in the genome, so both the consensus-reference rule and
  the intron refinement reason about it wrongly. It needs handling as its own
  case rather than as a multi-exon gene.
- **1-2 nt slips** remain the largest class, now concentrated in tRNAs, whose
  introns the GT...AY rule deliberately does not touch.
- **>100 nt errors** are wrong-ORF selection, not a boundary problem.
- The overlap-resolution pass (`%temp_store`) is unstable: moving one gene's
  coordinates can change which *other* features survive it, so a local fix shows
  up as unrelated features changing. That instability is worth addressing before
  chasing the remaining few points.
- Untouched from the original diagnostic: `exon_mods_mid` still emits `End + 1`
  for the 3' end of middle exons, `exon_mods_end` / `exon_mods_mid` still have no
  match-quality guard, the reverse-strand `$distance` in the short-exon search is
  still never updated, and circularity is not handled anywhere.

## Unit tests

```bash
perl unit_boundary.pl
```

Exercises `match_start`, `alt_start`, `match_end`, `translate`,
`best_match_frames` and `extend_rna_hit` directly. The subs are lifted out of
`match_orfs_to_blast_v2.2.pl` at run time, so there is no duplicated copy to
drift.

All 33 currently pass. A test may be marked `KNOWN-BAD` to pin a defect that is
understood but not yet fixed: it reports as a known failure rather than a
surprise, and when the fix lands the label comes off and the test must pass.
That is how the out-of-frame `match_start` bug was pinned before it was fixed.
