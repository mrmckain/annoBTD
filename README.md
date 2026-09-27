# annoBTD

Annotation of plastid (chloroplast) genomes from GenBank guide genomes chosen
for each target, with gene boundaries settled by the lineage's annotation
conventions rather than by any single reference. annoBTD is the annotator behind
[Verdant](https://verdant.iplantcollaborative.org).

On nine classic plastomes scored against their GenBank records, guides drawn from
the full database place 98.8% of features with 95.9% precision, and 96.5% of the
matched features have both boundaries exact. The measurements and how to
reproduce them are in [test/README.md](test/README.md).

## Installation

annoBTD is Perl and Bash and needs no Perl modules beyond the core distribution.

Requirements:

| what | why | install |
|---|---|---|
| Perl 5 | all scripts | present on macOS and Linux |
| BLAST+ (`blastn`, `tblastx`, `makeblastdb`) | guide alignment | `brew install blast`, `conda install -c bioconda blast`, `apt install ncbi-blast+` |
| a C compiler (`cc`) | builds ARAGORN once | Xcode command line tools, `apt install build-essential` |
| `curl` | fetches guide records from NCBI and the reference database | present on macOS and Linux |

Then:

```bash
git clone https://github.com/mrmckain/annoBTD.git
cd annoBTD
cc -O2 -w -o bin/aragorn third_party/aragorn/aragorn1.2.41.c -lm    # tRNA boundaries
./fetch_refdb.sh                                                    # reference database, ~350 MB download
```

BLAST+ must be on `PATH`, or set `BLAST_BIN=/path/to/ncbi-blast/bin`. ARAGORN
(Laslett and Canback 2004, GPL v3) is vendored as source in `third_party/aragorn`;
without the binary the annotator still runs but skips the tRNA boundary and strand
checks that depend on it.

Check the install with the built-in regression, which annotates one of the nine
test genomes from the other eight and scores it against its GenBank record (about a
minute, no database needed):

```bash
cd test
./run_regression.sh --only Nicotiana_tabacum --outdir smoke
```

The last lines report recall, precision and the share of features with both
boundaries exact. A run with guides from the reference database, the way real
annotation works, is:

```bash
./run_regression.sh --db ../refdb/reference_sketches.tsv --db2 ../refdb/reference_sketches_full.tsv --nguides 8 --only Nicotiana_tabacum --outdir smoke_db
```

## The reference database

Guides are chosen from a database of every complete plastid genome in GenBank,
and the profile that settles gene boundaries is mined from the same records. The
database is 1.5 GB and is not in git; `fetch_refdb.sh` downloads it as one
checksummed bundle into `refdb/`:

| file | what it is |
|---|---|
| `reference_sketches.tsv` | MinHash sketches of the curated pool: 7,882 of the 7,912 structurally verified plastomes in `test/taxonomy_cache.tsv` |
| `reference_sketches_full.tsv` | sketches of every complete plastome in GenBank (67,461 at the September 2026 build) |
| `taxonomy_all.tsv` | accession, organism, family and order of each record |
| `gene_expect.tsv`, `gene_profile.tsv` | per-gene and per-exon length expectations by genus, family and order, one vote per species and one vote per submitting group within a species |
| `species_units.tsv`, `species_lengths.tsv` | which records are one genome, which record represents each species, and each species' consensus gene lengths |
| `all_lengths.tsv` | every record's gene and exon lengths, for the guide consistency filter |
| `taxonomy_consistency.tsv` | records whose sequence contradicts their declared lineage (never used as guides) |
| `lineage_species_counts.tsv` | species per family and order |

`fetch_refdb.sh --version V` fetches another release, `--url` another source (a
mirror or a local `file:///` path), `--dir` another destination; `ANNOBTD_REFDB`
tells the scripts where the database lives if it is not in `refdb/`. Every table
can be rebuilt from GenBank with `test/refresh_references.sh` (see
[test/README.md](test/README.md)); `fetch_refdb.sh --build V` packs a rebuilt
directory into a new bundle.

## Annotating a genome

```bash
./annobtd.sh my_plastome.fasta Mygenus_species --family Poaceae --order Poales
```

The FASTA holds one circular plastome. Family and order choose the length
expectations and the guide filters; without them, annoBTD takes the lineage of the
three most similar genomes in the database and prints it. `--genus` adds the
genus stratum of the profile. `--db` picks the guide pool: `curated` (the 7,882
verified plastomes), `full` (everything), or the default `hybrid`, which uses the
curated pool when it holds three guides from the target's family and the full
pool otherwise. `-n` sets the number of guides (default 8).

What happens: the target is sketched and ranked against the pool; candidates whose
gene lengths follow another lineage's conventions, that duplicate a species, or
whose sequence contradicts their declared lineage are passed over; the guide
records are fetched from NCBI (cached in `refdb/gb_cache`); the pipeline finds
ORFs, aligns the guides to them and places boundaries; a post-annotation filter
then collapses duplicate RNA calls, names every tRNA by sequence against a
curated library, settles strands and splice junctions from the intron boundary
motif, and extends truncated starts to the lineage's settled length.

Output, in `Mygenus_species_annobtd/`:

| file | content |
|---|---|
| `<prefix>_annotation.txt` | `gene  start  end  strand`, 1-based inclusive; exons as `gene_exonN`, introns as `gene_intronN`, spacers as `geneA~geneB`, plus `LSC`, `SSC`, `IRA`, `IRB`, `FULL` region rows |
| `<prefix>_guides.tsv` | the guides used, their similarity, and why each passed or was backfilled |
| `<prefix>_post_filter.txt` | every decision the post filter took, one line each |
| `<prefix>_checks/` | protein sanity (`check_annotation.tsv`) and splice-site (`check_splice.tsv`) reports |
| `work/` | the complete run: BLAST output, ORFs, pipeline log |

## Layout

| | |
|---|---|
| `annobtd.sh` | annotate one genome (entry point) |
| `fetch_refdb.sh` | download or build the reference database bundle |
| `choose_guides.sh`, `select_guides.pl`, `fetch_guides.sh` | guide ranking, filtering and retrieval |
| `run_annotation_pipeline.sh` | the annotation pipeline proper: `annotate_plastome.pl` (ORFs), `run_multiblast.sh`, `identify_best_ref_for_orf_withscoring_v0.7.pl`, `match_orfs_to_blast_v2.2.pl` (boundaries), `post_filter_annotation.pl` |
| `build_*.pl`, `extract_*.pl`, `check_taxonomy_consistency.pl`, `assign_lineages.pl`, `mine_genbank_batch.sh` | building the reference database from GenBank |
| `trna_library.fasta`, `trna_window_residuals.tsv`, `trna_exon_lengths.tsv`, `gene_synonyms.tsv`, `plastid_trna_names.txt` | curated tables the annotator reads |
| `test/` | regression harness, truth sets, measured results |
| `third_party/aragorn` | ARAGORN 1.2.41 source |

## License

GPL v3 (see `LICENSE.txt`). ARAGORN is GPL v3, by Dean Laslett and Bjorn Canback.
