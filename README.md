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

**From bioconda** (installs the commands `annobtd` and `annobtd-fetch-refdb`, with BLAST+ and ARAGORN):

```bash
conda create -n annobtd -c conda-forge -c bioconda annobtd
conda activate annobtd
annobtd-fetch-refdb                 # reference database, ~310 MB download, into the environment
annobtd plastome.fasta PREFIX       # the same options as annobtd.sh below
```

Set `ANNOBTD_REFDB=/some/dir` before `annobtd-fetch-refdb` and `annobtd` to keep the database
outside the environment. Guide records downloaded from NCBI are cached in the database
directory (`gb_cache/`), so every run shares them.

**From a clone, with conda for the dependencies**:

```bash
git clone https://github.com/mrmckain/annoBTD.git
cd annoBTD
conda env create -f environment.yml
conda activate annobtd
./fetch_refdb.sh                    # reference database, ~310 MB download
```

**Without conda**, the requirements are:

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
./fetch_refdb.sh                                                    # reference database, ~310 MB download
```

BLAST+ must be on `PATH`, or set `BLAST_BIN=/path/to/ncbi-blast/bin`. ARAGORN
(Laslett and Canback 2004, GPL v3) is vendored as source in `third_party/aragorn`
and is also taken from `PATH` if it is installed there; without it the annotator
still runs but skips the tRNA boundary and strand checks that depend on it.

`conda/recipe/` is the source of the bioconda recipe; a new release is a tag plus a
version and checksum update there and in bioconda-recipes.

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
| `<prefix>_homopolymers.tsv` | genes whose reading frame shifts at a homopolymer run: the run's position, its length and the lengths that would restore the frame |
| `<prefix>_homopolymer_reads.tsv` | with `--reads`: each of those runs counted in the reads, with a verdict |
| `<prefix>_corrected.fasta`, `<prefix>_corrected/` | with `--reads`, when the reads call for another run length: the sequence with it applied, and that sequence's annotation |
| `<prefix>_checks/` | protein sanity (`check_annotation.tsv`) and splice-site (`check_splice.tsv`) reports |
| `work/` | the complete run: BLAST output, ORFs, pipeline log |

### Checking the result feature by feature

```bash
perl tools/validate_annotation.pl Mygenus_species_annobtd/Mygenus_species_annotation.txt my_plastome.fasta \
     --guides-dir Mygenus_species_annobtd/work/files --guide-seqs Mygenus_species_annobtd/work/sequenceFiles \
     --family Poaceae --order Poales [--genus Mygenus] [--val upload_files/Mygenus_species.val] --out validation.tsv
```

One row per CDS, tRNA and rRNA. A feature that agrees with the lineage profile,
the guides' copies of the gene and its own translation is `PASS`. Anything else
is re-derived from the guide protein with tblastn and given one of these verdicts:
`ARTIFACT` (the sequence carries a canonical form the annotator missed: a start
one codon off, an exon boundary within 3 nt, a start over-extended into the
neighbouring gene, a remnant copy at the IR boundary; the fix is given in
`tools/verdant_edit.pl` format), `REAL` (the guide protein aligns in two frames,
stops early, or part of the gene is gone: a pseudogene call with the evidence) or
`UNRESOLVED` (ambiguity codes in the span, an IR boundary inside the gene, an exon
that cannot be placed: look at the assembly). A gene whose reading frame shifts at
a homopolymer run and is intact on both sides of it is `FRAMESHIFT` (see below).
Canonical genes the guides carry but the annotation lacks are tested the same way.

Starts, junctions and RNA genes are not taken from one guide on trust, and a genome's own
GenBank record, when it is among the guides, is one vote like any other (ten *Lilium* GenBank
genomes annotated with and without their own record now differ in none of ~1,280 features):
- **Curated references.** `tables/nterm_refs.tsv` and `tables/junction_refs.tsv` hold the first
  30 aa of every protein and the peptides either side of every cis-splice junction of the 51
  hand-curated genomes, eight orders from *Marchantia* to the grasses and Asparagaceae
  (`refdb_build/build_nterm_refs.pl`, `build_junction_refs.pl`; `pipeline/ConservedStart.pm`).
  One order is one vote, however many genomes it has.
- **Starts.** A rare start codon (TTG, ATA, ...) is kept only when at least half of the guides
  begin their own copy with it; one beside an ATG where two orders start the protein moves
  there. A call starting more than 10 codons downstream of where most guides' N-termini align
  moves up (rpl22, ndhD, rpoB called inside the protein), unless two orders already start the
  protein where the call starts and none upstream. A start no order shares is trimmed to where
  two orders start (rbcL at MSPQTE..., not an upstream ATG some records use); a call inside the
  protein is extended (ycf1, which most Liliaceae RefSeq records keep only as a pseudogene).
- **Splice junctions.** A junction a codon or two off changes the protein at the junction, so
  it moves to the donor and acceptor that two or more orders' proteins agree on (petB exon 2
  is 648, 642 or 657 nt in Liliaceae GenBank; cytochrome b6 is MS|KVYDWF... everywhere). A
  first exon too short to align (petB, petD, rpl16) takes the references' length.
- **tRNAs and rRNAs.** Guide tRNAs are named and oriented from their sequence against the
  curated library when they are fetched (`guides/relabel_guide_trnas.pl`: GenBank swaps
  trnG-GCC/trnG-UCC and now and then a strand). The missing exon of a spliced tRNA is placed
  from the guides; an rRNA called on both strands keeps the strand of the rRNA beside it; rRNA
  ends follow the guides' majority. Intron rows are rebuilt from the exons.
- **IR-boundary remnants.** The partial second copy of a gene that crosses an IR boundary
  (ycf1, ~1.2 kb) is not called, as in the curated genomes.
- **Frameshifts at homopolymer runs.** When a gene's reading frame shifts at a run of 7 or
  more identical bases (cemA in the Liliaceae, 9 to 25 A right after the ATG; ndhG in the
  North American *Lilium*, 18 to 21 A inside the gene), the gene is called whole, from its
  start codon to its stop, in frame before the run and in frame after it, and it keeps its
  name. It is treated as functional: chloroplasts read through frameshifts in such runs
  (recoding at an 11-A run in atpB, Malinova et al. 2021, Plant Cell 33:1682; transcriptional
  slippage, Lin et al. 2015, PLoS One 10:e0129396), so the frameshift alone does not make
  the gene a pseudogene. `<prefix>_homopolymers.tsv` lists each one with the run's
  coordinates, length and the
  lengths that would restore the frame. The validator gives the row the verdict `FRAMESHIFT`
  with the run's position, and `tools/annobtd2tbl.pl` writes its CDS as two spans with
  /exception "ribosomal slippage" and a note naming the run, the form GenBank asks for when
  a frameshift is thought to be corrected by slippage. Older calls of such a gene (an
  internal ATG, an upstream TTG/ATA/AAT, a row named `<gene>_pseudogene`) come out
  `ARTIFACT` with the move or rename as the fix.
- **Counting the runs in the reads.** A long run is where an assembly is most likely to be
  off by a base, so check each reported run against the reads:

  ```bash
  annobtd.sh my_plastome.fasta Mygenus_species --reads R1.fastq.gz,R2.fastq.gz [other options]
  ```

  or, for an annotation that already exists,

  ```bash
  perl tools/check_homopolymers.pl my_plastome.fasta Mygenus_species_annobtd/Mygenus_species_homopolymers.tsv \
       --reads R1.fastq.gz R2.fastq.gz --out homopolymer_reads.tsv [--fix corrected.fasta]
  ```

  `<prefix>_homopolymer_reads.tsv` has one row per run: how many reads cross it, the run
  lengths they show, and a verdict. `confirmed`: the assembly's length is the commonest in
  the reads by at least 1.5 times the next. `reads differ`: another length is; it is
  suggested. `unresolved`: two or more lengths are about equally common; if one of them keeps
  the gene in frame it is suggested, since an intact gene is the better guess when the reads
  cannot decide. `too few reads`: fewer than 10 reads cross the run. When a length is
  suggested, `annobtd.sh --reads` writes `<prefix>_corrected.fasta` with it applied and
  annotates that sequence in `<prefix>_corrected/`; use the corrected sequence and its
  annotation from then on. Reads undercount long runs more often than they overcount them,
  and above about 18 bases the commonest length is itself uncertain.

With `--val`, table2asn's validator
messages for the same features are attached. `test/validate_truth.sh` scores the
verdicts against a project's hand edits (see `test/README.md`).

### Submitting to GenBank

`tools/annobtd2tbl.pl` writes the feature table and the modifier-bearing FASTA
that table2asn takes, translating every CDS first (`--check`); `tools/verdant_edit.pl`
applies hand edits (rename, move, add, delete) and rebuilds the spacer rows.

## Layout

| | |
|---|---|
| `annobtd.sh` | annotate one genome (the entry point) |
| `fetch_refdb.sh` | download or build the reference database bundle |
| `pipeline/` | the annotation pipeline: `run_annotation_pipeline.sh` drives `annotate_plastome.pl` (ORFs), `run_multiblast.sh`, `identify_best_ref_for_orf_withscoring_v0.7.pl`, `match_orfs_to_blast_v2.2.pl` (boundaries) and `post_filter_annotation.pl`; `check_annotation.pl` and `check_splice.pl` are its sanity reports |
| `guides/` | choosing guides: `choose_guides.sh` (filters and the hybrid rule), `select_guides.pl` (sketch ranking), `score_guide_consistency.pl`, `fetch_guides.sh` |
| `refdb_build/` | building the reference database from GenBank: the taxonomy table, the flat-file mine, sketches, record metadata, species consensus, the profile, the taxonomy check, lineage assignment, the tRNA library and truth recuration |
| `tables/` | curated tables the annotator reads: the tRNA sequence library and boundary conventions, gene synonyms, plastid tRNA identities |
| `tools/` | format converters and the de novo exon experiment |
| `test/` | regression harness, truth sets, measured results, `refresh_references.sh` |
| `conda/recipe/`, `environment.yml` | bioconda recipe and the run-time environment |
| `third_party/aragorn` | ARAGORN 1.2.41 source |
| `legacy/` | pre-2026 versions of the pipeline scripts, kept for reference |
| `refdb/`, `bin/` | the downloaded database and the built ARAGORN binary (not in git) |

## License

GNU Affero General Public License v3 (AGPL-3.0, see `LICENSE.txt`). The vendored ARAGORN source in `third_party/aragorn` keeps its own GPL v3 license (Dean Laslett and Bjorn Canback); GPL v3 and AGPL v3 are compatible.
