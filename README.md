<p align="center">
  <img src="app/assets/branding/pepscope-logo.png" alt="MSP1 PepScope logo" width="520">
</p>

# MSP1 PepScope (Africa)

> **Prototype disclaimer:** PepScope is a research prototype for exploring MSP1 sequence conservation and diversity across African samples in the MalariaGEN Pf8 dataset. This is an unofficial early-access release intended solely for research and testing (not for clinical, diagnostic, or public-health decision-making). For feedback and questions, please contact jacob.almagro-garcia@lstmed.ac.uk.

PepScope helps researchers explore amino-acid conservation and diversity across the
*Plasmodium falciparum* MSP1 protein. It uses thousands of African AA haplotypes to derive an 
interactive collection of block, amino-acid, mutation, and sliding-peptide-window metrics and summaries.

You can use the app to ask questions such as:

- Which parts of MSP1 are well called and conserved?
- Which amino-acid substitutions are common in Africa or in a particular African population or country?
- How many distinct complete haplotypes or peptides occur in a region?
- How much of a cohort is covered by the most common peptide sequence?
- Where does a peptide sit on the gene protein?

Note: this application contains **aggregate metrics and representative sequences only**. It does
not contain sample identifiers or sample-level haplotype data.

## Data source

The biological data source is the
[MalariaGEN *Plasmodium falciparum* Community Project Pf8 release](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/).
Pf8 is an open worldwide parasite-genome variation dataset. PepScope uses derived MSP1
amino-acid haplotype calls from QC-pass African samples, summarized for Africa as a whole,
African populations, and African countries.

Please cite the source dataset if using the data:

> MalariaGEN et al. (2025). *Pf8: an open dataset of Plasmodium falciparum genome variation in
> 33,325 worldwide samples*. Wellcome Open Research, 10:325.
> [https://doi.org/10.12688/wellcomeopenres.24031.1](https://doi.org/10.12688/wellcomeopenres.24031.1)

Below we use synthetic data to quickly explain how the main metrics are derived (but they are not real Pf8 sample sequences).

## How metrics are derived

The database is built from all QC-pass African samples, with one amino-acid haplotype per MSP1 block. The following toy dataset represents 20 samples in one hypothetical 20-amino-acid block. It includes common variants, mixed calls, incomplete calls, and three samples with poor-quality stretches (but still belonging to QC-pass samples).

The hypothetical 3D7/reference sequence is `MKTAYIAKQRQISFVKSHFS`.

| Sample | 20-locus block haplotype | Interpretation |
|---|---|---|
| S01 | `MKTAYIAKQRQISFVKSHFS` | Complete, matches the reference |
| S02 | `MKTAYIAKQRQISFVKSHFS` | Complete, matches the reference |
| S03 | `MKTAYIAKQRQISFVKSHFS` | Complete, matches the reference |
| S04 | `MKTAYIAKQRQISFVKSHFS` | Complete, matches the reference |
| S05 | `MKTAFIAKQRQISFVKSHFS` | Complete, Y5F variant |
| S06 | `MKTAFIAKQRQISFVKSHFS` | Complete, Y5F variant |
| S07 | `MKTAYIAKQKQISFVKSHFS` | Complete, R10K variant |
| S08 | `MKTAYIAKQKQISFVKSHFS` | Complete, R10K variant |
| S09 | `MKTAYIAKQRQISFIKSHFS` | Complete, V15I variant |
| S10 | `MKTAYIAKQRQISFIKSHFS` | Complete, V15I variant |
| S11 | `MKTAFIAKQKQISFVKSHFS` | Complete, Y5F and R10K variants |
| S12 | `MKTAYIAKQRQISFIKSYFS` | Complete, V15I and H18Y variants |
| S13 | `MKTA[Y/F]IAKQRQISFVKSHFS` | Mixed Y/F call at locus 5 |
| S14 | `MKTAYIAKQ[R/K]QISFVKSHFS` | Mixed R/K call at locus 10 |
| S15 | `MKTAYIAKQRQISF[V/I]KSHFS` | Mixed V/I call at locus 15 |
| S16 | `MKTAYIA-QRQISFVKSHFS` | One missing call at locus 8 |
| S17 | `MKTAYIAKQRQISFVKS!FS` | One unresolved call at locus 18 |
| S18 | `MKTAYIAKQR----------` | **Poor quality:** 10 consecutive missing calls |
| S19 | `-----IAKQRQISFVK----` | **Poor quality:** 9 missing calls in two stretches |
| S20 | `MKT-----QRQISF------` | **Poor quality:** 11 missing calls in two stretches |

### Reading the notation

- A single letter such as `A` or `V` is a called amino acid.
- `-` is a missing call. It is not treated as another amino acid.
- A bracketed value such as `[Y/F]` is a mixed call. It means that both residues are supported at
  that locus. The brackets count as one locus, so every example row still spans 20 loci.
- `!` and `*` mark unresolved calls or haplotypes (usually refering to the whole haplotype but here refering only to a call for simplicity, samples with unresolved haplotypes are rare and always excluded).

### Which haplotypes are complete?

A **complete haplotype** has one called amino acid at every locus: no `-`, no `!` or `*`, and no
mixed calls. In the example, S01–S12 are complete. S13–S20 are not.

The 12 complete rows contain six distinct complete haplotypes:

| Complete haplotype | Samples | Count | Frequency among complete rows |
|---|---|---:|---:|
| `MKTAYIAKQRQISFVKSHFS` | S01–S04 | 4 | 33.3% |
| `MKTAFIAKQRQISFVKSHFS` | S05–S06 | 2 | 16.7% |
| `MKTAYIAKQKQISFVKSHFS` | S07–S08 | 2 | 16.7% |
| `MKTAYIAKQRQISFIKSHFS` | S09–S10 | 2 | 16.7% |
| `MKTAFIAKQKQISFVKSHFS` | S11 | 1 | 8.3% |
| `MKTAYIAKQRQISFIKSYFS` | S12 | 1 | 8.3% |

Complete haplotypes are used for exact sequence comparisons. This avoids interpreting a missing
or mixed call as if it were a fully observed sequence. In this example:

- the top complete haplotype is the reference sequence, at `4 / 12 = 33.3%`
- there are six unique complete haplotypes
- two are singletons, so the singleton-haplotype fraction is `2 / 6 = 33.3%`
- the top three complete haplotypes cover `8 / 12 = 66.7%` of complete observations

### Callability and missingness

For callability, an amino acid or a mixed call is callable, `-`, `!`, and `*` are not.
Across the example there are 400 possible calls (`20 samples × 20 loci`):

- 368 are callable, giving an overall call rate of `368 / 400 = 92.0%`
- 31 are explicitly missing, giving a missing fraction of `31 / 400 = 7.75%`
- one is unresolved, accounting for the remaining `0.25%`
- four samples have at least one `-` call (`20%`)
- three samples have at least five missing calls (`15%`).

PepScope retains callability and missingness alongside diversity. **A region can look conserved
simply because difficult samples or loci are absent**, so the app lets you inspect both kinds of
evidence together.

### How mixed calls contribute to allele frequencies

Mixed calls contribute fractionally instead of becoming an extra amino acid. For
example, `[Y/F]` contributes `0.5` to Y and `0.5` to F. It contributes a total weight of one, just
like an unambiguous call.

At locus 5 in the example:

- 18 of 20 samples are callable, so the locus call rate is `90.0%`
- Y receives a total weight of `14.5`
- F receives a total weight of `3.5`
- the allele frequencies are therefore Y `14.5 / 18 = 80.6%` and F `3.5 / 18 = 19.4%`.

Missing and unresolved calls do not enter the allele denominator. This is why the denominator is
18 rather than 20.

### Reference, consensus, and top sequences

These labels answer related but different questions:

- **3D7/reference sequence:** the fixed reference sequence used for coordinates and comparisons.
- **Top complete haplotype or peptide:** the most frequent fully called sequence actually observed
  among complete sequences.
- **Allelic consensus:** a sequence assembled locus by locus from the most frequent callable amino
  acid. Mixed calls contribute fractionally. Missing and unresolved calls do not contribute.
- **All-call consensus:** a locus-by-locus consensus that also lets `-` contribute as an AA state. If
  missingness is the most frequent observation at a locus, the consensus contains `-`. Mixed calls
  still contribute to their component amino acids but the output never contains a bracketed mixed
  token.
- **Top haplotype (all):** the most frequent full observed pattern when incomplete and mixed
  strings are retained exactly as called (it could contain both missing and mixed calls).

Note: the allelic consensus does not have to be an observed haplotype: its most common residue at each
locus may come from different samples. The top complete haplotype, by contrast, is always an
observed complete sequence.

## Derived metrics and summaries

PepScope stores a collection of pre-computed metrics across different cohorts, loci, gene blocks, and window lengths, so the app can filter and plot
them quickly.

### Block summaries

MSP1 is divided into 17 predefined blocks. For every available geographic group, the database
stores block coordinates and reference sequence, callability and missingness, complete-haplotype
counts, sequence frequencies, diversity, pairwise distance, and top-haplotype coverage.

Some useful block metrics are:

| Metric | Explanation |
|---|---|
| Mean call rate | The average fraction of block loci that are callable in a sample |
| Complete haplotypes | Samples with one called amino acid at every locus in the block |
| Haplotype diversity (Hd) | The chance that two randomly chosen complete haplotypes differ, with a small-sample correction |
| Effective haplotype number | The number of equally common haplotypes that would give the observed Simpson diversity |
| Shannon diversity | Diversity that increases with both the number of haplotypes and how evenly they are represented |
| Top-N coverage | The fraction of complete observations represented by the 1, 2, 3, 5, or 10 most frequent haplotypes |
| Mean pairwise distance | The average number of amino-acid differences between pairs of complete haplotypes |

### Locus summaries

For each amino-acid position, PepScope stores many metrics, among them: call rate, missing and mixed fractions, reference and
consensus frequencies, the second-most-common residue frequency, observed-residue counts,
polymorphism indicators, expected heterozygosity, Shannon entropy, Gini inequality, and Pielou
evenness. These are the default signals shown in the locus tracks (can be expanded to show more metrics).

### Sliding peptide windows

Within each block, windows of 7–30 amino acids move one locus at a time. Each window is treated as
a short peptide haplotype and receives the same kinds of callability, sequence-frequency,
diversity, distance, and coverage summaries. Windows never cross a block boundary.

This lets the filtering page find regions that combine good callability with conservation, or
highlight regions where apparently strong conservation needs to be interpreted cautiously because
many calls are missing.

### Mutation summaries

Mutations use canonical labels such as `Y5F`: reference residue, MSP1 amino-acid coordinate, and
alternate residue. Their frequencies include fractional evidence from mixed calls and are shown
together with the missing fraction at the same locus. Tables can compare the currently selected
group (e.g., country) with the other groups present in the database (e.g., other countries in the same population or Africa as a whole). The app can also show the top complete haplotypes that carry a mutation.

### Sequence logos and protein structure

Sequence logos (mimicked with barplots) summarize residue frequencies at each locus. Mixed calls are divided across their
possible residues, missing and unresolved mass is displayed separately rather than as another
amino acid.

When a reviewed experimental structure or AlphaFold model is available for MSP1, the structure viewer
highlights the selected block or peptide and differences between the 3D7 sequence and the top
complete sequence. Remember that the structure is useful spatial context, not evidence that a substitution changes
function, immunity, or vaccine performance.

## Geographic groups

Every metric is calculated independently for the groups (cohorts) bundled in this release:

- Africa as a whole.
- African population groups.
- African countries.

Changing the cohort selector changes the samples contributing to frequencies and diversity, **it is
not just a display filter**. Small groups should be interpreted carefully because rare variants
and frequency estimates are less stable.

## What you can explore in the app (at a glance)

- **Block overview:** block metadata, quality, haplotype frequencies, diversity, sequence logos,
  mutations, and structure context.
- **Locus tracks:** compact tracks (interactive plots) for callability, allele frequency, polymorphism, and diversity
  along every amino-acid coordinate in a block.
- **Window tracks:** the same kind of overview for sliding peptide windows of a selected length.
- **Region explorer:** all windows overlapping a coordinate or coordinate range.
- **Filter and export:** conservation and quality filters, filter diagnostics, and the dataset with the results (you can download it).
- **About and limitations:** this page give details about the app and the underlying DB, including some advice on interpreting results and highlighting limitations.

## Details for this unofficial prototype release (MSP1)

- Gene: *P. falciparum* MSP1 (`PF3D7_0930300`)
- Reference: 3D7, amino-acid positions 1-1720
- Samples used for metrics: QC-pass African samples (Pf8)
- Geographic levels: Africa, African populations, and African countries
- Blocks: 17
- Peptide windows: 34,140 physical windows, with lengths 7-30.
- Storage: lazy DuckDB views over 73 Zstandard-compressed Parquet files
- Aggregate database size: ~60 MB

(The database and six local mmCIF structure models are read-only)

## Install and run

### Requirements

- Python 3.11
- Approximately 100 MB of project storage, excluding the Python environment
- Writable temporary storage below `/tmp` for bounded DuckDB spill files

From the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

The entry point discovers the adjacent `src/` directory automatically. `PYTHONPATH` is not
required.

### Run the checks

```bash
python -m pip install -r requirements-dev.txt
pytest
ruff check app src tests
```

The test suite validates the Parquet manifest and checksums, structure assets, aggregate counts,
representative database queries, the standalone import boundary, and all Streamlit pages with
Streamlit's headless `AppTest` framework.

## Project layout

```text
app/                 Streamlit entry point, JSON configuration, CSS, and images
data/database/       Immutable aggregate DuckDB-over-Parquet DB
data/structures/     Local mmCIF models
src/aa_metrics_app/  Minimal query, filtering, export, structure, and chart package
tests/               Database, application, and portability checks
```

## Limitations and interpretation

The public database stores aggregate metrics and sequence summaries, not sample-level identifiers.
The application is designed for exploration and hypothesis generation. Missingness, mixed calls,
incomplete haplotypes, cohort size, uneven geographic sampling, reference choice, and structure
coverage can all affect interpretation.

In particular:

- high conservation is not by itself proof of vaccine suitability
- a frequent peptide is not necessarily immunogenic or protective
- structure highlighting does not predict the functional effect of a mutation
- country and population summaries describe the samples available in Pf8, not every parasite
  circulating in those places.

## Data license

The MalariaGEN Pf8 source data
have their own terms and citation requirements, consult the
[official Pf8 data page](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v80/)
before reusing it.
