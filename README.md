# SEEDscraper 🌰

Download protein seed sequences from UniProtKB using interactive queries or a query
TSV. SEEDscraper is a standalone tool: it does not run GERMINATE or download a
reference search database.

## Installation

Use Python 3.10 or newer (tested with Python 3.12):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The pinned requirements include Requests and its dependencies. BeautifulSoup and
Biopython are no longer needed. The FASTA reader validates the UniProt headers and
wrapped protein sequences returned by this API; it is not a general FASTA importer.

## Interactive use

```bash
python src/SEEDscraper.py --out results/run1
```

Enter a UniProt query, then a unique gene/output label. For example:

```text
UniProt query: gene:murA AND taxonomy_id:2 AND reviewed:true
Gene/output label: murA
```

Enter more query/label pairs, then type `EXIT` at the query prompt to begin downloading.
Blank queries are ignored. Ctrl-C cancels the run. Labels must start with a letter
or digit and contain only letters, digits, underscores, periods, or hyphens.

The example restricts results to reviewed bacterial entries. Such restrictions are
not added automatically: the query controls the search, and downloaded hits still
require evaluation for suitability as seeds.

## Repeatable batch use

Create a tab-separated file with exactly the columns `gene` and `query`. A small
accession-based download example is provided as `example-data/queries.tsv`:

```bash
python src/SEEDscraper.py --queries example-data/queries.tsv --out results/example --per-gene
```

Each label must be unique, ignoring case. Use a combined UniProt query when you want
multiple search conditions associated with one label.

## Outputs and reruns

Every run requires a **new output directory**; existing directories are rejected.
Nothing is appended to an earlier database or overwritten.

| File | Contents |
|---|---|
| `all_sequences.fasta` | All unique accessions, sorted by accession, with original headers |
| `proteins.db` | SQLite tables for proteins, queries, and their associations |
| `manifest.json` | Completion status, UTC retrieval start, queries, counts, API URL, and UniProt release metadata when provided |
| `seeds/{gene}_seeds.faa` | Optional per-query export, enabled with `--per-gene` |

The `proteins` table stores `accession`, `entry_name`, full `protein_name`, full
`organism`, `sequence`, and original `header`. The `queries` and `query_proteins`
tables preserve each query and which accessions it found. Overlapping queries
share a single protein record. This is accession deduplication, not sequence-identity
clustering: different accessions with identical sequences are retained.

Per-gene exports can be supplied manually to GERMINATE or another tool. A zero-hit
query is reported explicitly and produces an empty per-gene FASTA. Check such files
before using them downstream; GERMINATE rejects empty seed inputs.

## Download reliability

- Follows pagination links and checks unique counts against UniProt totals when provided.
- Uses 10-second connection and 60-second read timeouts, with up to three retries
  for connection failures and retryable HTTP statuses, including rate limiting.
- Stops if release metadata changes during the run, counts disagree, a response is
  malformed, or an accession has conflicting records.
- Saves pages to SQLite as they arrive. Final outputs are published together only
  after all queries succeed.
- On failure or cancellation during downloading, a hidden `.NAME.partial-*`
  directory retains downloaded records and an `incomplete` manifest for inspection.
  It is not a completed dataset and cannot be automatically resumed. Rerun the
  command to download again. Unexpected process termination may leave a partial
  directory without a manifest.

The default output is `results`. Use `--out` to choose another unused directory.
The SQLite schema differs from the original script; existing databases are not migrated.
The original `example-data/all_sequences.fasta` and `example-data/proteins.db` are
historical artifacts. The old database contains truncated metadata and is not an
expected-output fixture for this version. Regenerate examples using the command above.

## Tests

```bash
python -m unittest discover -s tests -v
```

The offline suite checks metadata, wrapped FASTA, pagination, overlapping queries,
zero matches, failures, output protection, query validation, and retry settings.
It does not assess biological seed quality.

UniProt references: [query API](https://www.uniprot.org/help/api_queries) and
[FASTA headers](https://www.uniprot.org/help/fasta-headers).
