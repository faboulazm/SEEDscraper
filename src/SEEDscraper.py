#!/usr/bin/env python3
"""Download complete UniProtKB query results into a consistent, traceable dataset."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

API = 'https://rest.uniprot.org/uniprotkb/search'
HEADER = re.compile(r'^(?:sp|tr)\|([^|\s]+)\|([^\s]+) (.+?) OS=(.+?)(?= (?:OX|GN|PE|SV)=|$)')


def parse_fasta(text):
    """Parse wrapped UniProt FASTA, retaining complete names and original headers."""
    records = []
    header, chunks = None, []

    def finish():
        if header is None:
            return
        match = HEADER.match(header)
        sequence = ''.join(chunks)
        if not match or not sequence or not re.fullmatch(r'[A-Za-z*]+', sequence):
            raise ValueError('Invalid UniProt FASTA record: ' + header[:120])
        accession, entry_name, protein_name, organism = match.groups()
        records.append((accession, entry_name, protein_name, organism, sequence, header))

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            finish()
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError('Response is not FASTA.')
        else:
            chunks.append(line)
    finish()
    return records


def make_session():
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
                    allowed_methods=['GET'], respect_retry_after_header=True)
    session.mount('https://', HTTPAdapter(max_retries=retries))
    session.headers['User-Agent'] = 'SEEDscraper/0.2'
    return session


def read_queries(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['gene', 'query']:
            raise ValueError('Query TSV must have exactly two columns: gene and query.')
        rows = []
        for row in reader:
            if None in row or row.get('query') is None:
                raise ValueError('Every TSV row must contain a gene and query.')
            rows.append((row['gene'].strip(), row['query'].strip()))
    return rows


def validate_queries(queries):
    if not queries:
        raise ValueError('No queries supplied.')
    seen = set()
    for gene, query in queries:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', gene):
            raise ValueError('Gene labels must start with a letter or digit and use only letters, digits, _, -, or .')
        if gene.casefold() in seen:
            raise ValueError('Gene labels must be unique: ' + gene)
        if not query or any(ch in query for ch in '\r\n\t'):
            raise ValueError('Queries must be nonempty single lines.')
        seen.add(gene.casefold())


def write_fasta(path, records):
    with path.open('w', encoding='utf-8') as handle:
        for sequence, header in records:
            handle.write('>' + header + '\n')
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + '\n')


def download(queries, output, per_gene=False, session=None):
    validate_queries(queries)
    output = Path(output).absolute()
    if output.exists():
        raise ValueError('Output already exists; choose a new --out directory: ' + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.' + output.name + '.partial-', dir=output.parent))
    manifest = {'status': 'incomplete', 'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
                'api': API, 'queries': [], 'unique_proteins': 0}
    own_session = session is None
    session = session or make_session()
    conn = sqlite3.connect(staging / 'proteins.db')
    try:
        conn.executescript('''
            CREATE TABLE proteins (accession TEXT PRIMARY KEY, entry_name TEXT NOT NULL,
                protein_name TEXT NOT NULL, organism TEXT NOT NULL, sequence TEXT NOT NULL,
                header TEXT NOT NULL);
            CREATE TABLE queries (gene TEXT PRIMARY KEY, query TEXT NOT NULL);
            CREATE TABLE query_proteins (gene TEXT NOT NULL REFERENCES queries(gene),
                accession TEXT NOT NULL REFERENCES proteins(accession), PRIMARY KEY(gene, accession));
        ''')
        release_seen = None
        for gene, query in queries:
            conn.execute('INSERT INTO queries VALUES (?, ?)', (gene, query))
            info = {'gene': gene, 'query': query, 'pages': 0, 'records_received': 0, 'unique_proteins': 0}
            manifest['queries'].append(info)
            url, params = API, {'query': query, 'format': 'fasta', 'size': 500}
            visited = set()
            while url:
                parsed = urlparse(url)
                if parsed.scheme != 'https' or parsed.netloc != 'rest.uniprot.org':
                    raise ValueError('Unexpected pagination destination.')
                if url in visited:
                    raise ValueError('Repeated pagination link; download stopped.')
                visited.add(url)
                response = session.get(url, params=params, timeout=(10, 60))
                response.raise_for_status()
                records = parse_fasta(response.text)
                release = response.headers.get('X-UniProt-Release')
                if release_seen and release and release != release_seen:
                    raise ValueError('UniProt release changed during download; rerun for a consistent dataset.')
                release_seen = release or release_seen
                info['uniprot_release'] = release or info.get('uniprot_release')
                info['uniprot_release_date'] = response.headers.get('X-UniProt-Release-Date') or info.get('uniprot_release_date')
                total = response.headers.get('X-Total-Results')
                if total is not None:
                    total = int(total)
                    if 'expected_results' in info and info['expected_results'] != total:
                        raise ValueError('Result count changed during pagination.')
                    info['expected_results'] = total
                for record in records:
                    previous = conn.execute('SELECT * FROM proteins WHERE accession=?', (record[0],)).fetchone()
                    if previous and previous != record:
                        raise ValueError('Conflicting records for accession ' + record[0])
                    conn.execute('INSERT OR IGNORE INTO proteins VALUES (?, ?, ?, ?, ?, ?)', record)
                    conn.execute('INSERT OR IGNORE INTO query_proteins VALUES (?, ?)', (gene, record[0]))
                conn.commit()
                info['pages'] += 1
                info['records_received'] += len(records)
                url = response.links.get('next', {}).get('url')
                params = None
            count = conn.execute('SELECT count(*) FROM query_proteins WHERE gene=?', (gene,)).fetchone()[0]
            info['unique_proteins'] = count
            if 'expected_results' in info and count != info['expected_results']:
                raise ValueError('Downloaded count does not match UniProt total for ' + gene)
            print(f'{gene}: {count} unique proteins ({info["pages"]} pages)' + (' — no matches' if not count else ''))
        write_fasta(staging / 'all_sequences.fasta', conn.execute('SELECT sequence, header FROM proteins ORDER BY accession'))
        if per_gene:
            (staging / 'seeds').mkdir()
            for gene, _ in queries:
                write_fasta(staging / 'seeds' / (gene + '_seeds.faa'), conn.execute(
                    'SELECT p.sequence, p.header FROM proteins p JOIN query_proteins q USING(accession) WHERE q.gene=? ORDER BY p.accession', (gene,)))
        manifest['unique_proteins'] = conn.execute('SELECT count(*) FROM proteins').fetchone()[0]
        manifest['status'] = 'complete'
        conn.close()
        (staging / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        # Reserve the destination exclusively, including against another simultaneous run.
        output.mkdir()
        staging.rename(output)
        print(f'Saved {manifest["unique_proteins"]} unique proteins to {output}')
        return manifest
    except BaseException as exc:
        conn.close()
        manifest['status'] = 'incomplete'
        manifest['error'] = str(exc) or type(exc).__name__
        (staging / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        print(f'Incomplete download retained at {staging}', file=sys.stderr)
        raise
    finally:
        if own_session:
            session.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queries', type=Path, help='TSV with gene and query columns')
    parser.add_argument('--out', type=Path, default=Path('results'), help='New output directory (default: results)')
    parser.add_argument('--per-gene', action='store_true', help='Also export seeds/{gene}_seeds.faa')
    args = parser.parse_args(argv)
    try:
        if args.queries:
            queries = read_queries(args.queries)
        else:
            queries = []
            print('SEEDscraper: enter UniProt queries. Type EXIT to download, or Ctrl-C to cancel.')
            while True:
                query = input('UniProt query: ').strip()
                if query.upper() == 'EXIT':
                    break
                if not query:
                    continue
                gene = input('Gene/output label: ').strip()
                queries.append((gene, query))
        download(queries, args.out, args.per_gene)
        return 0
    except (ValueError, OSError, sqlite3.Error, requests.RequestException) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1
    except (KeyboardInterrupt, EOFError):
        print('\nCancelled.', file=sys.stderr)
        return 130


if __name__ == '__main__':
    sys.exit(main())
