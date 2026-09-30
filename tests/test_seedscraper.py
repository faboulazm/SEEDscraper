import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import requests

spec = importlib.util.spec_from_file_location('scraper', Path(__file__).resolve().parents[1] / 'src' / 'SEEDscraper.py')
scraper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scraper)
FASTA = '>sp|P12345|MURA_TEST UDP-N-acetylglucosamine 1-carboxyvinyltransferase OS=Streptococcus pneumoniae (strain D39) OX=373153 GN=murA PE=1 SV=1\nMPEP\nTIDE\n'
SECOND = FASTA.replace('P12345', 'P23456').replace('MURA_TEST', 'MURB_TEST')

class Response:
    def __init__(self, text=FASTA, next_url=None, total=1, error=None, release='2026_01'):
        self.text = text
        self.headers = {'X-Total-Results': str(total), 'X-UniProt-Release': release}
        self.links = {'next': {'url': next_url}} if next_url else {}
        self.error = error
    def raise_for_status(self):
        if self.error:
            raise self.error

class Session:
    def __init__(self, *responses):
        self.responses = iter(responses)
        self.calls = []
    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        response = next(self.responses)
        if isinstance(response, Exception):
            raise response
        return response

class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / 'results'
    def test_complete_metadata_and_wrapped_sequence(self):
        r = scraper.parse_fasta(FASTA)[0]
        self.assertEqual(r[:5], ('P12345', 'MURA_TEST', 'UDP-N-acetylglucosamine 1-carboxyvinyltransferase', 'Streptococcus pneumoniae (strain D39)', 'MPEPTIDE'))
    def test_pagination_overlap_and_exports(self):
        session = Session(Response(next_url=scraper.API + '?cursor=abc', total=2), Response(SECOND, total=2), Response())
        manifest = scraper.download([('murA', 'query one'), ('murB', 'query two')], self.out, True, session)
        self.assertEqual(manifest['unique_proteins'], 2)
        self.assertEqual(session.calls[1][1]['params'], None)
        self.assertEqual(session.calls[0][1]['timeout'], (10, 60))
        with sqlite3.connect(self.out / 'proteins.db') as conn:
            self.assertEqual(conn.execute('SELECT count(*) FROM query_proteins').fetchone()[0], 3)
        self.assertEqual((self.out / 'all_sequences.fasta').read_text().count('>'), 2)
        self.assertEqual((self.out / 'seeds/murB_seeds.faa').read_text().count('>'), 1)
        self.assertEqual(json.loads((self.out / 'manifest.json').read_text())['status'], 'complete')
    def test_zero_matches(self):
        manifest = scraper.download([('empty', 'query')], self.out, True, Session(Response('', total=0)))
        self.assertEqual(manifest['unique_proteins'], 0)
        self.assertEqual((self.out / 'seeds/empty_seeds.faa').read_text(), '')
    def test_failure_retains_partial_and_no_final_output(self):
        session = Session(Response(next_url=scraper.API + '?cursor=abc', total=2), requests.Timeout('timeout'))
        with self.assertRaises(requests.Timeout):
            scraper.download([('gene', 'query')], self.out, session=session)
        self.assertFalse(self.out.exists())
        partial = next(self.root.glob('.results.partial-*'))
        self.assertEqual(json.loads((partial / 'manifest.json').read_text())['status'], 'incomplete')
        with sqlite3.connect(partial / 'proteins.db') as conn:
            self.assertEqual(conn.execute('SELECT count(*) FROM proteins').fetchone()[0], 1)
    def test_http_error(self):
        with self.assertRaises(requests.HTTPError):
            scraper.download([('gene', 'query')], self.out, session=Session(Response(error=requests.HTTPError('400'))))
        self.assertFalse(self.out.exists())
    def test_existing_output_untouched(self):
        self.out.mkdir()
        (self.out / 'keep').write_text('original')
        with self.assertRaises(ValueError):
            scraper.download([('gene', 'query')], self.out, session=Session())
        self.assertEqual((self.out / 'keep').read_text(), 'original')
    def test_malformed_fasta(self):
        for data in ['<html>error</html>', '>bad\nMPEPTIDE', FASTA.split('\n')[0]]:
            with self.subTest(data=data), self.assertRaises(ValueError):
                scraper.parse_fasta(data)
    def test_bad_labels_and_empty_queries(self):
        for queries in [[], [('a', '')], [('../escape', 'q')], [('gene', 'q'), ('GENE', 'q')]]:
            with self.subTest(queries=queries), self.assertRaises(ValueError):
                scraper.validate_queries(queries)
    def test_count_mismatch(self):
        with self.assertRaisesRegex(ValueError, 'count'):
            scraper.download([('gene', 'q')], self.out, session=Session(Response(total=2)))
    def test_release_change(self):
        with self.assertRaisesRegex(ValueError, 'release changed'):
            scraper.download([('gene', 'q')], self.out, session=Session(Response(next_url=scraper.API + '?cursor=1', total=2), Response(SECOND, total=2, release='2026_02')))
    def test_query_file(self):
        path = self.root / 'queries.tsv'
        path.write_text('gene\tquery\nmurA\tgene:murA AND reviewed:true\n')
        self.assertEqual(scraper.read_queries(path), [('murA', 'gene:murA AND reviewed:true')])
    def test_cli_failure_exit(self):
        with patch.object(scraper, 'read_queries', return_value=[('gene', 'q')]), patch.object(scraper, 'download', side_effect=requests.Timeout('timeout')):
            self.assertEqual(scraper.main(['--queries', 'anything.tsv']), 1)
    def test_retry_policy(self):
        with scraper.make_session() as session:
            retry = session.get_adapter(scraper.API).max_retries
            self.assertEqual(retry.total, 3)
            self.assertIn(429, retry.status_forcelist)
            self.assertTrue(retry.respect_retry_after_header)

if __name__ == '__main__':
    unittest.main()
