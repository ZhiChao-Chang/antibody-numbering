"""Public failure-mode regressions for the private development candidate.

Expected mappings come from frozen receipts. The engine receives sequences only;
these fixtures are never installed as runtime framework assets.
"""
from pathlib import Path
import json
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from antibody_numbering.engine import AnchorNumberer,pack_residues,packed_substitution_score

class ExpandedRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine=AnchorNumberer(ROOT/'data/reference/independent_numbering_public_20261109/framework_templates.json')
        cls.cases=json.loads((ROOT/'tests/fixtures/expanded_public_regressions.json').read_text())

    def assert_category(self, category):
        for case in self.cases:
            if case['category']!=category:continue
            with self.subTest(case=case['id']):
                actual=self.engine.number(case['sequence'])
                if category=='negative':
                    self.assertEqual(actual['status'],'unclassified')
                    continue
                for scheme,expected in case['positions'].items():
                    self.assertEqual(actual['schemes'][scheme]['positions'],expected,scheme)

    def test_internal_cdr3_motif_does_not_hide_later_fr4(self):
        self.assert_category('internal_fr4_like_cdr3')

    def test_framework_indels_keep_downstream_anchors(self):
        self.assert_category('framework_indel')

    def test_aromatic_c2_variants_keep_variable_domain(self):
        self.assert_category('aromatic_c2')

    def test_short_heavy_h1_keeps_fr2_boundary(self):
        self.assert_category('short_h1')

    def test_ig_like_and_tcr_controls_remain_unclassified(self):
        self.assert_category('negative')

    def test_packed_scoring_matches_explicit_residue_scoring(self):
        import random
        alphabet='ACDEFGHIKLMNPQRSTVWYX';rng=random.Random(20261117)
        pairs=[(a,b)for a in alphabet for b in alphabet]
        pairs+=[(''.join(rng.choices(alphabet,k=n)),''.join(rng.choices(alphabet,k=n)))for n in range(41)for _ in range(20)]
        for reference,query in pairs:
            expected=sum(0 if 'X'in(a,b)else 3 if a==b else -1 for a,b in zip(reference,query))
            self.assertEqual(packed_substitution_score(pack_residues(reference),pack_residues(query)),expected)

    def test_truncated_domain_ending_at_cysteine_keeps_framework_labels(self):
        public=json.loads((ROOT/'tests/fixtures/independent_public_6vja_L.json').read_text())
        end=public['reference']['schemes']['kabat']['positions'].index('L88')+1
        result=self.engine.number(public['sequence'][:end])
        self.assertEqual(result['status'],'numbered')
        for scheme,reference in public['reference']['schemes'].items():
            self.assertEqual(result['schemes'][scheme]['positions'],reference['positions'][:end])

    def test_cdr_substitutions_do_not_require_a_query_answer_cache(self):
        # A single conservative substitution at a central CDR3 position creates
        # a new sequence while preserving every framework and loop length.
        case=next(c for c in self.cases if c['category']=='internal_fr4_like_cdr3')
        k=case['positions']['kabat'];i=k.index('H100')
        seq=case['sequence'];replacement='S'if seq[i]!='S'else'T'
        result=self.engine.number(seq[:i]+replacement+seq[i+1:])
        for scheme,expected in case['positions'].items():
            self.assertEqual(result['schemes'][scheme]['positions'],expected,scheme)

if __name__=='__main__':unittest.main()
