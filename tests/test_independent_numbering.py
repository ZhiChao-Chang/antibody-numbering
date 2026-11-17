from pathlib import Path
import json
import random
import socket
import sys
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from antibody_numbering.engine import AnchorNumberer,align,centered,loop_labels

ASSETS=ROOT/'data/reference/independent_numbering_public_20261109'


class IndependentNumberingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine=AnchorNumberer(ASSETS/'framework_templates.json')
        cls.public=json.loads((ROOT/'tests/fixtures/independent_public_6vja_H.json').read_text())
        cls.heavy=cls.public['sequence']

    def test_alignment_traceback_preserves_insertions(self):
        _,pairs=align('ABCDE','ABXCDE')
        self.assertIn((None,2),pairs)
        self.assertEqual([j for _,j in pairs if j is not None],list(range(6)))

    def test_unknown_truncated_terminal_is_retained(self):
        k=self.public['reference']['schemes']['kabat']['positions'];c=k.index('H92')
        sequence=self.heavy[:c+3]+'X';result=self.engine.number(sequence)
        self.assertEqual(result['schemes']['kabat']['positions'][-1],'H95')
        self.assertEqual(result['schemes']['imgt']['positions'][-1],'H117')
        self.assertEqual(result['schemes']['aho']['positions'][-1],'H109')

    def test_conserved_anchor_substitution_does_not_force_rejection(self):
        k=self.public['reference']['schemes']['kabat']['positions']
        for label,aa in [('H22','F'),('H36','Q'),('H92','A')]:
            i=k.index(label);seq=self.heavy[:i]+aa+self.heavy[i+1:];r=self.engine.number(seq)
            self.assertEqual(r['status'],'numbered')
            self.assertEqual(r['schemes']['kabat']['positions'],k)

    def test_heavy_h1_insertion_has_scheme_specific_site(self):
        k=loop_labels('H',1,12,'kabat');c=loop_labels('H',1,12,'chothia')
        self.assertIn('35B',k);self.assertIn('31B',c)
        self.assertEqual(len(set(k)),12);self.assertEqual(len(set(c)),12)

    def test_long_imgt_loop_preserves_decimal_insertion_order(self):
        self.assertEqual(centered(105,117,14),[str(x)for x in range(105,112)]+['112.1']+[str(x)for x in range(112,118)])

    def test_aho_legacy_partial_output_is_explicit(self):
        labels=loop_labels('H',2,19,'aho')
        self.assertEqual(labels[:6],['57','58','59','60','61','65']);self.assertEqual(labels[6:],['']*13)

    def test_long_light_loop_preserves_kappa_gap_positions(self):
        kappa=loop_labels('L',1,17,'aho',light_subtype='kappa')
        self.assertEqual(kappa,['24','25','26']+[str(x)for x in range(29,43)])
        self.assertNotIn('27',kappa);self.assertNotIn('28',kappa)

    def test_lambda_and_kappa_light_loops_have_distinct_gaps(self):
        kappa=loop_labels('L',1,13,'aho',light_subtype='kappa')
        lam=loop_labels('L',1,13,'aho',light_subtype='lambda')
        self.assertNotEqual(kappa,lam)
        self.assertIn('27',lam);self.assertNotIn('27',kappa)

    def test_n_terminal_truncation_retains_original_standard_labels(self):
        public=json.loads((ROOT/'tests/fixtures/independent_public_6vja_L.json').read_text())
        result=self.engine.number(public['sequence'][1:])
        for scheme,reference in public['reference']['schemes'].items():
            self.assertEqual(result['schemes'][scheme]['positions'],reference['positions'][1:])

    def test_c_terminal_framework_truncation_does_not_jump_to_endpoint(self):
        public=json.loads((ROOT/'tests/fixtures/independent_public_6vja_L.json').read_text())
        result=self.engine.number(public['sequence'][:-3])
        for scheme,reference in public['reference']['schemes'].items():
            self.assertEqual(result['schemes'][scheme]['positions'],reference['positions'][:-3])

    def test_fifteen_residue_kappa_loop_preserves_right_flank(self):
        labels=loop_labels('L',1,15,'aho',light_subtype='kappa')
        self.assertEqual(labels,['24','25','26']+[str(x)for x in range(29,36)]+[str(x)for x in range(38,43)])

    def test_inference_is_offline_and_has_no_query_reference_argument(self):
        with patch.object(socket,'socket',side_effect=AssertionError('network forbidden')):
            result=self.engine.number(self.heavy)
        self.assertEqual(result['status'],'numbered')
        self.assertNotIn('anarcii',sys.modules);self.assertNotIn('anarci',sys.modules)

    def test_framework_assets_mask_all_loop_answers(self):
        data=json.loads((ASSETS/'framework_templates.json').read_text())
        for t in data['templates'].values():
            framework={i for a,b in t['segments'].values()for i in range(a,b)}
            for i,aa in enumerate(t['sequence']):
                if i not in framework:
                    self.assertEqual(aa,'X')
                    self.assertTrue(all(labels[i]==''for labels in t['labels'].values()))

    def test_non_receptor_controls_are_not_assigned_positions(self):
        rng=random.Random(20261109)
        seqs=['A'*120,'ACDEFGHIKLMNPQRSTVWY'*6]+[''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=120))for _ in range(100)]
        for sequence in seqs:self.assertEqual(self.engine.number(sequence)['status'],'unclassified')



if __name__=='__main__':unittest.main()
