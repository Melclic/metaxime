import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys
import csv

sys.path.insert(0, '../..')

from metaxime import rpExtractSink

class TestRPextractSink(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        #load a rpSBML file
        #self.data = json.load(open(os.path.join('data', 'rpselenzyme', 'data.json'), 'r'))
        self.maxDiff = None

    #def _convertToCobra

    #def _reduce_model

    #def _removeDeadEnd

    def test_genSink(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpExtractSink.genSink(os.path.join('data', 'rpextractsink', 'gem.xml'), os.path.join(tmp_output_folder, 'tmp.csv'), remove_dead_end=True)
            result_dict = {}
            with open(os.path.join(tmp_output_folder, 'tmp.csv'), mode='r') as in_file:
                csv_reader = csv.reader(in_file)
                result_dict = {rows[0]:rows[1] for rows in csv_reader}
            self.assertIn('MNXM1860', result_dict)
            self.assertEqual(result_dict['MNXM1860'], 'InChI=1S/C10H15N4O15P3/c15-5-3(1-26-31(22,23)29-32(24,25)28-30(19,20)21)27-9(6(5)16)14-2-11-4-7(14)12-10(18)13-8(4)17/h2-3,5-6,9,15-16H,1H2,(H,22,23)(H,24,25)(H2,19,20,21)(H2,12,13,17,18)/p-4/t3-,5-,6-,9-/m1/s1')
            self.assertIn('MNXM7343', result_dict)
            self.assertEqual(result_dict['MNXM7343'], 'InChI=1S/C2H6O3S/c1-2-6(3,4)5/h2H2,1H3,(H,3,4,5)/p-1')


if __name__ == '__main__':
    unittest.main()
