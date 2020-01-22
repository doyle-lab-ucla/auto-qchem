import os
import unittest


class TestOpenBabel(unittest.TestCase):

    def test_openbabel_datadir(self):
        ob_datadir = os.environ['BABEL_DATADIR']
        self.assertTrue(os.path.exists(ob_datadir))
