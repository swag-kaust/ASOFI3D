import os
import shutil
import subprocess
import unittest


class Test01(unittest.TestCase):
    def test_01(self):
        SRC_DIR = 'src'
        OLD_MODEL = 'src/hh_elastic.c'
        BAK_MODEL = 'src/hh_elastic.c.bak'
        TEST_MODEL = 'tests/fixtures/test_01/hh_elastic.c'
        # Preserve old model.
        os.rename(OLD_MODEL, BAK_MODEL)
        # Copy test model.
        shutil.copy(TEST_MODEL, SRC_DIR)
        # Compile code.
        subprocess.call('make sofi3D', shell=True)
        # Run code.
        subprocess.call('./run_sofi3D.sh')
        # Read the files.
        # Compare with the old output.


        # Teardown
        os.rename(BAK_MODEL, OLD_MODEL)
