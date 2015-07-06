#!/usr/bin/python
from nose.tools import assert_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal
from pynrrd import *

import os.path as path

TEST_PATH = 'testdata'
def test_nhdr_open():
    r = NrrdReader()
    hdr, b = r.load(path.join(TEST_PATH, 'seed.nhdr'), get_raw=False)

    assert hdr._data
    assert_equal(hdr['type'], 'short')
    assert_equal(hdr['dimension'], '3')
    assert_equal(hdr['space'], 'left-posterior-superior')
    assert_array_equal(hdr['sizes'], [256,256,55])
    assert_array_equal(hdr['space directions'], [[0.9375000000000002, 0.0, 0.0], [0.0, -0.9375000000000002, 0.0], [0.0, 0.0, 3.000000000000001]])

    assert_array_equal(hdr['sizes'], b.shape)

def test_nrrd_open():
    r = NrrdReader()
    hdr, b = r.load(path.join(TEST_PATH, 'seed.nrrd'), get_raw=False)

    assert hdr._data
    assert_equal(hdr['type'], 'short')
    assert_equal(hdr['dimension'], '3')
    assert_equal(hdr['space'], 'left-posterior-superior')
    assert_array_equal(hdr['sizes'], [256,256,55])
    assert_array_almost_equal(hdr['space directions'], [[0.9375000000000002, 0.0, 0.0], [0.0, -0.9375000000000002, 0.0], [0.0, 0.0, 3.000000000000001]])

    assert_array_equal(hdr['sizes'], b.shape)

def test_dwi():
    r = NrrdReader()
    hdr, b = r.load(path.join(TEST_PATH, 'dwi.nhdr'), get_raw=False)    
    assert hdr.isDTMR()
    assert_equal(hdr.getBval(), 1000.)
    print hdr.getBvals()

if __name__ == "__main__":
    import nose
    nose.runmodule()