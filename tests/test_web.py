import pytest
import os

from comfy_BED.web import getLrgId, checkLrgExists

def test_getLrgId():
    # lrg ID and HGNC gene name
    assert getLrgId('LRG_1') == 'LRG_1'
    assert getLrgId('COL1A1') == 'LRG_1'

    # ensembl transcripts
    assert getLrgId('ENSG00000108821') == 'LRG_1'
    assert getLrgId('ENSP00000225964.5') == 'LRG_1'
    assert getLrgId('ENST00000225964.9') == 'LRG_1'

    # refseq transcripts
    assert getLrgId('NM_000088.3') == 'LRG_1'
    assert getLrgId('NG_007400.1') == 'LRG_1'
    assert getLrgId('NP_000079.2') == 'LRG_1'

    # Invalid inputs should throw an assertion error or value error
    with pytest.raises((ValueError, AssertionError)):
        getLrgId('invalid_input')


def test_CheckLrgExists():
    # LRG ID that exists
    assert checkLrgExists('LRG_1') == True

    # LRG ID that doesn't exist - should throw assertion error
    with pytest.raises(AssertionError):
        checkLrgExists('not_an_lrg')


def test_getLrgStatus():
    pass

def test_getLrgXml():
    pass

def test_getLrgWeb():
    pass