import pytest
import os
import hashlib

from comfy_BED.comfy_BED_web import getLrgId, checkLrgExists, checkCurrentLrgStatus, getLrgXml


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


def test_checkLrgExists():
    # LRG ID that exists
    assert checkLrgExists('LRG_1') == True

    # LRG ID that doesn't exist - should throw assertion error
    with pytest.raises(AssertionError):
        checkLrgExists('invalid_input')


def test_checkCurrentLrgStatus():
    # Public and pending LRGs - returns tuple of status and status error, check that status matches
    assert checkCurrentLrgStatus('LRG_5')[0] == 'public'
    assert checkCurrentLrgStatus('LRG_9')[0] == 'pending'

    # LRG that doesn't exist (throws index error as it trys to parse json output that wasnt retrieved)
    with pytest.raises(IndexError):
        checkCurrentLrgStatus('invalid_input')


def test_getLrgXml():
    '''
    NOTE: This test works by comparing the md5 checksum of a locally 
    saved xml, downloaded through the LRG website, to the md5 checksum
    of the output from the web API.

    Since the LRGs are regularly updated, this test might fail. If it 
    does, download the latest LRG_1.xml file from the LRG website, 
    replace tests/test_data/LRG_1.xml and re-run the tests.

    I've used a new LRG from the different tests in case updating the 
    file affects the other tests, so we should keep LRG_1 for this test only
    '''
    # open local file and calculate md5 checksum
    with open('tests/test_data/LRG_1.xml') as LRG_1:
        LRG_1_file = LRG_1.read()
        LRG_1_file_md5 = hashlib.md5(LRG_1_file.encode('utf-8')).hexdigest()

    # get xml from web api and calculate md5 checksum
    LRG_1_web = getLrgXml('LRG_1', 'public')
    LRG_1_web_md5 = hashlib.md5(LRG_1_web.encode('utf-8')).hexdigest()

    # check that the checksums are the same
    assert LRG_1_file_md5 == LRG_1_web_md5
