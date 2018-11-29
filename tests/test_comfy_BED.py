import pytest
import os
import xml.etree.ElementTree as ET

from comfy_BED.comfy_BED import getGenomeMapping

'''
def test_get_args(self):
    args = get_args(['-t', 't1', '.xml'])
    assert os.path.isfile(args.path), 'The input is not a file.'
    assert args.path.endswith('.xml'), 'The input file is not an xml file.'
'''


def test_getGenomeMapping():
    # Setup - make a root xml tree for each test set:
    # Public LRG with gene in 3' -> 5' orientation
    tree_LRG_5 = ET.parse(os.path.abspath("tests/test_data/LRG_5.xml"))
    root_LRG_5 = tree_LRG_5.getroot()

    # Pending LRG with gene in 5' -> 3' orientation
    tree_LRG_9 = ET.parse(os.path.abspath("tests/test_data/LRG_9.xml"))
    root_LRG_9 = tree_LRG_9.getroot()

    # Public LRG with gene in 5' -> 3' orientation
    tree_LRG_293 = ET.parse(os.path.abspath("tests/test_data/LRG_293.xml"))
    root_LRG_293 = tree_LRG_293.getroot()

    # LRG_293 with mapping section removed - should throw error
    tree_LRG_293_mapping_removed = ET.parse(os.path.abspath("tests/test_data/LRG_293_mapping_removed.xml"))
    root_LRG_293_mapping_removed = tree_LRG_293_mapping_removed.getroot()

    # Test that genome mapping is parsed correctly for each test set, with both genome builds
    assert getGenomeMapping(root_LRG_5, 'GRCh37') == ('chr1', 43210006, 43237755, '-1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_5, 'GRCh38') == ('chr1', 42744335, 42772084, '-1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_9, 'GRCh37') == ('chr11', 111952571, 111992353, '1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_9, 'GRCh38') == ('chr11', 112081847, 112121630, '1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_293, 'GRCh37') == ('chr13', 32884617, 32975809, '1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_293, 'GRCh38') == ('chr13', 32310480, 32401672, '1'), "Failed to obtain mapping information"
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_mapping_removed, 'GRCh37')
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_mapping_removed, 'GRCh38')
