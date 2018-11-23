import pytest
import os

from comfy_BED import *

'''
def test_get_args(self):
    args = get_args(['-t', 't1', '.xml'])
    assert os.path.isfile(args.path), 'The input is not a file.'
    assert args.path.endswith('.xml'), 'The input file is not an xml file.'
'''

def test_getGenomeMapping():
    #make the root for test sets
    #use one 5->3 and one 3->5 set
    #two genome_builds
    tree_LRG_293 = ET.parse(os.path.abspath("LRG_293.xml"))
    root_LRG_293 = tree_LRG_293.getroot()
    tree_LRG_293_bad = ET.parse(os.path.abspath("LRG_293_bad.xml"))
    root_LRG_293_bad = tree_LRG_293_bad.getroot()
    tree_LRG_5 = ET.parse(os.path.abspath("LRG_5.xml"))
    root_LRG_5 = tree_LRG_5.getroot()

    assert getGenomeMapping(root_LRG_293, 'GRCh37') == ('chr13', 32884617, 32975809, '1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_293, 'GRCh38') == ('chr13', 32310480, 32401672, '1'), "Failed to obtain mapping information"
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_bad, 'GRCh37')
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_bad, 'GRCh38')
    #assert getGenomeMapping(root_LRG_5, 'GRCh37') ==
    #assert getGenomeMapping(root_LRG_5, 'GRCh38') ==
    #it outputs these: chr, genome_start, genome_end, genome_strand


