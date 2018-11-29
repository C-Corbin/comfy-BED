import pytest
import xml.etree.ElementTree as ET
import os
from comfy_BED import *

'''
def test_get_args(self):
    args = get_args(['-t', 't1', '.xml'])
    assert os.path.isfile(args.path), 'The input is not a file.'
    assert args.path.endswith('.xml'), 'The input file is not an xml file.'
'''

def test_getGenomeMapping():
    '''
    make roots for two test LRGs, one which is 5->3 and one which is 3->5
    test sets on two genome_builds
    '''
    tree_LRG_293 = ET.parse(os.path.abspath("LRG_293.xml"))
    root_LRG_293 = tree_LRG_293.getroot()
    tree_LRG_293_bad = ET.parse(os.path.abspath("LRG_293_bad.xml"))
    root_LRG_293_bad = tree_LRG_293_bad.getroot()
    #tree_LRG_5 = ET.parse(os.path.abspath("LRG_5.xml"))
    #root_LRG_5 = tree_LRG_5.getroot()

    assert getGenomeMapping(root_LRG_293, 'GRCh37') == ('chr13', 32884617, 32975809, '1'), "Failed to obtain mapping information"
    assert getGenomeMapping(root_LRG_293, 'GRCh38') == ('chr13', 32310480, 32401672, '1'), "Failed to obtain mapping information"
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_bad, 'GRCh37')
    with pytest.raises(UnboundLocalError):
        getGenomeMapping(root_LRG_293_bad, 'GRCh38')
    #assert getGenomeMapping(root_LRG_5, 'GRCh37') ==
    #assert getGenomeMapping(root_LRG_5, 'GRCh38') ==
    #it outputs these: chr, genome_start, genome_end, genome_strand

def test_getLrgExons():
    '''
    Test that getLrgExons takes transcript and lrg_id, and makes a correct dictionary
    (in which the key is exon label (string), and the value is a tuple of 'start' (int) and 'end' (int) positions).
    Try 2 LRG IDs and some different transcripts.
    '''

    #LRG_293 variables
    #LRG_293 has only t1. Strand = 1.
    tree_LRG_293 = ET.parse(os.path.abspath("LRG_293.xml"))
    root_LRG_293 = tree_LRG_293.getroot()
    LRG_ID_293 = "LRG_293"
    LRG_293_t1 = "t1" #Valid transcript for LRG_293
    LRG_293_invalid_transcript = "t9" #This transcript doesn't exist for LRG_293

    #LRG_293 tests
    for transcript_LRG_293 in root_LRG_293.iter('transcript'):
        name_transcript_LRG_293 = str(transcript_LRG_293.get('name'))
        if name_transcript_LRG_293 == LRG_293_t1: #if this is the valid transcript t1
            assert getLrgExons(transcript_LRG_293, LRG_ID_293) == {"exon_1": (5001, 5188), "exon_2": (5943, 6048), "exon_3": (8598, 8846), "exon_4": (14597, 14705), "exon_5": (15622, 15671), "exon_6": (15763, 15803), "exon_7": (16020, 16134), "exon_8": (18964, 19013), "exon_9": (20440, 20551), "exon_10": (21793, 22908), "exon_11": (25786, 30717), "exon_12": (34079, 34174), "exon_13": (36348, 36417), "exon_14": (44382, 44809), "exon_15": (45949, 46130), "exon_16": (47263, 47450), "exon_17": (52044, 52214), "exon_18": (52700, 53054), "exon_19": (59923, 60078), "exon_20": (60477, 60621), "exon_21": (66191, 66312), "exon_22": (68838, 69036), "exon_23": (69271, 69434), "exon_24": (69528, 69666), "exon_25": (84210, 84454), "exon_26": (86419, 86565), "exon_27": (87683, 89193) }, "Failed to obtain correct LRG exon dictionary from XML (LRG_293, strand=1, t1)"
        with pytest.raises(UnboundLocalError):
            assert getLrgExons(transcript_LRG_293, LRG_293_invalid_transcript)

    #LRG_5 variables
    #LRG_5 has t1, t2 and t3. t3 was added after LRG_5 was made public. Strand = -1.
    tree_LRG_5 = ET.parse(os.path.abspath("LRG_5.xml"))
    root_LRG_5 = tree_LRG_5.getroot()
    LRG_ID_5 = "LRG_5"
    LRG_5_t1 = "t1" #Valid transcript
    LRG_5_t3 = "t3" #Another valid transcript
    LRG_5_invalid_transcript = "t10" #This transcript doesn't exist for LRG_5

    #LRG_5 tests
    for transcript_LRG_5 in root_LRG_5.iter('transcript'):
        name_transcript_LRG_5 = str(transcript_LRG_5.get('name'))
        if name_transcript_LRG_5 == LRG_5_t1: #if this is valid transcript t1
            assert getLrgExons(transcript_LRG_5, LRG_ID_5) == {"exon_1": (5001, 5578), "exon_2": (9610, 9762), "exon_3": (12695, 12884), "exon_4": (13102, 13233), "exon_5": (14163, 14302), "exon_6": (16448, 16537), "exon_7": (16868, 16920), "exon_8": (17095, 17216), "exon_9": (19421, 19548), "exon_10": (19716, 19811), "exon_11": (21749, 21899), "exon_12": (23768, 23885), "exon_13": (24287, 24362), "exon_14a": (24673, 24813), "exon_15": (25233, 25750)}, "Failed to obtain correct LRG exon dictionary from XML (LRG_5, strand=-1, t1)"
        if name_transcript_LRG_5 == LRG_5_t3:
            assert getLrgExons(transcript_LRG_5, LRG_ID_5) == {"exon_1": (5001, 5578), "exon_2": (9610, 9762), "exon_3": (12695, 12884), "exon_4": (13102, 13233), "exon_5": (14163, 14302), "exon_6": (16448, 16537), "exon_7": (16868, 16920), "exon_8": (17095, 17216), "exon_9": (19421, 19548), "exon_10": (19716, 19811), "exon_11": (21749, 21899), "exon_12": (23768, 23885), "exon_13": (24287, 24362), "exon_14c": (24673, 25750)}, "Failed to obtain correct LRG exon dictionary from XML (LRG_5, strand=-1, t3)"
        with pytest.raises(UnboundLocalError):
            assert getLrgExons(transcript_LRG_5, LRG_5_invalid_transcript)

test_getLrgExons()