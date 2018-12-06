import pytest
import os
import xml.etree.ElementTree as ET

from comfy_BED.comfy_BED.comfy_BED import checkValidTranscripts, calculateGenomicPositions, getLrgExons

'''
def test_get_args(self):
    args = get_args(['-t', 't1', '.xml'])
    assert os.path.isfile(args.path), 'The input is not a file.'
    assert args.path.endswith('.xml'), 'The input file is not an xml file.'
'''

def test_checkValidTranscripts():
    # Setup - make a root xml tree for each test set
    # Setup - need mock arg.transcripts
    # 5' 3' public
    tree_LRG_5 = ET.parse(os.path.abspath("test_data/LRG_5.xml"))
    root_LRG_5 = tree_LRG_5.getroot()

    # 3' 5' public
    tree_LRG_9 = ET.parse(os.path.abspath("test_data/LRG_9.xml"))
    root_LRG_9 = tree_LRG_9.getroot()

    # Pending LRG with gene in 5' -> 3' orientation
    tree_LRG_9 = ET.parse(os.path.abspath("test_data/LRG_9.xml"))
    root_LRG_9 = tree_LRG_9.getroot()

    common_transcript = "t1"
    common_transcript_list = "t1,t2"
    rarer_transcript = "t4"
    invalid_transcript = "t5"

    assert checkValidTranscripts(common_transcript, root_LRG_5) == True, "LRG_5 transcript 1 wrongly marked as invalid"


def test_getGenomeMapping():
    # Setup - make a root xml tree for each test set:
    # Public LRG with gene in 3' -> 5' orientation
    tree_LRG_5 = ET.parse(os.path.abspath("test_data/LRG_5.xml"))
    root_LRG_5 = tree_LRG_5.getroot()

    # Pending LRG with gene in 5' -> 3' orientation
    tree_LRG_9 = ET.parse(os.path.abspath("test_data/LRG_9.xml"))
    root_LRG_9 = tree_LRG_9.getroot()

    # Public LRG with gene in 5' -> 3' orientation
    tree_LRG_293 = ET.parse(os.path.abspath("test_data/LRG_293.xml"))
    root_LRG_293 = tree_LRG_293.getroot()

    # LRG_293 with mapping section removed - should throw error
    tree_LRG_293_mapping_removed = ET.parse(os.path.abspath("test_data/LRG_293_mapping_removed.xml"))
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

        
def test_getLrgExons():
    '''
    Test that getLrgExons takes transcript and lrg_id, and makes a correct dictionary
    (in which the key is exon label (string), and the value is a tuple of 'start' (int) and 'end' (int) positions).
    Try 2 LRG IDs and some different transcripts.
    '''

    #LRG_293 variables
    #LRG_293 has only t1. Strand = 1.
    tree_LRG_293 = ET.parse(os.path.abspath("test_data/LRG_293.xml"))
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
    tree_LRG_5 = ET.parse(os.path.abspath("test_data/LRG_5.xml"))
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

           
def test_calculateGenomicPositions():
    # Setup - define transcript dicts
    # LRG_5 - three transcripts
    dict_LRG_5_t1 = {'exon_13': (24287, 24362), 'exon_12': (23768, 23885), 'exon_11': (21749, 21899), 'exon_10': (19716, 19811), 
                     'exon_15': (25233, 25750), 'exon_14a': (24673, 24813), 'exon_9': (19421, 19548), 'exon_8': (17095, 17216), 
                     'exon_3': (12695, 12884), 'exon_2': (9610, 9762), 'exon_1': (5001, 5578), 'exon_7': (16868, 16920), 
                     'exon_6': (16448, 16537), 'exon_5': (14163, 14302), 'exon_4': (13102, 13233)}
    dict_LRG_5_t2 = {'exon_13': (24287, 24362), 'exon_12': (23768, 23885), 'exon_11': (21749, 21899), 'exon_10': (19716, 19811), 
                     'exon_15': (25233, 25750), 'exon_14b': (24673, 24832), 'exon_9': (19421, 19548), 'exon_8': (17095, 17216), 
                     'exon_3': (12695, 12884), 'exon_2': (9610, 9762), 'exon_1': (5001, 5578), 'exon_7': (16868, 16920), 
                     'exon_6': (16448, 16537), 'exon_5': (14163, 14302), 'exon_4': (13102, 13233)}
    dict_LRG_5_t3 = {'exon_13': (24287, 24362), 'exon_12': (23768, 23885), 'exon_11': (21749, 21899), 'exon_10': (19716, 19811), 
                     'exon_14c': (24673, 25750), 'exon_9': (19421, 19548), 'exon_8': (17095, 17216), 'exon_3': (12695, 12884), 
                     'exon_2': (9610, 9762), 'exon_1': (5001, 5578), 'exon_7': (16868, 16920), 'exon_6': (16448, 16537), 
                     'exon_5': (14163, 14302), 'exon_4': (13102, 13233)}

    # LRG_9 - one transcript
    dict_LRG_9_t1 = {'exon_3': (7021, 7165), 'exon_2': (6011, 6127), 'exon_1': (4978, 5113), 'exon_4': (12959, 13955)}

    # LRG_293 - one transcript
    dict_LRG_293_t1 = {'exon_22': (68838, 69036), 'exon_23': (69271, 69434), 'exon_20': (60477, 60621), 'exon_21': (66191, 66312), 
                       'exon_26': (86419, 86565), 'exon_27': (87683, 89193), 'exon_24': (69528, 69666), 'exon_25': (84210, 84454), 
                       'exon_13': (36348, 36417), 'exon_12': (34079, 34174), 'exon_11': (25786, 30717), 'exon_10': (21793, 22908), 
                       'exon_17': (52044, 52214), 'exon_16': (47263, 47450), 'exon_15': (45949, 46130), 'exon_14': (44382, 44809), 
                       'exon_19': (59923, 60078), 'exon_18': (52700, 53054), 'exon_9': (20440, 20551), 'exon_8': (18964, 19013), 
                       'exon_3': (8598, 8846), 'exon_2': (5943, 6048), 'exon_1': (5001, 5188), 'exon_7': (16020, 16134), 
                       'exon_6': (15763, 15803), 'exon_5': (15622, 15671), 'exon_4': (14597, 14705)}


    # Test calculateGenomicPositions produces the correct output - a list of tuples
    # input is transcript_dict, chrom, gen_start, gen_end, gen_strand
    assert calculateGenomicPositions(dict_LRG_5_t1, 'chr1', 43210006, 43237755, '-1') == [
        ('chr1', 43213394, 43213469, 'exon_13'), ('chr1', 43213871, 43213988, 'exon_12'), ('chr1', 43215857, 43216007, 'exon_11'), 
        ('chr1', 43217945, 43218040, 'exon_10'), ('chr1', 43212006, 43212523, 'exon_15'), ('chr1', 43212943, 43213083, 'exon_14a'), 
        ('chr1', 43218208, 43218335, 'exon_9'), ('chr1', 43220540, 43220661, 'exon_8'), ('chr1', 43224872, 43225061, 'exon_3'), 
        ('chr1', 43227994, 43228146, 'exon_2'), ('chr1', 43232178, 43232755, 'exon_1'), ('chr1', 43220836, 43220888, 'exon_7'), 
        ('chr1', 43221219, 43221308, 'exon_6'), ('chr1', 43223454, 43223593, 'exon_5'), ('chr1', 43224523, 43224654, 'exon_4')]

    assert calculateGenomicPositions(dict_LRG_5_t2, 'chr1', 43210006, 43237755, '-1') == [
        ('chr1', 43213394, 43213469, 'exon_13'), ('chr1', 43213871, 43213988, 'exon_12'), ('chr1', 43215857, 43216007, 'exon_11'), 
        ('chr1', 43217945, 43218040, 'exon_10'), ('chr1', 43212006, 43212523, 'exon_15'), ('chr1', 43212924, 43213083, 'exon_14b'), 
        ('chr1', 43218208, 43218335, 'exon_9'), ('chr1', 43220540, 43220661, 'exon_8'), ('chr1', 43224872, 43225061, 'exon_3'), 
        ('chr1', 43227994, 43228146, 'exon_2'), ('chr1', 43232178, 43232755, 'exon_1'), ('chr1', 43220836, 43220888, 'exon_7'), 
        ('chr1', 43221219, 43221308, 'exon_6'), ('chr1', 43223454, 43223593, 'exon_5'), ('chr1', 43224523, 43224654, 'exon_4')]

    assert calculateGenomicPositions(dict_LRG_5_t3, 'chr1', 43210006, 43237755, '-1') == [
        ('chr1', 43213394, 43213469, 'exon_13'), ('chr1', 43213871, 43213988, 'exon_12'), ('chr1', 43215857, 43216007, 'exon_11'), 
        ('chr1', 43217945, 43218040, 'exon_10'), ('chr1', 43212006, 43213083, 'exon_14c'), ('chr1', 43218208, 43218335, 'exon_9'), 
        ('chr1', 43220540, 43220661, 'exon_8'), ('chr1', 43224872, 43225061, 'exon_3'), ('chr1', 43227994, 43228146, 'exon_2'), 
        ('chr1', 43232178, 43232755, 'exon_1'), ('chr1', 43220836, 43220888, 'exon_7'), ('chr1', 43221219, 43221308, 'exon_6'), 
        ('chr1', 43223454, 43223593, 'exon_5'), ('chr1', 43224523, 43224654, 'exon_4')]

    assert calculateGenomicPositions(dict_LRG_9_t1, 'chr11', 111952571, 111992353, '1') == [
        ('chr11', 111959591, 111959735, 'exon_3'), ('chr11', 111958581, 111958697, 'exon_2'), ('chr11', 111957548, 111957683, 'exon_1'), 
        ('chr11', 111965529, 111966525, 'exon_4')]

    assert calculateGenomicPositions(dict_LRG_293_t1, 'chr13', 32884617, 32975809, '1') == [
        ('chr13', 32953454, 32953652, 'exon_22'), ('chr13', 32953887, 32954050, 'exon_23'), ('chr13', 32945093, 32945237, 'exon_20'), 
        ('chr13', 32950807, 32950928, 'exon_21'), ('chr13', 32971035, 32971181, 'exon_26'), ('chr13', 32972299, 32973809, 'exon_27'), 
        ('chr13', 32954144, 32954282, 'exon_24'), ('chr13',32968826, 32969070, 'exon_25'), ('chr13', 32920964, 32921033, 'exon_13'), 
        ('chr13', 32918695, 32918790, 'exon_12'), ('chr13', 32910402, 32915333, 'exon_11'), ('chr13', 32906409, 32907524, 'exon_10'), 
        ('chr13', 32936660, 32936830, 'exon_17'), ('chr13', 32931879, 32932066, 'exon_16'), ('chr13', 32930565, 32930746, 'exon_15'), 
        ('chr13', 32928998, 32929425, 'exon_14'), ('chr13', 32944539, 32944694, 'exon_19'), ('chr13', 32937316, 32937670, 'exon_18'), 
        ('chr13', 32905056, 32905167, 'exon_9'), ('chr13', 32903580, 32903629, 'exon_8'), ('chr13', 32893214, 32893462, 'exon_3'), 
        ('chr13', 32890559, 32890664, 'exon_2'), ('chr13', 32889617, 32889804, 'exon_1'), ('chr13', 32900636, 32900750, 'exon_7'), 
        ('chr13', 32900379, 32900419, 'exon_6'), ('chr13', 32900238, 32900287, 'exon_5'), ('chr13', 32899213, 32899321, 'exon_4')]


    # Invalid strand inputs - should throw ValueError
    with pytest.raises(ValueError):
        calculateGenomicPositions(dict_LRG_293_t1, 'chr13', 32884617, 32975809, 1)
    with pytest.raises(ValueError):
        calculateGenomicPositions(dict_LRG_293_t1, 'chr13', 32884617, 32975809, '')
    with pytest.raises(ValueError):
        calculateGenomicPositions(dict_LRG_293_t1, 'chr13', 32884617, 32975809, '2')

