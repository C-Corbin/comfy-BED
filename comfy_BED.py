import argparse
import textwrap
import os
import xml.etree.ElementTree as ET


# load arguments
def get_args():
    """
    Use argparse package to take arguments from the command line. 
    See descriptions for full detail of each argument.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent(
        '''
        summary:
        Takes an LRG xml file and returns a BED file of exon co-ordinates.
        '''
    ))

    # path to xml - test is a file, xml ending
    parser.add_argument(
        'input_LRG', action='store', 
        help='Filepath to input LRG xml file. REQUIRED.'
    )

    # transcript options
    parser.add_argument(
        '-t', '--transcripts', action='store', 
        help=textwrap.dedent(
        '''
        List of transcript to include
        '''
    ))

    return parser.parse_args()

def getLrgExons(transcript, lrg_id):
    '''
    Make a dictionary (key is string exon label, value is tuple of 'start' and 'end' 
    positions as integers)
    '''
    transcript_dict = {}
    for exon in transcript.iter('exon'):
        name = 'exon_{}'.format(exon.get('label'))
        #from coordinates
        for coordinate in exon.iter('coordinates'):
            if coordinate.get('coord_system') == lrg_id:
                start = int(coordinate.get('start'))
                end = int(coordinate.get('end'))
        transcript_dict[name] = (start, end)
    return(transcript_dict)


def main():
    args = get_args()
    assert os.path.isfile(args.input_LRG), 'The input is not a file.'
    assert args.input_LRG.endswith('.xml'), 'The input file is not an xml file.'

    # create xml element tree object
    # test that file is an lrg (root.tag)
    tree = ET.parse(os.path.abspath(args.input_LRG))
    root = tree.getroot()
    assert root.tag.upper() == "LRG", 'The input file is not an LRG file'
    for levels in root.iter('fixed_annotation'):
        lrg_id = levels.find('id').text
    # extract the exon boundries - lrg numbering - make into python dict
    for transcript in root.iter('transcript'):
        if str(transcript.get('name')) in args.transcripts:
            transcript_dict = getLrgExons(transcript, lrg_id)
            print(transcript_dict)

    # extract chr, start, end, strand from mapping region of xml
    # calculate genomic coordinates, depending on stand orientation
        # get strand: 1 or -1
        # add 5000 for 5' or 2000 for 3'
        # add/subtract overall genomic coords to exon numbers (depending on strand)
    # output in tab delimted text file - chr, start, end, exon_no (.bed)


if __name__ == '__main__':
    main()
