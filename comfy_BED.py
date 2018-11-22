import argparse
import textwrap
import os


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



# create xml element tree object
  # test that file is an lrg (root.tag)

# extract transcripts - for each transcript:
  # extract the exon boundries - lrg numbering - make into python dict
  # extract chr, start, end, strand from mapping region of xml
  # calculate genomic coordinates, depending on stand orientation
    # get strand: 1 or -1
    # add 5000 for 5' or 2000 for 3'
    # add/subtract overall genomic coords to exon numbers (depending on strand)
  # output in tab delimted text file - chr, start, end, exon_no (.bed)

def main():
    args = get_args()
    assert os.path.isfile(args.input_LRG), 'The input is not a file.'
    assert args.input_LRG.endswith('.xml'), 'The input file is not an xml file.'

if __name__ == '__main__':
    main()
