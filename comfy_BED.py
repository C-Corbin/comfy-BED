import argparse
import textwrap
import os
import xml.etree.ElementTree as ET


# load arguments
def getArgs():
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

    # genome build options
    parser.add_argument(
        '-g', '--genome_build', action='store', 
        choices=['GRCh37', 'GRCh38'], default='GRCh37',
        help=textwrap.dedent(
        '''
        Select genome build from 'GRCh37' or 'GRCh38'. Defaults to GRCh37.
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


def getGenomeMapping(root, genome_build):
    '''
    Parses chromosome number, strand direction and LRG start and end 
    postions on the genome build (default GCRh37). The start and end
    positions will be used to convert between LRG numbering and 
    genome numbering.
    '''
    for mapping in root.iter('mapping'):
        if str(mapping.get('coord_system')).startswith(genome_build):
            chr = 'chr{}'.format(mapping.get('other_name'))

            for m_span in mapping.iter('mapping_span'):
                start = int(m_span.get('other_start'))
                end = int(m_span.get('other_end'))
                strand = str(m_span.get('strand'))
    return(chr, start, end, strand)


def calculateGenomicPositions(chr, transcript_dict, genome_start, genome_strand):
    '''
    add the genome_start number to the LRG exon boundaries, to convert them
    to genomic exon boundaries
    '''
    list_of_exons = []
    if genome_strand == '1':
        for item in transcript_dict.iteritems():
            exon_label = str(item[0])
            lrg_start = int(item[1][0])
            lrg_end = int(item[1][1])
            gen_exon_start = genome_start + lrg_start
            gen_exon_end = genome_start + lrg_end
            list_of_exons.append((chr, gen_exon_start, gen_exon_end, exon_label))
    #else
    return list_of_exons


def main():
    args = getArgs()
    assert os.path.isfile(args.input_LRG), 'The input is not a file.'
    assert args.input_LRG.endswith('.xml'), 'The input file is not an xml file.'

    # create xml element tree object
    # test that file is an lrg (root.tag)
    tree = ET.parse(os.path.abspath(args.input_LRG))
    root = tree.getroot()
    assert root.tag.upper() == "LRG", 'The input file is not an LRG file'

    # get lrg id
    for levels in root.iter('fixed_annotation'):
        lrg_id = levels.find('id').text

    # extract chr, start, end, strand from mapping region of xml
    chr, genome_start, genome_end, genome_strand = getGenomeMapping(root, args.genome_build)

    # extract the exon boundries - lrg numbering - make into python dict
    for transcript in root.iter('transcript'):
        if str(transcript.get('name')) in args.transcripts:
            transcript_dict = getLrgExons(transcript, lrg_id)
            exon_genomic_positions = calculateGenomicPositions(chr, transcript_dict, genome_start, genome_strand)



    # calculate genomic coordinates, depending on stand orientation

        # get strand: 1 or -1
        # add 5000 for 5' or 2000 for 3'
        # add/subtract overall genomic coords to exon numbers (depending on strand)
    # output in tab delimted text file - chr, start, end, exon_no (.bed)


if __name__ == '__main__':
    main()
