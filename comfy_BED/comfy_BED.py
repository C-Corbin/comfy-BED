import argparse
import textwrap
import os
import csv
import xml.etree.ElementTree as ET
import datetime
import six
import logging

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


def setUpLogs(args, now):
    '''
    Makes a log file to help with spotting errors in LRG-to-BED conversion
    Log file name ends with .log and will contain the current date, LRG ID, and 'comfy_BED'
    Log file created by day and appends to day
    '''
    log_filename = now.strftime("%Y-%m-%d") + "_comfy_BED" + ".log"
    logging.basicConfig(filename=log_filename, level=logging.DEBUG)


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
    logging.info("Fetched start and end coordinates of LRG exons")
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
    logging.info("Obtained the LRG and genomic coordinates of the start and end of the selected LRG gene")
    return(chr, start, end, strand)


def calculateGenomicPositions(transcript_dict, chrom, gen_start, gen_end, strand):
    '''
    Convert LRG exon boundry positions into genomic positions.

    Input -
    transcript_dict: Dictionary from getLrgExons function, key is exon
        label, value is tuple of LRG start and LRG end.
    chrom (string):  The chromosome that the LRG is in, parsed from
        the LRG using getGenomeMapping.
    gen_start (int): The start position of the LRG within the genome,
        parsed from the LRG using getGenomeMapping.
    gen_end (int):   The end position of the LRG within the genome,
        parsed from the LRG using getGenomeMapping.
    strand (string): Whether the LRG is orientated in the genome from
        5' -> 3' (strand is '1'), or from 3' -> 5' (strand is '-1').
        Any value other than 1 or -1 will throw an error. Strand is
        parsed from the LRG using getGenomeMapping.

    Output -
    list_of_exons: List of tuples, one per exon in the transcript_dict.
        Each tuple contains the chromosome, genome start coordinate,
        genome end coordinate and exon label.
    '''
    # open empty list to add each exon tuple to
    list_of_exons = []

    # loop through exon dictionary, extract exon name, lrg start and end
    # use six library for iterating as it has support for both python 2 and 3
    for transcript in six.iteritems(transcript_dict):
        exon_label = str(transcript[0])
        lrg_start = int(transcript[1][0])
        lrg_end = int(transcript[1][1])

        # if strand is 5' -> 3'
        if strand == '1':
            gen_exon_start = gen_start + lrg_start - 1
            gen_exon_end = gen_start + lrg_end - 1
            list_of_exons.append((chrom, gen_exon_start, gen_exon_end, exon_label))

        # if strand is 3' -> 5'
        elif strand == '-1':
            exon_length = lrg_end - lrg_start
            gen_exon_end = gen_end - lrg_end + 1
            gen_exon_start = gen_exon_end + exon_length
            # start and end are switch round because bed files should have the smallest value first
            list_of_exons.append((chrom, gen_exon_end, gen_exon_start, exon_label))
        else:
            # raise a value error if strand is anything other than 1 or -1
            raise ValueError('Cannot determine strand')
    logging.info("Converted LRG start-and-end coordinates, to genomic coordinates, for the user-selected genome build and transcript")
    return list_of_exons


def writeToFile(data_list, file_name, now):
    '''
    Take a list of tuples and look through and write as a tab seperated file
    '''
    with open(file_name, 'wb') as out:
        out.writelines('#BED file generated at: ' + now.strftime("%Y-%m-%d %H:%M") + '\n')
        writer = csv.writer(out, delimiter='\t')
        for row in data_list:
            writer.writerow(row)
    logging.info("Wrote exon start-and-end coordinates, for the user-selected genome build and transcript, to BED file")
    logging.info("The BED file is named: " + file_name)


def main():
    args = getArgs()
    now = datetime.datetime.now()
    # set up logs
    setUpLogs(args, now)
    logging.info("comfy_BED started running at: " + str(now))

    #quick file checks
    assert os.path.isfile(args.input_LRG), 'The input is not a file.'
    if os.path.isfile(args.input_LRG) == False:
        logging.error("The input file is not a file")
    assert args.input_LRG.endswith('.xml'), 'The input file is not an xml file.'
    if args.input_LRG.endswith('.xml') == False:
        logging.error("The input file is not an .xml file")

    # create xml element tree object
    # test that file is an lrg (root.tag)
    tree = ET.parse(os.path.abspath(args.input_LRG))
    root = tree.getroot()
    assert root.tag.upper() == "LRG", 'The input file is not an LRG file'
    if root.tag.upper() != "LRG":
        logging.error("The input file is not an LRG file")

    # get lrg id
    for levels in root.iter('fixed_annotation'):
        lrg_id = levels.find('id').text
        logging.info("LRG_ID: " + lrg_id)

    # extract chr, start, end, strand from mapping region of xml
    logging.info("Genome build: " + str(args.genome_build))
    chrom, genome_start, genome_end, genome_strand = getGenomeMapping(root, args.genome_build)
    logging.info("Chromosome: " + chrom)
    logging.info("Strand: " + genome_strand)
    logging.info("Start position of gene on " + args.genome_build + ": " + str(genome_start))
    logging.info("End position of gene on " + args.genome_build + ": " + str(genome_end))    # extract the exon boundries - lrg numbering - make into python dict
    # calculate genomic coordinates, depending on strand orientation
    for transcript in root.iter('transcript'):
        transcript_name = str(transcript.get('name'))
        if transcript_name in args.transcripts:
            logging.info("Started BED production for transcript: " + transcript_name)
            transcript_dict = getLrgExons(transcript, lrg_id)
            exon_genomic_positions = calculateGenomicPositions(transcript_dict, chrom, genome_start, genome_end, genome_strand)

            # output in tab delimted text file
            #TODO add header, option to change filename, sorting
            file_name = '{}_{}.bed'.format(lrg_id, transcript_name)
            writeToFile(exon_genomic_positions, file_name, now)
            logging.info("Completed BED production for transcript: " + transcript_name)
    logging.info("comfy_BED run complete")

if __name__ == '__main__':
    main()
