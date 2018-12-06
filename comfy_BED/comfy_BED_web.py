import requests
import logging
import xml.etree.ElementTree as ET


def getLrgId(input_text):
    '''
    Get the ID number of an LRG

    Input -
    input_text: String. Name to extract LRG ID from, can be an LRG 
      number, a HGNC gene symbol or a transcript name

    Output -
    lrg_id: String. The LRG ID in the format LRG_<number>. If the 
      LRG ID can't be calculated from the input, an error will be thrown.
    '''
    logging.info('Web query input: {}'.format(input_text))

    # if input is an lrg number, save the variable
    if input_text.startswith('LRG_'):
        lrg_id = input_text

    # if input isn't an lrg number, try to query by name to find lrg number
    else:
        logging.info('Querying webservices to get LRG ID and check that it is valid')

        # try to query by HGNC name
        try:
            name_query_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=name:{}'.format(input_text)
            name_query_response = requests.get(name_query_url)

            if name_query_response.status_code != 200:
                logging.error('Could not query the API, check your connection and try again.')
            assert name_query_response.status_code == 200, 'Could not query the API, check your connection and try again.'

            # loop through the response xml
            root = ET.fromstring(name_query_response.text)
            for child in root.iter('hitCount'):
                # check that there is exactly 1 entry returned
                if child.text == '1':
                    pass
                elif child.text == '0':
                    logging.error('There were no hits for {}, check the input'.format(input_text))
                else:
                    logging.error('Expected one hit but there were multiple, check the input')
                assert child.text == '1'
            # if so, extract the lrg id and save as a variable
            for child in root.iter('entry'):
                lrg_id = child.get('id')
                logging.info('Found LRG ID for {}: {}'.format(input_text, lrg_id))


        # try to query by other references
        except AssertionError:
            name_query_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query={}'.format(input_text)
            name_query_response = requests.get(name_query_url)
                
            if name_query_response.status_code != 200:
                logging.error('Could not query the API, check your connection and try again.')
            assert name_query_response.status_code == 200, 'Could not query the API, check your connection and try again.'

            # loop through the response xml
            root = ET.fromstring(name_query_response.text)
            for child in root.iter('hitCount'):
                # check that there is exactly 1 entry returned
                if child.text == '1':
                    pass
                elif child.text == '0':
                    logging.error('There were no hits for {}, check the input'.format(input_text))
                else:
                    logging.error('Expected one hit but there were multiple, check the input')
                assert child.text == '1'
            # if so, extract the lrg id and save as a variable
            for child in root.iter('entry'):
                lrg_id = child.get('id')
                logging.info('Found LRG ID for {}: {}'.format(input_text, lrg_id))

        # throw error if both queries fail
        except:
            logging.error('Cannot find the LRG file from the given input.')
            raise ValueError('Cannot find the LRG file from the given input.')

    return(lrg_id)


def checkLrgExists(lrg_id):
    '''
    Check that an LRG exists

    Input -
    lrg_id: String. An LRG ID to check whether or not it exists. 
      Must be in the format LRG_<number>

    Output -
    Boolean. If the function runs without hitting an assertion error, 
    the function will return as true, since the LRG ID can be found 
    through the API and therefore it exists.
    '''
    # query api, returns xml that says whether lrg exists or not
    lrg_query_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query={}'.format(lrg_id)
    lrg_query_response = requests.get(lrg_query_url)
    assert lrg_query_response.status_code == 200, 'Could not query the API, check your connection and try again.'

    # parse the section that says if lrg exists or not, throw assertion error if it doesn't
    root = ET.fromstring(lrg_query_response.text)
    for child in root.iter('hitCount'):
        assert child.text == '1', 'LRG does not exist'
    
    # if no assertion error is thrown, return true
    return True


def checkCurrentLrgStatus(lrg_id):
    '''
    Checks the CURRENT status of the user-provided LRG ID on LRG website, returns information to the BED header and log file
    'Public' LRGs have a 'fixed annotation' section which has been fully finalised
    'Pending' LRGs do NOT have a finalised 'fixed annotation' section
    WARNING: User-provided XMLs contain no information as to whether they were public or private at time of download
    The end user should always download their XMLs very shortly before use
    '''
    #get data from webservice
    url_p1 = "https://www.ebi.ac.uk/ebisearch/ws/rest/lrg/entry/"
    url_p3 = "?fields=status&format=json"
    url_full = url_p1 + str(lrg_id) + url_p3
    logging.info("Checking status with webservice: " + str(url_full))
    data_return = requests.get(url_full)
    parsed_data_return = data_return.json()

    # parse the returned data, return status and log message
    lrg_status_return = parsed_data_return['entries'][0]['fields']['status'][0]
    if lrg_status_return == "public":
        lrg_status_message = "The LRG is currently marked 'public' on the LRG website: note that the user-provided file could have been downloaded before the LRG going public"
    else:
        if lrg_status_return == "pending":
            lrg_status_message = "The LRG is currently marked 'pending' on the LRG website: the fixed annotation is not yet finalised, so it should be interpreted with caution"
    assert (lrg_status_return == "public") or (lrg_status_return == "pending"), "The LRG status could not be resolved as public or private"
    if (lrg_status_return != "public") and (lrg_status_return != "pending"):
        logging.error("The LRG status could not be resolved as public or pending")
        lrg_status_message = "ERROR: The LRG status could not be resolved with the webservice as public or pending"

    logging.info("LRG status is: " + lrg_status_return)
    logging.info(lrg_status_message)
    return lrg_status_return, lrg_status_message


def getLrgXml(lrg_id, lrg_status):
    '''
    Get the XML data for an LRG

    Input -
    lrg_id: String. An LRG ID to get the xml for. Must be in the 
      format LRG_<number>
    lrg_status: String. The release status of the LRG, either public or pending.

    Output -
    lrg_xml: String. String of entire LRG XML file.
    '''
    # get api, address is different depending on whether lrg is public or pending
    if lrg_status == 'public':
        xml_url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/{}.xml'.format(lrg_id)
    if lrg_status == 'pending':
        xml_url = 'http://ftp.ebi.ac.uk/pub/databases/lrgex/pending/{}.xml'.format(lrg_id)
    xml_response = requests.get(xml_url)
    assert xml_response.status_code == 200, 'Could not query the API, check your connection and try again.'

    # return response as string
    lrg_xml = xml_response.text
    return(lrg_xml)


def getLrgFromWeb(input_text):
    '''
    Main web API script, calls the functions above to:
    - if gene name/ transcript name, get lrg id, if already an lrg id then save lrg_id
    - check lrg exists
    - get status of the lrg
    - get xml file from the api
    '''
    lrg_id = getLrgId(input_text)
    assert checkLrgExists(lrg_id)
    lrg_status = checkCurrentLrgStatus(lrg_id)[0]
    lrg_xml = getLrgXml(lrg_id, lrg_status)

    return(lrg_xml)
