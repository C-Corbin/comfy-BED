import requests
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
    # if input is an lrg number, save the variable
    if input_text.startswith('LRG_'):
        lrg_id = input_text

    # if input isn't an lrg number, try to query by name to find lrg number
    else:
        try:
            name_query_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query=name:{}'.format(input_text)
            name_query_response = requests.get(name_query_url)
            assert name_query_response.status_code == 200, 'Could not query the API, check your connection and try again.'

            # loop through the response xml
            root = ET.fromstring(name_query_response.text)
            for child in root.iter('hitCount'):
                # check that there is exactly 1 entry returned
                assert child.text == '1'
            # if so, extract the lrg id and save as a variable
            for child in root.iter('entry'):
                lrg_id = child.get('id')

        except AssertionError:
            name_query_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg?query={}'.format(input_text)
            name_query_response = requests.get(name_query_url)
            assert name_query_response.status_code == 200, 'Could not query the API, check your connection and try again.'

            # loop through the response xml
            root = ET.fromstring(name_query_response.text)
            for child in root.iter('hitCount'):
                # check that there is exactly 1 entry returned
                assert child.text == '1'
            # if so, extract the lrg id and save as a variable
            for child in root.iter('entry'):
                lrg_id = child.get('id')
                
        except:
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


def getLrgStatus(lrg_id):
    '''
    Get the release status of an LRG (either public or pending). If an 
    LRG is pending, thrown a warning to use the data with caution since
    the contents of the LRG are yet to be confirmed.

    Input -
    lrg_id: String. An LRG ID to check the status of. Must be in the 
      format LRG_<number>

    Output -
    lrg_status: String. Pulled from the API response, either public or pending.
    '''
    # query api, returns xml with status nested within it
    status_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/lrg/entry/{}?fields=status'.format(lrg_id)
    status_response = requests.get(status_url)
    assert status_response.status_code == 200, 'Could not query the API, check your connection and try again.'

    # extract status from the xml
    root = ET.fromstring(status_response.text)
    for child in root.iter('value'):
        lrg_status = child.text

    #TODO Add warning if pending

    return(lrg_status)


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
    pass


def getLrgWeb(input_text):
    '''
    - if gene name/ transcript name, get lrg id, if already an lrg id then save lrg_id
    - check lrg exists
    - get status of the lrg
    - get xml file from the api
    '''
    pass

