

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
    pass


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
    pass


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
    pass


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

