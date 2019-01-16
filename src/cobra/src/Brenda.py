#!/usr/bin/python
import string
import hashlib
# from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file

# setup some constants -
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
password = "<your password goes here>".encode('utf-8')
encoded_password = hashlib.sha256(password).hexdigest()
client = SOAPProxy(endpointURL)

# load the ec number file -
input_text_file = open("<path to file with ec numbers>", "r")
for line_raw in input_text_file:

    # Remove the \n, and split along the tab -
    value_array = line_raw.strip().split('\t')

    # element 1 = ec:number, so we need to split along the :
    ec_number = value_array[1].split(':')[1]

    # setup the parameters for the call -
    # ecNumber*1.1.1.1#organism*Homo sapiens#"
    parameters = "<your username goes here>,"+encoded_password+",ecNumber*"+ec_number+"#Escherichia coli#"

    # make the call -
    resultString = client.getTurnoverNumber(parameters)

    # message to user -
    msg = "Completed "+ec_number
    print(msg)

    # if we have data, then dump to disk -
    if len(resultString) != 0:
        # finally, we need to open a file, and dump the kinetic data to disk -
        file_path = "<path to output files>/ec_data."+ec_number+".dat"
        with open(file_path, 'w') as f:
            output_string = resultString.encode('utf-8')
            f.write(output_string)

# close the dile -
input_text_file.close()
