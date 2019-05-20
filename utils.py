# Filename: utils.py

from datetime import datetime
import os
import fnmatch
import re
import yaml
import logging
import smtplib
from orderedset import OrderedSet
from collections import OrderedDict

from email.mime.text import MIMEText

# Constants
# Sending an email via smtplib occassionally throws a '[Errno 111] Connection refused' - the constants below controls
# how many times the code should re-try re-sending the email, and how long to sleep between each attempt (secs)
MAXIMUM_NUMBER_OF_EMAIL_ATTEMPTS = 5
EMAIL_RETRY_INTERVAL = 5

# HCA to MAGETAB translation errors class
class HCA2MagetabTranslationError(Exception):
    pass

def unix_time_millis(dt):
    epoch = datetime.utcfromtimestamp(0)
    return int((dt - epoch).total_seconds() * 1000)

# Returns dict containing process_name's config in the same environment as this script
def get_config(process_name):
    with open(os.path.join(process_name + ".yml"), 'r') as stream:
        return yaml.load(stream, Loader=yaml.FullLoader)

# Get value for a given key in config
def get_val(config, key):
    return config[key]

# Create a logger writing to file file_path
def create_logger(data_dir, process_name, mode):
    file_path = os.path.join(data_dir, 'logs', process_name + '.' + mode + '.' + datetime.now().strftime('%Y-%m-%d') + '.log')
    # Need to give getLogger file_path (name) argument - so that it creates a unique logger for a given file_path
    logger = logging.getLogger(file_path)
    logger.setLevel(logging.INFO)
    hdlr = logging.FileHandler(file_path)
    formatter = logging.Formatter('%(asctime)s %(levelname)s : %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    return logger

# Translate hca_value to its equivalent in magetab, specific to magetab_label. If no such equivalent was found in config, return hca_value
def get_magetab_equivalent(magetab_label, hca_value, config):
    hca2mtab_dict = get_val(config, 'cv_translate')
    if not isinstance(magetab_label, list) and magetab_label in hca2mtab_dict.keys():
        sdrf_col_specific_translations = hca2mtab_dict.get(magetab_label)
        if hca_value in sdrf_col_specific_translations.keys():
            magetab_value = sdrf_col_specific_translations.get(hca_value)
        elif 'default' in sdrf_col_specific_translations.keys():
            # Used for SDRF columns with binary values, e.g. 'Comment[LIBRARY_LAYOUT]'
            magetab_value = sdrf_col_specific_translations.get('default')
        else:
            magetab_value = hca_value
    else:
        magetab_value = hca_value
    return magetab_value

# Retrieve value that corresponds to magetab_label in hca_data_structure, populated from json files in study_uuid's bundle_uuid.
# To access that value traverse hca_data_structure recursively using keys in list: hca_schema_path in turn.
def get_hca_value(hca_schema_path, hca_data_structure, logger, config, warn_of_missing_fields_in_hca_json, magetab_label = None, context = None):
    leaf = None
    struct = hca_data_structure
    for key in hca_schema_path:
        if key in struct:
            if isinstance(struct[key], dict):
                # Continue traversing if struct[key] is a Dict
                struct = struct[key]
            else:
                # leaf could be string or a list, but never a Dict
                leaf = struct[key]
        else:
            leaf = get_val(config, 'notfound')
            if warn_of_missing_fields_in_hca_json:
                if context and magetab_label:
                    context_info = " for accession: %s, magetab label: '%s' from HCA json for study: %s bundle (row): %s" % (context[0], magetab_label, context[1], context[2])
                else:
                    context_info=""
                logger.warning("Key: '%s' (in %s) is missing%s" % (key, hca_schema_path, context_info))
            break
    return get_magetab_equivalent(magetab_label, leaf, config)

# Returns True if version1 is the same or lower than version2
def hca_versions_in_asc_order(version1, version2):
    return int(version1.replace('.','')) <= int(version2.replace('.',''))

# Add characteristic->value for characteristic to dict: characteristic_values
# N.B. The values in characteristic_values are sets, as this dict is later used in deciding which characteristic should
# also be a factor and thus we're only interested if count(distinct values) is greater than 1
def store_characteristic_value(characteristic, value, characteristic_values):
    if characteristic not in characteristic_values.keys():
        characteristic_values[characteristic] = set([value])
    else:
        characteristic_values[characteristic].add(value)

# 1. Append column name: sdrf_column_header to sdrf_column_headers
# 2. If sdrf_column_header is a characteristic with a meaningful value (i.e. != get_val(config, 'curate') which is an instruction for Atlas curators to fix it),
#    add it to characteristic_values (in turn used later to decide which characteristic should also be a factor)
# 3. If the value != get_val(config, 'notfound'), register the corresponding column as non-empty
def add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, sdrf_column_header, value, row, characteristic_values, config):
    sdrf_column_headers.append(sdrf_column_header)
    if re.search("Characteristic", sdrf_column_header) and value not in [get_val(config, 'curate'), get_val(config, 'notfound')]:
        store_characteristic_value(sdrf_column_header, value, characteristic_values)
    if value == get_val(config, 'notfound'):
        # On curators request, never output the get_val(config, 'notfound'), but empty string instead
        value = ''
    if value != '':
        # Record that sdrf column index == len(sdrf_column_headers) is not empty
        indexes_of_non_empty_sdrf_columns.add(len(sdrf_column_headers)-1)
    row.append(value)

# Return the url for the next page of results - having extract it via 'Link' key from dict: headers (itself return by the previous api call)
def get_api_url_for_next_page(headers):
    url = None
    if 'Link' in headers.keys() and re.search(r"rel=\"next\"$", headers['Link']):
        m = re.search(r'^.*\<(https.*?)\>.*$', headers['Link'])
        if m:
            url = m.group(1)
    return url

# Retrieve HCA project uuid from idf_file_path
def get_project_uuid_from_idf(idf_file_path):
    project_uuid = None
    with open(idf_file_path, 'r') as file:
        for line in file:
            if re.search(r'^Comment\[SecondaryAccession\]', line):
                project_uuid = re.sub('\n','',line.split('\t')[-1])
                break
    return project_uuid


# Return a list of unique project uuids already imported from HCA
def get_previously_imported_projects(data_dir):
    project_uuids = set([])
    gxa_acc_regex_obj = re.compile('(E-\w{4}-\d+)+\.idf\.txt')
    for file_name in os.listdir(data_dir):
        if gxa_acc_regex_obj.match(file_name):
            project_uuid = get_project_uuid_from_idf(os.path.join(data_dir, file_name))
            project_uuids.add(project_uuid)
    return list(project_uuids)

# Return a mapping between hca project uuids and accessions - either already existing in gxa_dir (we're updating an
# experiment previously imported from HCA) or newly minted E-CAND-* accession (for experiments about to be imported from HCA for the first time)
def resolve_gxa_accession_for_project_uuid(gxa_dir, project_uuid2accession):
    # An auxiliary list used to enable minting of E-CAND-* accessions
    candidate_experiments = []

    # First go through all the idf files in gxa_dir to find uuids of projects that were previously imported from HCA
    # NB. Each magetab file name can consist of more than one consecutive gxa accession. This can happen only for multiple-technology HCA projects that correspond
    # to more than one (single-technology) _existing_ gxa experiment.
    gxa_acc_regex_obj = re.compile('(E-\w{4}-\d+)+\.idf\.txt')
    candidate_acc_regex_obj = re.compile('E-CAND-\d+\.idf\.txt')
    for file_name in os.listdir(gxa_dir):
        if gxa_acc_regex_obj.match(file_name):
            project_uuid = get_project_uuid_from_idf(os.path.join(gxa_dir, file_name))
            if project_uuid in project_uuid2accession.keys():
                if project_uuid2accession[project_uuid] == None:
                    # Since gxa accession was not found (by get_project_uuid2accessions_for_technology() call) in project_uuid's json,
                    # assign to project_uuid the accession of the existing idf file in which this project_uuid was found.
                    project_uuid2accession[project_uuid] = re.sub(r'\.idf\.txt','', file_name)
            if candidate_acc_regex_obj.match(file_name):
                candidate_experiments.append(file_name)

    # Now mint accessions for all the new experiments
    if len(candidate_experiments) > 0:
        maximum_candidate_exp_num = int(re.sub('\.idf\.txt|E\-CAND\-', '', sorted(candidate_experiments)[-1]))
    else:
        maximum_candidate_exp_num = 0
    for project_uuid in project_uuid2accession.keys():
        if project_uuid2accession[project_uuid] == None:
            project_uuid2accession[project_uuid] = "E-CAND-%d" % (maximum_candidate_exp_num + 1)
            maximum_candidate_exp_num += 1

    return project_uuid2accession

# Return True if it would be 'valid' to insert col_name to the end of sdrf_column_headers.
# The position is deemed 'valid' if:
# 1. The relevant config doesn't say anything to the contrary (i.e. col_name not in sdrf_colkey_after_colvalues.keys()), or:
# 2. The column that comes just before (sdrf_column_headers[-1]) is explicitly listed as an allowed predecessor in sdrf_colkey_after_colvalues[col_name]
def position_valid_for_sdrf_column(col_name, sdrf_column_headers, config):
    sdrf_colkey_after_colvalues = get_val(config, 'sdrf_colkey_after_colvalues')
    if col_name not in sdrf_colkey_after_colvalues.keys() or sdrf_column_headers[-1] in sdrf_colkey_after_colvalues[col_name]:
        return True
    return False

# Remove elements of row corresponding to indexes in indexes_of_empty_columns
def remove_empty_columns(row, indexes_of_empty_columns):
    position_shift = 0
    for idx in sorted(indexes_of_empty_columns):
        row.pop(idx - position_shift)
        position_shift += 1

# Find out number of protocols per each type
def get_number_of_protocols_per_type(protocol_type2protocols):
    protocol_type2num_protocols = OrderedDict()
    for protocol_type in protocol_type2protocols.keys():
        protocol_type2num_protocols[protocol_type] = len(protocol_type2protocols[protocol_type])
    return protocol_type2num_protocols

# It is possible for more than one protocol of a given type to be:
# 1. used in a single experiment-technology
# 2. be applied to the same cell/sample
# This function expands the number of 'Protocol REF' columns in sdrf to be able to express the above potential variability by representing one protocol per 'Protocol REF' column
def expand_protocol_columns(row, headers, protocol_type2num_protocols, logger):
    protocol_type_idx = 0
    col_idx = 0
    protocol_type2values = {}
    protocol_types = list(protocol_type2num_protocols.keys())
    # We will be expanding headers list - hence we need to be iterating over its (original) copy
    original_headers = headers.copy()
    
    for header in original_headers:
        if header == 'Protocol REF':
            # Assumption: the first 'Protocol REF' sdrf column corresponds to protocol_types[0], the second to protocol_types[1], and so on.
            protocol_type = protocol_types[protocol_type_idx]
            # We subtract 1 for the one column we already have for that protocol_type
            num_new_columns_to_be_added = protocol_type2num_protocols[protocol_type] - 1
            # Retrieve all values for protocol_type into a list: protocol_type2values[protocol_type]
            if row:
                # row[col_idx] is a comma-separated string of values - one per each protocol (of protocol_type) that was applied to sample/cell represented by the current sdrf row
                protocol_type2values[protocol_type] = row[col_idx].split(',')
                
            while num_new_columns_to_be_added > 0 : 
                if row:
                    # We're expanding protocol columns in an sdrf data line
                    # If any values still remain in protocol_type2values[protocol_type], insert the first of them into the current 'Protocol REF' column
                    if len(protocol_type2values[protocol_type]) > 0:
                        row.insert(col_idx, protocol_type2values[protocol_type][0])
                        protocol_type2values[protocol_type].pop(0)
                else:
                    # We're expanding protocol columns in the sdrf headers line
                    headers.insert(col_idx, 'Protocol REF')
                col_idx += 1
                num_new_columns_to_be_added -= 1
                
            # Insert the last remaining protocol value (or '' otherwise) into the final column for this protocol_type    
            if row:
                if len(protocol_type2values[protocol_type]) > 0:
                    row[col_idx] = protocol_type2values[protocol_type][0]
                    protocol_type2values[protocol_type].pop(0)
                else:
                    row[col_idx] = ''
                    
            protocol_type_idx += 1
        col_idx += 1

# Email report on all experiments imported from HCA DCC        
def email_report(subject, body, sender, email_recipients):
    msg = MIMEText(body, "plain", "utf-8")
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = email_recipients

    for attempt in range(MAXIMUM_NUMBER_OF_EMAIL_ATTEMPTS):
        try:
            # Send the message via the local smtp server
            s = smtplib.SMTP('localhost')
            s.sendmail(msg['From'], email_recipients, msg.as_string())
            s.quit()
            break
        except EnvironmentError as exc:
            if exc.errno == errno.ECONNREFUSED:
                time.sleep(EMAIL_RETRY_INTERVAL)
            else:
                raise # re-raise otherwise
    else: # we never broke out of the for loop
        raise RuntimeError("Maximum number of unsuccessful attempts to email validation report reached for imported experiments report")

# Retrieve ArrayExpress single-cell technology name, corresponding to hca_technology
def get_gxa_technology(hca_technology, config):
    # technology_mtab2hca list contains the mapping between a single-cell technlogy valid in ArrayExpress
    # (c.f. https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?spaceKey=GXA&title=6.+Single-Cell+Curation+Guide)
    # and the list of its synonyms in HCA (for now we're assuming there could be more than one synonym)
    technology_mtab2hca = get_val(config, 'technology_mtab2hca')
    gxa_technology = None
    for key in technology_mtab2hca:
        for val in technology_mtab2hca[key]:
            hca_technology_regex = re.compile(val)
            if hca_technology_regex.match(hca_technology):
                gxa_technology = key
                break
    return gxa_technology
        
version = '0.1'
# End of utils.py