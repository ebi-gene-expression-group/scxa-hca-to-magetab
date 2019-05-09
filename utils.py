# Filename: utils.py

from datetime import datetime
import os
import fnmatch
import re
import yaml
import logging
import json
import smtplib
from orderedset import OrderedSet
from collections import OrderedDict
import requests, requests.packages.urllib3.util.retry
from email.mime.text import MIMEText

# Sending an email via smtplib occassionally throws a '[Errno 111] Connection refused' - the constants below controls
# how many times the code should re-try re-sending the email, and how long to sleep between each attempt (secs)
MAXIMUM_NUMBER_OF_EMAIL_ATTEMPTS = 5
EMAIL_RETRY_INTERVAL = 5

# Experiment errors class
class HCA2MagetabTranslationError(Exception):
    pass

# For HTTP semantics of re-try logic required see: https://dss.integration.data.humancellatlas.org/
class RetryPolicy(requests.packages.urllib3.util.retry.Retry):
    def __init__(self, retry_after_status_codes={301}, **kwargs):
        super(RetryPolicy, self).__init__(**kwargs)
        self.RETRY_AFTER_STATUS_CODES = frozenset(retry_after_status_codes | requests.packages.urllib3.util.retry.Retry.RETRY_AFTER_STATUS_CODES)

retry_policy = RetryPolicy(read=5, status=5, status_forcelist=frozenset({500, 502, 503, 504}))
s = requests.Session()
a = requests.adapters.HTTPAdapter(max_retries=retry_policy)
s.mount('https://', a)

# Constants

def unix_time_millis(dt):
    epoch = datetime.utcfromtimestamp(0)
    return int((dt - epoch).total_seconds() * 1000)

# Return data directory
def get_data_dir(process_name):
    return os.path.join(os.environ['ATLAS_PROD'], 'singlecell', 'experiment', process_name)

# Returns dict containing process_name's config in the same environment as this script
def get_config(process_name):
    with open(os.path.join(process_name + ".yml"), 'r') as stream:
        return yaml.load(stream, Loader=yaml.FullLoader)

# Get value for a given key in config
def get_val(config, key):
    return config[key]

# Create a logger writing to file file_path
def create_logger(process_name, mode):
    file_path = os.path.join(get_data_dir(process_name), 'logs', process_name + '.' + mode + '.' + datetime.now().strftime('%Y-%m-%d') + '.log')
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
# If no value is found, check if the required HCA schema location was affected by a migration - according to property_migrations
def get_hca_value(accession, magetab_label, hca_schema_path, hca_data_structure, project_bundle_uuids, logger, config, property_migrations, warn_of_missing_fields_in_hca_json):
    value = get_hca_value_for_path(hca_schema_path, hca_data_structure, logger, config, warn_of_missing_fields_in_hca_json, accession, magetab_label, project_bundle_uuids)
    if value == get_val(config, 'notfound'):
        # Now check to see if key may be missing due to a HCA property migration
        # Example describedBy value: https://schema.humancellatlas.org/type/project/9.0.3/project
        describedBy = get_val(config, 'hca_schema_version_field_name')
        if hca_schema_path[0] in hca_data_structure:
            # We need describedBy_version to be able to check if a schema migration applies; if hca_schema_path[0] cannot be found in hca_data_structure, we don't have describedBy_version.
            # In such cases, this bundle does not have data in schema: hca_schema_path[0]
            describedBy_version = hca_data_structure[hca_schema_path[0]][describedBy].split('/')[-2]
            migrated_schema_path = get_migrated_location(property_migrations, describedBy_version, hca_schema_path, logger)
            if migrated_schema_path:
                value = get_hca_value_for_path(migrated_schema_path, hca_data_structure, logger, config, warn_of_missing_fields_in_hca_json, accession, magetab_label, project_bundle_uuids)
    return value

# Retrieve value that corresponds to magetab_label in hca_data_structure, populated from json files in study_uuid's bundle_uuid.
# To access that value traverse hca_data_structure recursively using keys in list: hca_schema_path in turn.
def get_hca_value_for_path(hca_schema_path, hca_data_structure, logger, config, warn_of_missing_fields_in_hca_json, accession = None, magetab_label = None, project_bundle_uuids = None):
    leaf = None
    for key in hca_schema_path:
        if key in list(hca_data_structure.keys()):
            if isinstance(hca_data_structure[key], dict):
                # Continue traversing if hca_data_structure[key] is a Dict
                hca_data_structure = hca_data_structure[key]
            else:
                # leaf could be string or a list, but never a Dict
                leaf = hca_data_structure[key]
        else:
            leaf = get_val(config, 'notfound')
            if warn_of_missing_fields_in_hca_json:
                if accession and magetab_label and project_bundle_uuids:
                    context_info = " for accession: %s, magetab label: '%s' from HCA json for study: %s bundle (row): %s" % (accession, magetab_label, project_bundle_uuids[0], project_bundle_uuids[1])
                else:
                    context_info=""
                logger.warning("Key: '%s' (in %s) is missing%s" % (key, hca_schema_path, context_info))
            break
    return get_magetab_equivalent(magetab_label, leaf, config)

# Returns True if version1 is the same or lower than version2
def hca_versions_in_asc_order(version1, version2):
    return int(version1.replace('.','')) <= int(version2.replace('.',''))

# Retrieve from property_migrations a new path to which hca_schema_path was migrated - if one exists for version on or later than describedBy_version
def get_migrated_location(property_migrations, describedBy_version, hca_schema_path, logger):
    hca_key = hca_schema_path[0]
    source_schema = re.sub(r"\_\d+$", "", hca_key)
    property = '.'.join(hca_schema_path[1:])
    m = re.search(r"\_\d+$", hca_key)
    postfix = ''
    if m != None:
        postfix = m.group(0)
    migrated_schema_path = None
    for entry in property_migrations['migrations']:
        if source_schema == entry['source_schema'] and property == entry['property'] and hca_versions_in_asc_order(entry['effective_from'], describedBy_version):
            migrated_schema_path = [entry['target_schema'] + postfix] +  entry['replaced_by'].split('.')
            logger.warning("Replacing %s with %s due to property migration from %s onwards" % ('.'.join(hca_schema_path), '.'.join(migrated_schema_path), entry['effective_from']))
            break
    return migrated_schema_path

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

# Retrieve json from url via method (passing data in the call, if method == 'post')
# Returns a tuple: (<result (json)>, <returned headers dict (for 'post' only)>)
def get_remote_json(url, logger, method = 'get', data = None):
    result = None
    headers = None
    err_msg = None
    try:
        if method == 'get':
            r = s.get(url)
        elif method == 'post':
            r = s.post(url, data = json.dumps(data), headers = { "Accept" : "application/json", "Content-Type" : "application/json" })
        else:
            err_msg = 'Unknown HTTP request method: "' + url + '" : ' + method
        result = json.loads(r.text)
        headers = r.headers
    except requests.HTTPError as e:
        err_msg = 'HTTPError when retrieving url: "' + url + '" : ' + str(e.code)
        logger.error(err_msg)
    except requests.ConnectionError as e:
        err_msg = 'ConnectionError when retrieving url: "' + url + '" : ' + str(e.reason)
        logger.error(err_msg)
    if err_msg:
        raise HCA2MagetabTranslationError(err_msg)
    return (result, headers)

# Return the url for the next page of results - having extract it via 'Link' key from dict: headers (itself return by the previous api call)
def get_api_url_for_next_page(headers):
    url = None
    if 'Link' in headers.keys() and re.search(r"rel=\"next\"$", headers['Link']):
        m = re.search(r'^.*\<(https.*?)\>.*$', headers['Link'])
        if m:
            url = m.group(1)
    return url

# Retrieve the list of unique project uuids, corresponding to all HCA projects of type: technology
def get_project_uuid2accessions_for_technology(hca_api_url_root, logger, technology, ncbi_taxon_id, config):
    project_uuid2accession = {}
    smart_regex = re.compile('smart-.*$')
    tenxV2_regex = re.compile('10xV2')
    if smart_regex.match(technology):
        data = { "es_query": {
                     "query": {
                         "bool": {
                             "must": [
                                 # Smart - seq2
                                 { "match": { "files.library_preparation_protocol_json.library_construction_approach.ontology": "EFO:0008931" } },
                                 { "match": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": ncbi_taxon_id } }
                             ],
                             "should": [
                                 # fluorescence - activated cell sorting
                                 { "match": { "files.dissociation_protocol_json.dissociation_method.ontology": "EFO:0009108" } },
                                 # enzymatic dissociation
                                 { "match": { "files.dissociation_protocol_json.dissociation_method.ontology": "EFO:0009128"}},
                                 { "match": { "files.dissociation_protocol_json.dissociation_method.text": "mouth pipette" } }
                             ], "must_not": [
                                 { "match": { "files.analysis_process_json.process_type.text": "analysis" } },
                                 { "range": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": { "lt": ncbi_taxon_id } } },
                                 { "range": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": { "gt": ncbi_taxon_id } } }
                             ] } } } }
    elif tenxV2_regex.match(technology):
        data = { "es_query": {
                     "query": {
                         "bool": {
                             "must": [
                                 # 10X v2 sequencing
                                 # TODO: what is the EFO number for 10xV1* and 10xV3?
                                 { "match": { "files.library_preparation_protocol_json.library_construction_approach.ontology": "EFO:0009310" } },
                                 { "match": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": ncbi_taxon_id } }
                             ],
                             "must_not": [
                                 { "match": { "files.analysis_process_json.process_type.text": "analysis" } },
                                 { "range": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": { "lt": ncbi_taxon_id } } },
                                 { "range": { "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": { "gt": ncbi_taxon_id } } }
                             ] } } } }
    else:
        raise HCA2MagetabTranslationError('Unknown technology: %s' % technology)
    # page size: 500 is the maximum allowed by HCA
    url = '%s/%s' % (hca_api_url_root, 'search?output_format=raw&replica=aws&per_page=500')
    json = get_remote_json(url, logger, 'post', data)[0]
    for result in json['results']:
        for project_json in get_hca_value_for_path(get_val(config, 'hca_project_json_path'), result, logger, config, True):
            project_uuid = get_hca_value_for_path(get_val(config, 'hca_project_uuid_path'), project_json, logger, config, True)
            project_title = get_hca_value_for_path(get_val(config, 'hca_project_title_path'), project_json, logger, config, True)
            gxa_accessions = OrderedSet([])
            if re.search('Test *$', project_title):
                # Skip test data sets
                continue
            # E-AAAA-00 appears to be the default value HCA uses when no ArrayExpress accession is available
            hca_old_arrayexpress_label = get_val(config, 'hca_old_arrayexpress_label')
            hca_new_arrayexpress_label = get_val(config, 'hca_new_arrayexpress_label')
            if hca_old_arrayexpress_label in project_json.keys():
                if project_json[hca_old_arrayexpress_label] != 'E-AAAA-00':
                    gxa_accessions.add(project_json[hca_old_arrayexpress_label])
            elif hca_new_arrayexpress_label in project_json.keys():
                # For the reason for the loop below see a comment near hca_old_arrayexpress_label in hca2mtab.yml
                for gxa_accession in project_json[hca_new_arrayexpress_label]:
                    gxa_accessions.add(gxa_accession)
            hca_supplementary_links_label = get_val(config, 'hca_supplementary_links_label')
            if hca_supplementary_links_label in project_json.keys():
                for url in project_json[hca_supplementary_links_label]:
                    m = re.search(r'^.*?\/(E-\w{4}-\d+).*$', url)
                    if m:
                        gxa_accessions.add(m.group(1))
            # If HCA project uuid corresponds to multiple existing gxa accessions, we join them into a single label - the decision on how to split them into separate experiments is left to gxa curators
            if len(gxa_accessions) > 0:
                accession = ''.join(gxa_accessions)
            else:
                accession = None
            project_uuid2accession[project_uuid] = accession
    return project_uuid2accession

# Retrieve the set of unique project uuids, corresponding to all smart-seq2 and 10x projects in HCA.
def get_all_hca_project_uuids(hca_api_url_root, config, logger):
    project_uuid2accession = {}
    for technology in get_val(config, 'technology_mtab2hca').keys():
        for ncbi_taxon_id in [9606, 10090]:
            # human and mouse
            project_uuid2accession.update(get_project_uuid2accessions_for_technology(hca_api_url_root, logger, technology, ncbi_taxon_id, config))
    return project_uuid2accession

# Retrieve HCA project uuid from idf_file_path
def get_project_uuid_from_idf(idf_file_path):
    project_uuid = None
    with open(idf_file_path, 'r') as file:
        for line in file:
            if re.search(r'^Comment\[SecondaryAccession\]', line):
                project_uuid = re.sub('\n','',line.split('\t')[-1])
                break
    return project_uuid

# Return a mapping between hca project uuids and accessions - either already existing in gxa_dir (we're updating an
# experiment previously imported from HCA) or newly minted E-CAND-* accession (for experiments about to be imported from HCA for the first time)
def get_gxa_accession_for_project_uuid(gxa_dir, project_uuid2accession):
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

# For a given project uuid, retrieve dict: bundle uuid->list of (file_uuid, file_name) tuples
# N.B. if mode == 'test', retrieve only the first page of results
def get_bundle2metadata_files_for_project_uuid(project_uuid, hca_api_url_root, mode, logger, config):
    bundle2metadata_files = {}
    hca_project_uuid_elasticsearch_path = get_val(config, 'hca_project_uuid_elasticsearch_path')
    data = { "es_query": { "query": { "match": { hca_project_uuid_elasticsearch_path : project_uuid }}}}
    # Page size = 10 appears to be the maximum currently allowed
    url = '%s/%s' % (hca_api_url_root, 'search?output_format=raw&replica=aws&per_page=10')
    bundle_cnt = 0
    while url:
        bundle_cnt += 10
        print('.', end = '', flush = True)
        (json, headers) = get_remote_json(url, logger, 'post', data)
        for result in json['results']:
            bundle_url = result['bundle_url']
            metadata_files = []
            m = re.search(r'/bundles/([\w\d\-]+)\?', bundle_url)
            if m:
                bundle_uuid = m.group(1)
            else:
                err_msg = "Failed to retrieve bundle uuid for project uuid: %s from url: %s" % (project_uuid, bundle_url)
                logger.error(err_msg)
                raise HCA2MagetabTranslationError(err_msg)
            analysis_bundle = False
            for file_json in get_hca_value_for_path(get_val(config, 'hca_file_json_path'), result, logger, config, True):
                file_name = get_hca_value_for_path(get_val(config, 'hca_file_json_name_path'), file_json, logger, config, True)
                if re.search(r"" + get_val(config, 'hca_analysis_file_regex'), file_name):
                    analysis_bundle = True
                elif re.search(r'\.json$', file_name):
                    metadata_files.append((file_json['uuid'], file_name))
            if not analysis_bundle:
                # Exclude analysis bundles from being loaded into gxa
                bundle2metadata_files[bundle_uuid] = metadata_files
        url = get_api_url_for_next_page(headers)
        if mode == 'test' and bundle_cnt >= get_val(config, 'test_max_bundles'): 
            break
    return bundle2metadata_files

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
def email_report(body, subject, from, email_recipients):
    msg = MIMEText(body, "plain", "utf-8")
    msg['Subject'] = subject
    msg['From'] = from
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

# Retrieve into hca_json/hca_json_cache all the json files for bundle_uuid
def retrieve_hca_json_for_bundle(accession, bundle2metadata_files, project_bundle_uuids, hca_json, hca_json_cache, hca_api_url_root, config, logger):
    # Iterate over (file_uuid, file_name) tuples corresponding to bundle_uuid
    for (file_uuid, file_name) in bundle2metadata_files[project_bundle_uuids[1]]:
        file_type = file_name.split('.')[0]
        # Retrieve the json for file_uuid from cache if it's there; otherwise retrieve it from HCA DCC and add to the cache
        if file_uuid not in hca_json_cache.keys():
            file_json_url = '%s/files/%s?replica=aws' % (hca_api_url_root, file_uuid)
            logger.debug('Study uuid: %s ; bundle uuid: %s ; Accession: %s - About to retrieve file (of type: %s) : %s ' % (project_bundle_uuids[0], project_bundle_uuids[1], accession, file_type, file_json_url))
            row_data = get_remote_json(file_json_url, logger)[0]
            if not re.search(r"" + get_val(config, 'hca_json_files_excluded_from_cache_regex'), file_name):
                hca_json_cache[file_uuid] = row_data
            # Store retrieved json in hca_json[file_type]
            hca_json[file_type] = row_data
        else:
            # json is in cache - retrieve from there
            hca_json[file_type] = hca_json_cache[file_uuid]

        logger.debug("%s : %s --> %s" % (accession, file_name, hca_json[file_type]))


        
version = '0.1'
# End of utils.py