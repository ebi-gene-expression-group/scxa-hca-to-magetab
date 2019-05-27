# -*- coding: utf-8 -*-
# Filename: hcadam.py

# A data access module to HCA DCC - used for HCA project data discovery and json content retrieval

import requests, requests.packages.urllib3.util.retry
import re
import utils
import json
from datetime import datetime
from orderedset import OrderedSet

# Data retrieval from HCA error class
class HCARetrievalError(Exception):
    pass

# Cache of all json retrieved from HCA - to avoid repeating search queries
# Structure: project_uuid -> bundle_url -> file_type -> list of json objects
hca_json_cache = {}

# For HTTP semantics of re-try logic required see: https://dss.integration.data.humancellatlas.org/
class RetryPolicy(requests.packages.urllib3.util.retry.Retry):
    def __init__(self, retry_after_status_codes={301}, **kwargs):
        super(RetryPolicy, self).__init__(**kwargs)
        self.RETRY_AFTER_STATUS_CODES = frozenset(retry_after_status_codes | requests.packages.urllib3.util.retry.Retry.RETRY_AFTER_STATUS_CODES)

retry_policy = RetryPolicy(read=5, status=5, status_forcelist=frozenset({500, 502, 503, 504}))
s = requests.Session()
a = requests.adapters.HTTPAdapter(max_retries=retry_policy)
s.mount('https://', a)

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
        raise HCARetrievalError(err_msg)
    return (result, headers)

# Retrieve the list of unique project uuids, corresponding to all HCA projects of type: technology
def get_hca_projects_for_technology(hca_api_url_root, technology, ncbi_taxon_id, project_uuids_filter, already_imported_project_uuids, config, mode, logger):
    logger.info("About to retrieve HCA projects and their json content for technology: %s ncbi_taxon_id: %d" % (technology, ncbi_taxon_id))
    time_start = utils.unix_time_millis(datetime.now())
    project_uuids = set([])
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
        raise HCARetrievalError('Unknown technology: %s' % technology)
    # Page size = 10 appears to be the maximum currently allowed
    url = '%s/%s' % (hca_api_url_root, 'search?output_format=raw&replica=aws&per_page=10')
    bundle_cnt = 0
    json_files_cnt = 0
    while url:
        (json, headers) = get_remote_json(url, logger, 'post', data)
        if 'results' not in json:
            logger.info("No HCA projects were retrieved for technology: %s and ncbi_taxon_id: %d" % (technology, ncbi_taxon_id))
            break
        for bundle in json['results']:
            project_json_path = utils.get_val(config, 'hca_files_path') + utils.get_val(config, 'hca_project_json')
            project_json = get_hca_structure_for_path(project_json_path, bundle)[0]
            project_title = utils.get_hca_value(utils.get_val(config, 'hca_project_title_path'), project_json, logger, config, True)
            project_uuid = utils.get_hca_value(utils.get_val(config, 'hca_project_uuid_path'), project_json, logger, config, True)
            if re.search('Test *$', project_title) or (project_uuids_filter and project_uuid not in project_uuids_filter):
                # Skip test data sets, or sets not in the list specifically requested to be imported
                continue
            if project_uuid not in already_imported_project_uuids:
                # In new_only mode already_imported_project_uuids may not be empty - we only include project_uuids (and cache their json)
                # if they have not been imported already.
                project_uuids.add(project_uuid)
                json_files_cnt += add_bundle_to_json_cache(bundle, project_uuid, config, logger)
        bundle_cnt += 10
        if mode == 'test' and bundle_cnt >= utils.get_val(config, 'test_max_bundles'):
            break
        if bundle_cnt % 500 == 0:
            print("%d bundles retrieved (projects so far: %d)" % (bundle_cnt, len(project_uuids)), flush = True)
        url = get_api_url_for_next_page(headers)
    time_end = utils.unix_time_millis(datetime.now())
    duration = (time_end - time_start) / 1000 / 60
    logger.info("Populating json cache for technology %s ncbi_taxon_id %d took: %d mins (%d bundles; %d json objects)" % (technology, ncbi_taxon_id, duration, bundle_cnt, json_files_cnt))
    return project_uuids

# Return the url for the next page of results - having extract it via 'Link' key from dict: headers (itself return by the previous api call)
def get_api_url_for_next_page(headers):
    url = None
    if 'Link' in headers.keys() and re.search(r"rel=\"next\"$", headers['Link']):
        m = re.search(r'^.*\<(https.*?)\>.*$', headers['Link'])
        if m:
            url = m.group(1)
    return url

# Retrieve the set of unique project uuids, corresponding to all smart-seq2 and 10x projects in HCA.
def get_hca_project_uuid_to_import(hca_api_url_root, config, mode, project_uuids_filter, already_imported_project_uuids, logger):
    project_uuids = set([])
    for technology in utils.get_val(config, 'technology_mtab2hca').keys():
        for ncbi_taxon_id in [9606, 10090]:
            # human and mouse
            project_uuids = project_uuids.union(get_hca_projects_for_technology(hca_api_url_root, technology, ncbi_taxon_id, project_uuids_filter, already_imported_project_uuids, config, mode, logger))
    for project_uuid in project_uuids:
        logger.info("Populated json cache for %s with %d bundles" % (project_uuid, len(list(hca_json_cache[project_uuid].keys()))))
    return project_uuids

# Populate hca_json_cache with all the json objects in bundle (within project_uuid); returns number of json files cached
def add_bundle_to_json_cache(bundle, project_uuid, config, logger):
    json_files_cnt = 0
    if project_uuid not in hca_json_cache:
        hca_json_cache[project_uuid] = {}
    hca_files_path = utils.get_val(config, 'hca_files_path')
    bundle_url_schema_path = utils.get_val(config, 'hca_bundle_url')
    bundle_url = utils.get_hca_value(bundle_url_schema_path, bundle, logger, config, True)
    if bundle_url in hca_json_cache[project_uuid]:
        return
    else:
        hca_json_cache[project_uuid][bundle_url] = {}
    file_type2json_list = get_hca_structure_for_path(hca_files_path, bundle)
    if not file_type2json_list:
        err_msg = "Failed to retrieve json content for project: %s and bundle url: %s using path: %s" % (project_uuid, bundle_url, '.'.join(hca_files_path))
        logger.error(err_msg)
        raise HCARetrievalError(err_msg)

    for file_type in file_type2json_list.keys():
        if re.search(r'\_json$', file_type):
            schema_type = re.sub(r"\_json$", "", file_type)
            if re.search(r"" + utils.get_val(config, 'hca_analysis_file_regex'), schema_type):
                # This is an analysis bundle - don't add it to cache
                hca_json_cache[project_uuid].pop(bundle_url)
                break
            hca_json_cache[project_uuid][bundle_url][schema_type] = []
            for json_dict in file_type2json_list[file_type]:
                hca_json_cache[project_uuid][bundle_url][schema_type].append(json_dict)
                json_files_cnt += 1
        # Schema type-level assumption checking
        if violates_assumption_hca_schemas_with_one_json_per_bundle_expected(file_type, hca_json_cache[project_uuid][bundle_url][schema_type], config):
            logger.warning("File type: %s for project uuid: %s ; bundle url: %s ; violates the assumption that only one json file for that type should exist in a bundle" % (file_type, project_uuid, bundle_url))

    # Bundle-level assumption checking
    if violates_assumption_hca_schema_types_in_every_bundle(set(hca_json_cache[project_uuid][bundle_url].keys()), config):
        file_types_in_every_bundle = ','.join(utils.get_val(config, 'hca_file_types_in_every_bundle'))
        err_msg = "Bundle: %s in project uuid: %s violates the assumption that each bundle should contain _all_ of the following file types: %s" % (bundle_url, project_uuid, file_types_in_every_bundle)
        logger.error(err_msg)
        raise HCARetrievalError(err_msg)

    return json_files_cnt

# Retrieve a list of all accessions corresponding to project uuid - from project_json of one of project_uuid's bundles in hca_json_cache
def get_gxa_accession_for_project_uuid(project_uuid, config):
    gxa_accessions = set([])
    if project_uuid in hca_json_cache:
        first_bundle_url = list(hca_json_cache[project_uuid].keys())[0]
        project_json = hca_json_cache[project_uuid][first_bundle_url][utils.get_val(config, 'hca_project')][0]
        # E-AAAA-00 appears to be the default value HCA uses when no ArrayExpress accession is available
        hca_old_arrayexpress_label = utils.get_val(config, 'hca_old_arrayexpress_label')
        hca_new_arrayexpress_label = utils.get_val(config, 'hca_new_arrayexpress_label')
        if hca_old_arrayexpress_label in project_json.keys():
            if project_json[hca_old_arrayexpress_label] != 'E-AAAA-00':
                gxa_accessions.add(project_json[hca_old_arrayexpress_label])
        elif hca_new_arrayexpress_label in project_json.keys():
            # For the reason for the loop below see a comment near hca_old_arrayexpress_label in hca2mtab.yml
            for gxa_accession in project_json[hca_new_arrayexpress_label]:
                gxa_accessions.add(gxa_accession)
        hca_supplementary_links_label = utils.get_val(config, 'hca_supplementary_links_label')
        if hca_supplementary_links_label in project_json.keys():
            for url in project_json[hca_supplementary_links_label]:
                m = re.search(r'^.*?\/(E-\w{4}-\d+).*$', url)
                if m:
                    gxa_accessions.add(m.group(1))
        # If HCA project uuid corresponds to multiple existing gxa accessions, we join them into a single label - the decision on how to split them into separate experiments is left to gxa curators
        if len(gxa_accessions) > 0:
            return ''.join(gxa_accessions)
    return None

# Retrieve all HCA json content for project_uuid
def get_json_for_project_uuid(project_uuid):
    return hca_json_cache[project_uuid]

# Return True if file_type violates our assumption that there should exist only one json file for tha type in a HCA bundle
def violates_assumption_hca_schemas_with_one_json_per_bundle_expected(schema_type, json_list, config):
    if schema_type in utils.get_val(config, 'hca_schemas_with_one_json_per_bundle_expected'):
        return len(json_list) > 1
    return False

# Return True if not all HCA file types in hca_file_types_in_every_bundle config can be found in set_of_file_types_in_bundle; Return False otherwise.
def violates_assumption_hca_schema_types_in_every_bundle(set_of_file_types_in_bundle, config):
    hca_file_types_in_every_bundle = utils.get_val(config, 'hca_schema_types_in_every_bundle')
    return len(set_of_file_types_in_bundle.intersection(hca_file_types_in_every_bundle)) != len(hca_file_types_in_every_bundle)


# Retrieve a dict in hca_data_structure that corresponds to hca_schema_path. If no such dict is found, return an empty dict
def get_hca_structure_for_path(hca_schema_path, hca_data_structure):
    struct = hca_data_structure
    for key in hca_schema_path:
        if key in struct:
            struct = struct[key]
        else:
            struct = {}
    return struct

version = '0.1'
# End of utils.py