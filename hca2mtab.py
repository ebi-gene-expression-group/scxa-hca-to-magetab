# -*- coding: utf-8 -*-
# Filename: hca2mtab.py

# A module to convert single-cell experimental metadata from HCA JSON to Magetab format

import re
import os
import sys
import unicodecsv as csv
from datetime import datetime
from collections import OrderedDict
from orderedset import OrderedSet
# For checking existence of magetab files
import glob
# Local modules
import utils
import hcadam

process_name = 'hca2mtab'

# Converts to magetab files bundles belonging to HCA uuid in the first element in each tuple in worklist.
# The second element in each tuple is the accession to be used on the gxa side.
# mode can be 'test' o 'prod':
# - 'test' retrieves just N bundle per project id, where N is defined in test_max_bundles field in config
# - 'prod' is the norma production run
# data_dir - directory in which the generated magetab files should be placed
def convert_hca_json_to_magetab(mode, data_dir, project_uuids_filter = None, new_only = True, sender = None, email_recipients = None ):
    # Retrieve the HCA Json to MAGETAB translation config
    config = utils.get_config(process_name)
    idf_config = utils.get_val(config, 'idf')
    sdrf_config = utils.get_val(config, 'sdrf')

    logger = utils.create_logger(data_dir, process_name, mode)
    hca_api_url_root = utils.get_val(config, 'hca_api_url_root')
    # already_imported_project_uuids will be excluded from the import (and their json will not be cached)
    if new_only:
        already_imported_project_uuids = utils.get_previously_imported_projects(data_dir)
    else:
        already_imported_project_uuids = []

    project_uuids = hcadam.get_hca_project_uuid_to_import(hca_api_url_root, config, mode, project_uuids_filter, already_imported_project_uuids, logger)

    # project_uuid2gxa_accession dict forms the worklist of experiments to be imported from HCA
    project_uuid2gxa_accession = {}
    for project_uuid in project_uuids:
        project_uuid2gxa_accession[project_uuid] = hcadam.get_gxa_accession_for_project_uuid(project_uuid, config)
    project_uuid2gxa_accession = utils.resolve_gxa_accession_for_project_uuid(data_dir, project_uuid2gxa_accession)

    # Experiments imported from HCA DCC - for email report
    imported_experiments = []
    
    # Log experiments to be imported
    logger.info("About to import from HCA DCC the following experiments:")
    for project_uuid in project_uuid2gxa_accession.keys():
        logger.info("%s -> %s" % (project_uuid, project_uuid2gxa_accession[project_uuid]))

    # Metadata retrieve starts here
    for project_uuid in project_uuid2gxa_accession.keys():
        time_start = utils.unix_time_millis(datetime.now())
        accession = project_uuid2gxa_accession.get(project_uuid)
        if new_only:
            # N.B. if new_only is True, HCA projects for which an idf file in data_dir doesn't exist will be imported
            idf_file_path = '%s/%s.idf.txt*' % (data_dir, accession)
            if glob.glob(idf_file_path):
                logger.info("Not importing %s as %s already exists (new_only mode: %s)" % (accession, idf_file_path, str(new_only)))
                continue
        else:
            logger.info('About to translate json for HCA study uuid: %s to magetab for gxa accession: %s' % (project_uuid, accession) )

        # Retrieve all HCA json content for project_uuid
        hca_json_for_project_uuid = hcadam.get_json_for_project_uuid(project_uuid)

        # Initialise SDRF-related data structures and flags
        # Set of technologies found in bundles for a given project uuid. The presence of a technology name in that set acts as a flag that sdrf column headers have been collected for that technology.
        technologies_found = set([])
        
        # List of SDRF column headers (per technology) that will be output in each (technology-specific) sdrf file
        technology2sdrf_column_headers = {}
        # List value corresponding to each technology key in technology2rows dict will be used to accumulate rows of data to be output into the generated (technology-specific) SDRF file
        technology2rows = {}
        # For a given technology key, after all the bundles for a given project have been seen, the value (set) indexes of sdrf columns that are empty for this technology
        # (and therefore will be removed before the sdrf matrix is output into the sdrf file)
        # N.B. Before any bundles are seen, all columns are assumed to be empty until at least one value is encountered for each.
        technology2indexes_of_empty_columns = {}
        # Initialise IDF-related data structures (for data such as protocols - that need to be collected from all the bundles)
        # technology2protocol_type2protocols is used to store all protocol names - to be populated later inside IDF file
        technology2protocol_type2protocols = {}
        # technology2protocol_type2max_protocol_num_per_sample stores maximum number of protocols per technology-protocol_type in any given sample/bundle.
        # This number will dictate how many 'Protocol REF' columns should be output for that protocol_type in sdrf file for that technology
        technology2protocol_type2max_protocol_num_per_sample = {}

        # characteristic_values_in_bundle dict stores sets of (unique) values per characteristic found - in order
        # to later automatically generate the corresponding Factors - for all characteristics for which the values change across the experiment (across all technologies).
        # N.B. A simplifying assumption is made here that in a multi-technology experiment, each technology-specific portion will get the same Factors
        characteristic_values = OrderedDict()

        # Auxiliary counter - used to limit number of HCA bundles processed during testing
        bundle_cnt = 0
        for bundle_url in hca_json_for_project_uuid.keys():
            # We want to warn of missing fields for the first bundle (since each bundle will contain some technology), the test below
            # effectively checks if we're dealing with the first bundle or not
            warn_of_missing_fields = not technologies_found
            hca_json_for_bundle = hca_json_for_project_uuid[bundle_url]
            context = (accession, project_uuid, bundle_url)
            ####################################################
            #### Collect protocols for IDF from bundle_uuid ####
            ####################################################
            protocol_type2protocols_in_bundle = OrderedDict([])
            for protocol_key in utils.get_val(config, 'protocol_types'):
                protocol_type2protocols_in_bundle[protocol_key] = OrderedSet([])
                for schema_type in list(hca_json_for_bundle.keys()):
                    if re.search(r"" + protocol_key, schema_type):
                        for protocol_json in hca_json_for_bundle[schema_type]:
                            protocol_name = utils.get_hca_value(utils.get_val(config, 'hca_protocol_name_path'), protocol_json, logger, config, False, 'Protocol Name', context)
                            if protocol_name != utils.get_val(config, 'notfound'):
                                protocol_description = utils.get_hca_value(utils.get_val(config, 'hca_protocol_description_path'), protocol_json, logger, config, False, 'Protocol Description', context)
                                protocol_type = utils.get_hca_value(utils.get_val(config, 'hca_protocol_type_path'), protocol_json, logger, config, False, 'Protocol Type', context)
                                protocol_type2protocols_in_bundle[protocol_key].add((protocol_name, protocol_description, protocol_type))

            ##################
            ###### SDRF ######
            ##################
            technology = None
            # Set of indexes of sdrf columns with non-empty sdrf columns for the current bundle_uuid
            indexes_of_non_empty_sdrf_columns = set([])
            # Output one SDRF row per each donor - sequence file tuple in the bundle
            # Assumptions relating to donors:
            # 1. Every HCA bundle has at least one json object for both: donor_organism and cell_suspension
            # 2. When multiple donor_organism and cell_suspension json objects exist, in the lists of json objects for donor_organism and
            #    cell_suspension respectively, the first JSON in donor_organism list corresponds to the first JSON in the cell_suspension list, and so on.
            #    However, in multi-donor samples with just one cell_suspension json object (e.g. project_uuid: d96c2451-6e22-441f-a3e6-70fd0878bb1b,
            #    bundle_url: https://dss.data.humancellatlas.org/v1/bundles/fb64e4f9-9a24-4a6a-856f-2b7c0d4f309d?version=2019-01-03T153203.452910Z&replica=aws
            #    that single cell_suspension json is assumed to apply to all donor_organism json objects in that bundle.
            donor_json_list = hca_json_for_bundle[utils.get_val(config, 'hca_donor_organism')]
            cell_suspension_json_list = hca_json_for_bundle[utils.get_val(config, 'hca_cell_suspension')]
            if len(cell_suspension_json_list) != len(donor_json_list) and len(cell_suspension_json_list) != 1:
                err_msg = " Project: %s bundle: %s contain multiple donor_organism and cell_suspension jsons, but their number is not the same" % (project_uuid, bundle_url)
                logger.error(err_msg)
                raise utils.HCA2MagetabTranslationError(err_msg)

            i = 0
            while i < len(donor_json_list):
                donor_json = donor_json_list[i]
                if len(cell_suspension_json_list) > 1:
                    cell_suspension_json = cell_suspension_json_list[i]
                else:
                    cell_suspension_json = cell_suspension_json_list[0]
                i += 1

                for datafile_json in hca_json_for_bundle[utils.get_val(config, 'hca_sequence_file')]:
                    sdrf_column_headers = []

                    current_row = []
                    for line_item in sdrf_config:
                        magetab_label = line_item[0]
                        hca_path = line_item[1]
                        if isinstance(hca_path, list):
                            if not utils.position_valid_for_sdrf_column(magetab_label, sdrf_column_headers, config):
                                # Skip sdrf columns if the position in which they would be inserted would not be valid given the column just before: sdrf_column_headers[-1]
                                continue
                            elif magetab_label in ['Characteristics[geographical location]', 'Characteristics[genotype]']:
                                # Special handling/parsing - geographical location - multiple json files need checking for field presence
                                value = utils.get_val(config, 'notfound')
                                regex = hca_path[0]
                                for schema_type in list(hca_json_for_bundle.keys()):
                                    if re.search(r"" + regex, schema_type):
                                        for json_dict in hca_json_for_bundle[schema_type]:
                                            value = utils.get_hca_value(hca_path[1:], json_dict, logger, config, warn_of_missing_fields, magetab_label, context)
                                            if value != utils.get_val(config, 'notfound'):
                                                break
                                utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, value, current_row, characteristic_values, config)
                            elif magetab_label == 'Protocol REF':
                                protocol_type = hca_path[0]
                                # TODO:
                                # Note that before sdrf is output, we will need to split protocol_ids into separate columns, but not before processing all the bundles in the project - we have to wait till the
                                # end to we know how many protocols per technology-protocol type we have. _In theory_ we could have 3 different enrichment protocols for a given technology in a project, e.g.
                                # FACS3, FACS5 and FACS8, and for the current bundle protocol_ids = 'FACS3,FACS8'. Then before outputting sdrf we will have to 'explode' the 'Protocol REF' column corresponding
                                # to the enrichment protocol into 3 (tab-delimited) new columns - and these columns for the current bundle_uuid will have values: 'FACS3\tFACS8\t' and
                                # headers: 'Protocol REF\tProtocol REF\tProtocol REF'
                                protocol_ids = ','.join([ x[0] for x in list(protocol_type2protocols_in_bundle[protocol_type]) ])
                                utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, protocol_ids, current_row, characteristic_values, config)
                            elif len(hca_path) > 0 and re.search(r"" + utils.get_val(config, 'hca_protocol_schema_regex'), hca_path[0]):
                                protocol_type = hca_path[0]
                                # Special handling/parsing - for a given protocol_type, various protocol-related information needs to be collected from potentially multiple HCA json files
                                values = set([])
                                for schema_type in [x for x in hca_json_for_bundle.keys() if x == protocol_type]:
                                    for json_dict in hca_json_for_bundle[schema_type]:
                                        value = utils.get_hca_value(hca_path[1:], json_dict, logger, config, warn_of_missing_fields, magetab_label, context)
                                        if value != utils.get_val(config, 'notfound'):
                                            if magetab_label == 'Comment[library construction]':
                                                # Capture technology for the current bundle
                                                hca_technology = value.lower()
                                                technology = utils.get_gxa_technology(hca_technology, config)
                                                value = technology
                                            values.add(str(value))
                                utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, ', '.join(values), current_row, characteristic_values, config)
                            elif magetab_label == 'Comment[HCA bundle url]':
                                utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, bundle_url, current_row, characteristic_values, config)
                            elif magetab_label in ['Comment[RUN]', 'Comment[FASTQ_URI]', 'Scan Name', 'Comment[technical replicate group]', 'Comment[HCA file uuid]']:
                                # Special handling/parsing - Comment[RUN] - datafile_key json file need checking for field presence
                                value = utils.get_hca_value(hca_path, datafile_json, logger, config, warn_of_missing_fields, magetab_label, context)
                                if magetab_label == 'Comment[RUN]':
                                    # NB. We're stripping e.g. _2.fastq.gz from the end - to retain just the core file name
                                    # Tested on the following types of file names:
                                    # "FCA7167226_I1.fastq.gz", "MantonBM7_HiSeq_4_S19_L005_R2_001.fastq.gz", "E18_20160930_Neurons_Sample_57_S054_L005_I1_010.fastq.gz", "FCA7167226.fastq.gz"
                                    value = re.sub(r"(\_\w\d|\_\w\d\_\d+|\_\d)*\.f\w+\.gz", "", value)
                                utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, value, current_row, characteristic_values, config)
                            else:
                                schema_type = hca_path[0]
                                if schema_type in hca_json_for_bundle:
                                    if schema_type == utils.get_val(config, 'hca_donor_organism'):
                                        json_dict = donor_json
                                    elif schema_type == utils.get_val(config, 'hca_cell_suspension'):
                                        json_dict = cell_suspension_json
                                    else:
                                        # Retrieving the first element below follows the assumption of one single json object in schema_type in a bundle
                                        # (all the special cases were handled above)
                                        json_dict = hca_json_for_bundle[schema_type][0]
                                    value = utils.get_hca_value(hca_path[1:], json_dict, logger, config, warn_of_missing_fields, magetab_label, context)
                                else:
                                    value = utils.get_val(config, 'notfound')
                                if magetab_label in \
                                    ['Characteristics[organism]', 'Characteristics[disease]', 'Characteristics[cell subtype]', 'Characteristics[ethnic group]','Characteristics[strain]'] \
                                    and value != utils.get_val(config, 'notfound'):
                                    # Special handling/parsing - organism, disease - could be multiple according to HCA schema
                                    utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, ','.join([x['text'] for x in value]), current_row, characteristic_values, config)
                                else:
                                    # magetab_label is not a list or a special case
                                    utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, str(value), current_row, characteristic_values, config)
                        else:
                            # hca_path is not a list - add to the row as is
                            utils.add_to_row(indexes_of_non_empty_sdrf_columns, sdrf_column_headers, magetab_label, hca_path, current_row, characteristic_values, config)
                    # At least one bundle has been seen - the SDRF columns have now been determined
                    if technology:
                        # Append current_row to the list of rows in the SDRF file being generated
                        if technology not in technology2rows.keys():
                            technology2rows[technology] = []
                        technology2rows[technology].append(current_row)
                        # The presence of a technology name in that set acts as a flag that sdrf column headers have been collected for that technology.
                        if technology not in technologies_found:
                            technology2sdrf_column_headers[technology] = sdrf_column_headers
                            # To start off with assume all columns are empty
                            technology2indexes_of_empty_columns[technology] = range(len(sdrf_config))
                            # Initialise technology2protocol_type2protocols with new technology
                            technology2protocol_type2protocols[technology] = OrderedDict()
                        technologies_found.add(technology)
                        # Store (without duplicates) for technology the protocols found for bundle_uuid (i.e. those in protocol_type2protocols_in_bundle)
                        for protocol_type in protocol_type2protocols_in_bundle.keys():
                            num_protocols_in_bundle = len(protocol_type2protocols_in_bundle[protocol_type])
                            if num_protocols_in_bundle > 0:
                                if technology not in technology2protocol_type2max_protocol_num_per_sample.keys():
                                    technology2protocol_type2max_protocol_num_per_sample[technology] = OrderedDict({ protocol_type : num_protocols_in_bundle })
                                elif protocol_type not in technology2protocol_type2max_protocol_num_per_sample[technology].keys():
                                    technology2protocol_type2max_protocol_num_per_sample[technology][protocol_type] = num_protocols_in_bundle
                                else:
                                    technology2protocol_type2max_protocol_num_per_sample[technology][protocol_type] = max(num_protocols_in_bundle, technology2protocol_type2max_protocol_num_per_sample[technology][protocol_type])
                                if protocol_type not in technology2protocol_type2protocols[technology].keys():
                                    technology2protocol_type2protocols[technology][protocol_type] = OrderedSet([])
                                # Merge set: protocol_type2protocols_in_bundle[protocol_type] into set already in technology2protocol_type2protocols[technology][protocol_type]
                                technology2protocol_type2protocols[technology][protocol_type] |= protocol_type2protocols_in_bundle[protocol_type]
                    else:
                        err_msg = "Failed to retrieve valid technology from value: \"%s\" in bundle: %s" % (hca_technology, bundle_url)
                        logger.error(err_msg)
                        raise utils.HCA2MagetabTranslationError(err_msg)

            # Now remove from technology2indexes_of_empty_columns[technology] all column indexes we found non-empty values for, for the current bundle_uuid
            technology2indexes_of_empty_columns[technology] = [x for x in technology2indexes_of_empty_columns[technology] if x not in indexes_of_non_empty_sdrf_columns]

            # Number of bundles processed per study - test mode cut-off
            if mode == 'test' and bundle_cnt >= utils.get_val(config, 'test_max_bundles'): 
                break

        # Now work out which Characteristics should be auto-generated as Factors also
        technology2factors = {}
        # Assumption - in experiments imported from HCA DCC, only one column for a unique characteristic name will be output in the resulting SDRF file
        technology2factor2characteristic_colnum = {}
        for technology in technologies_found:
            technology2factors[technology] = []
            technology2factor2characteristic_colnum[technology] = {}
            for characteristic in characteristic_values:
                if characteristic in technology2sdrf_column_headers[technology] and len(characteristic_values[characteristic]) > 1:
                    factor = re.sub("Characteristics","FactorValue", characteristic)
                    technology2factors[technology].append(factor)
                    technology2sdrf_column_headers[technology].append(factor)
                    # Store index (in each sdrf row) of the characteristic corresponding factor, so that we know where to get the value from
                    # when populating factor values in sdrf later
                    technology2factor2characteristic_colnum[technology][factor] = technology2sdrf_column_headers[technology].index(characteristic)

            # Add Factor for single cell identifier (smart-seq2 experiments only)
            smart_regex = re.compile('smart-.*$')
            if smart_regex.match(technology):
                factor = 'FactorValue[single cell identifier]'
                technology2sdrf_column_headers[technology].append(factor)
                technology2factors[technology].append(factor)
                technology2factor2characteristic_colnum[technology][factor] = technology2sdrf_column_headers[technology].index('Source Name')

        # For each technology, write out the generated SDRF file.
        # N.B. IF the HCA project is multi-technology, append the technology label to the end of the sdrf file name
        multi_technology_hca_project = len(technologies_found) > 1
        for technology in technologies_found:
            sdrf_file_name = "%s.sdrf.txt" % accession
            if multi_technology_hca_project:
                sdrf_file_name = "%s.%s" % (sdrf_file_name, technology)
            with open(os.path.join(data_dir, sdrf_file_name), 'wb') as f:
                csvwriter = csv.writer(f, delimiter = '\t', encoding='utf-8', escapechar='\\', quotechar='', lineterminator='\n', quoting=csv.QUOTE_NONE)
                
                # Remove from technology2sdrf_column_headers[technology] headers of columns that are empty for this technology
                utils.remove_empty_columns(technology2sdrf_column_headers[technology], technology2indexes_of_empty_columns[technology])
                # Expand protocol column headers to account for multiple protocols per protocol_type, if applicable
                expanded_headers = technology2sdrf_column_headers[technology].copy()
                utils.expand_protocol_columns(None, expanded_headers, technology2protocol_type2max_protocol_num_per_sample[technology], logger)

                # Write out sdrf header line
                csvwriter.writerow(expanded_headers)
                
                for row in technology2rows[technology]:
                    # Append to row values for all the auto-generated factors
                    for factor in technology2factors[technology]:
                        row.append(row[technology2factor2characteristic_colnum[technology][factor]])

                    # Remove from row elements in positions corresponding to columns that are empty for this technology
                    utils.remove_empty_columns(row, technology2indexes_of_empty_columns[technology])
                    # Expand protocol values into multiple columns to account for multiple protocols per protocol_type, if applicable
                    utils.expand_protocol_columns(row, technology2sdrf_column_headers[technology], technology2protocol_type2max_protocol_num_per_sample[technology], logger)
                    
                    # Write out sdrf data line
                    csvwriter.writerow(row)

        #################
        ###### IDF ######
        #################
        for technology in technologies_found:
            idf_file_name = "%s.idf.txt" % accession
            if multi_technology_hca_project:
                idf_file_name = "%s.%s" % (idf_file_name, technology)                    
            with open(os.path.join(data_dir, idf_file_name), 'wb') as f:
                csvwriter = csv.writer(f, delimiter = '\t', encoding='utf-8', escapechar='\\', quotechar='', lineterminator='\n', quoting=csv.QUOTE_NONE)
                for line_item in idf_config:
                    magetab_label = line_item[0]
                    hca_path = line_item[1]
                    if isinstance(hca_path, list):
                        if magetab_label in ['Term Source Name','Term Source File']:
                            # Special handling/parsing - hca_path is a list of literal values, rather than locations in HCA json files
                            csvwriter.writerow([magetab_label] + hca_path)
                            continue
                        if hca_path:
                            # Note the assumption that only one project_json object exists per bundle
                            # (c.f. hca_schemas_with_one_json_per_bundle_expected in hca2mtab.yml)
                            json_dict = hca_json_for_bundle[hca_path[0]][0]
                        value = utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label, context)
                        if magetab_label in ['Public Release Date'] and value != utils.get_val(config, 'notfound'):
                            # Special handling/parsing - Public Release date, Comment[HCALastUpdateDate], Comment[HCAReleaseDate]
                            m = re.search(r'^(\d{4}\-\d{2}\-\d{2}).*$', value)
                            if m:
                                value = m.group(1)
                            else:
                                logger.error("Failed to parse date out of: %s" % value)
                                value = ''
                            csvwriter.writerow([magetab_label, value])
                        elif magetab_label in ['Comment[ExpressionAtlasAccession]', 'SDRF File']:
                            # Special handling/parsing - use previously derived accession
                            value = accession
                            if magetab_label == 'SDRF File':
                                # SDRF file name - derive from experiment accession
                                value = re.sub(r'\.idf\.', '.sdrf.', idf_file_name)
                            candidate_acc_regex_obj = re.compile('E-CAND-\d+')                                
                            if magetab_label == 'SDRF File' or (magetab_label == 'Comment[ExpressionAtlasAccession]' and not candidate_acc_regex_obj.match(accession)):
                                csvwriter.writerow([magetab_label, value])
                        elif magetab_label in ['Comment[HCALastUpdateDate]']:
                            csvwriter.writerow([magetab_label, datetime.now().strftime("%Y-%m-%d")])
                        elif magetab_label == 'Comment[SecondaryAccession]':
                            # Special handling - secondary accessions
                            secondary_accessions = OrderedSet([])
                            for label in utils.get_val(config, 'hca_old_secondary_accessions_labels'):
                                hca_project_json = hca_json_for_bundle[utils.get_val(config, 'hca_project')]
                                if label in hca_project_json:
                                    secondary_accessions.add(hca_project_json[label])
                            # For the reason for the loop below see a comment near hca_old_secondary_accessions_labels in hca2mtab.yml
                            for label in utils.get_val(config, 'hca_new_secondary_accessions_labels'):
                                if label in hca_project_json:
                                    for secondary_accession in hca_project_json[label]:
                                        secondary_accessions.add(secondary_accession)
                            # Now append the HCA study uuid
                            secondary_accessions.add(project_uuid)
                            if len(secondary_accessions) > 0:
                                csvwriter.writerow(['Comment[SecondaryAccession]'] + list(secondary_accessions))
                        elif magetab_label in ['Experimental Factor Name','Experimental Factor Type']:
                            # Special handling - populate factors that where auto-generated in SDRF above
                            idf_line = [magetab_label]
                            for factor in technology2factors[technology]:
                                m = re.search(r'\[(.*)\]', factor)
                                if m:
                                    idf_line.append(m.group(1))
                                else:
                                    err_msg = "Failed to extract Factor name from %s" % factor
                                    logger.error(err_msg)
                                    raise utils.HCA2MagetabTranslationError(err_msg)
                            csvwriter.writerow(idf_line)
                        elif isinstance(magetab_label, list):
                            if re.search('Person Last Name', magetab_label[0]):
                                # Special handling/parsing - Contributors
                                contact_rows = OrderedDict()
                                for row_label in magetab_label:
                                    contact_rows[row_label] = []
                                for contact in utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label, context):
                                    contact_name_arr = contact['contact_name'].split(',')
                                    contact_rows['Person Last Name'].append(contact_name_arr[0])
                                    contact_rows['Person First Name'].append(contact_name_arr[-1].lstrip())
                                    if len(contact_name_arr) == 3:
                                        contact_rows['Person Mid Initials'].append(contact_name_arr[1])
                                for contact in utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label, context):
                                    email = utils.get_hca_value(['email'], contact, logger, config, True, magetab_label, context)
                                    contact_rows['Person Email'].append(email if email != utils.get_val(config, 'notfound') else '')
                                    contact_rows['Person Affiliation'].append(contact['institution'])
                                for contact in utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label, context):
                                    address = utils.get_hca_value(['address'], contact, logger, config, True, magetab_label, context)
                                    contact_rows['Person Address'].append(address if address != utils.get_val(config, 'notfound') else '')
                                for key in list(contact_rows.keys()):
                                    csvwriter.writerow([key] + contact_rows[key])
                            elif 'Protocol Name' == magetab_label[0]:
                                # Special handling/parsing - Protocols
                                protocol_rows = OrderedDict()
                                for row_label in magetab_label:
                                    protocol_rows[row_label] = []
                                for protocol_type in technology2protocol_type2protocols[technology].keys():
                                    # Traverse through protocol tuples in alphabetic order - by protocol name
                                    for protocol_tuple in sorted(technology2protocol_type2protocols[technology][protocol_type], key=lambda x: x[0]):
                                        protocol_rows['Protocol Name'].append(protocol_tuple[0])
                                        protocol_rows['Protocol Description'].append(protocol_tuple[1] if protocol_tuple[1] != utils.get_val(config, 'notfound') else '')
                                        protocol_rows['Protocol Type'].append(protocol_tuple[2] if protocol_tuple[2] != utils.get_val(config, 'notfound') else '')
                                for key in list(protocol_rows.keys()):
                                    csvwriter.writerow([key] + protocol_rows[key])
                            elif re.search('Publication Title', magetab_label[0]):
                                if utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label[0], context) == utils.get_val(config, 'notfound'):
                                    # Skip the publications-related idf config
                                    continue
                                # Special handling/parsing - Publications
                                publication_rows = OrderedDict()
                                for row_label in 'Publication Title', 'Publication Author List', 'PubMed ID', 'Publication DOI':
                                    publication_rows[row_label] = []
                                for publication in utils.get_hca_value(hca_path[1:], json_dict, logger, config, True, magetab_label, context):
                                    publication_rows['Publication Title'].append(
                                        utils.get_hca_value(utils.get_val(config, 'hca_publication_title_path'), publication, logger, config, True, magetab_label, context))
                                    publication_rows['Publication Author List'].append(
                                        ', '.join(utils.get_hca_value(utils.get_val(config, 'hca_publication_authors_path'), publication, logger, config, True, magetab_label, context)))
                                    pubmed_id = utils.get_hca_value(utils.get_val(config, 'hca_publication_pmid_path'), publication, logger, config, True, magetab_label, context)
                                    publication_rows['PubMed ID'].append(str(pubmed_id) if str(pubmed_id) != utils.get_val(config, 'notfound') else '')
                                    publication_doi = utils.get_hca_value(utils.get_val(config, 'hca_publication_doi_path'), publication, logger, config, True, magetab_label, context)
                                    publication_rows['Publication DOI'].append(publication_doi if publication_doi != utils.get_val(config, 'notfound') else '')
                                        
                                for key in list(publication_rows.keys()):
                                    csvwriter.writerow([key] + publication_rows[key])
                        else:
                            # magetab_label is not a list or a special case
                            csvwriter.writerow([magetab_label, value])
                            if magetab_label == 'Investigation Title':
                                imported_experiments.append("%s (%s - %d bundles): %s" % (accession, technology, len(hca_json_for_project_uuid.keys()), value))
                    else:
                        # hca_path is not a list
                        csvwriter.writerow(line_item)
        time_end = utils.unix_time_millis(datetime.now())
        duration = (time_end - time_start)/1000/60
        logger.info("Processing HCA study uuid: %s for gxa accession: %s took %d mins" % (project_uuid, accession, duration))
    if imported_experiments and sender and email_recipients:
        utils.email_report("New experiments imported from HCA DCC", '\n'.join(imported_experiments), sender, email_recipients)

if __name__ == '__main__':
    # Capture call arguments
    if len(sys.argv) < 3:
        print('Call arguments needed, e.g. : ')
        print(sys.argv[0] + ' test /magetab/output/path (e8642221-4c2c-4fd7-b926-a68bce363c88,d96c2451-6e22-441f-a3e6-70fd0878bb1b hca2mtab@yourorg.com smith@ebi.ac.uk,jones@ebi.ac.uk')
        sys.exit(1)

    mode = sys.argv[1]
    data_dir = sys.argv[2]
    sender, email_recipients, project_uuids_filter = None, None, None
    if len(sys.argv) > 3 and sys.argv[3].strip() != 'None':
        # A comma-separated (no spaces) list of HCA project uuids to be imported
        project_uuids_filter = sys.argv[3].split(',')
    if len(sys.argv) > 4:
        sender = sys.argv[4]
        email_recipients = sys.argv[5]
    try:
        convert_hca_json_to_magetab(mode, data_dir, project_uuids_filter, True, sender, email_recipients )
    except utils.HCA2MagetabTranslationError as exc:
        utils.email_report("%s error" % process_name, "%s has crashed with the following error: %s" % (process_name, str(exc)), sender, email_recipients)
