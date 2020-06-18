#!/usr/bin/env python2.7

import os
from collections import defaultdict


###################################### USAGE ######################################

### Steps
# 1: download GTF  file (e.g. "Mus_musculus.GRCm38.90.gtf" )
# 2; put file in current working directory
# 3: run this script
# 4: use the output (e.g. Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz) as mapping file

### INPUT GTF files
# HUMAN (current Ensembl): ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
# MOUSE (current Ensembl): ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/

###################################### INTRODUCTION ######################################

### General application and purpose

# Mapping rate
	# In order to maximize mapping rate to/from gene names and Ensembl IDs, we need seperate files.
	# This is because we can always map Ensembl IDs to their gene name, but we cannot uniquely map the gene name to Ensembl ID.
# Gene version numbers
	# Furthermore, ensembl gene names often contain version numbers which would like to eliminate in the final mapping - if possible.
	# So when mapping from to/from Ensembl IDs and gene names, we only want to map to the latest gene version.


### Mapping uniqueness
# 1) Ensembl IDs ALWAYS map to only one gene:  (ens_id, gene_name) is unique with a count of one.
# 2) Gene names does *NOT* map to only one Ensembl ID [that is Multiple Ensembl IDs may map to the same gene name. That is, mapping FROM gene names TO Ensembl IDs is *not* unique.]

### Mix notes
# Ensembl ID and gene_name versions is a MUCH bigger problem for human genes - and less a problem for mouse genes.
# Reference on Ensembl 'stable IDs': https://www.ensembl.org/info/genome/stable_ids/index.html


###################################### OUTPUT FILES ######################################


### file_out_ensembl2gene_name_version
# This file maps FROM ensembl_gene_id TO gene names. 
# Ensembl IDs ALWAYS map to only one gene name
# All lines in the file are unique
# Col1: ensembl_gene_id 					Ensembl gene ID
# Col2:	gene_name_optimal			  		Gene names 'final mapping'. That is, gene name version for ALL genes EXCEPT the latest version of the gene.
# Col3:	gene_name_base				 		Gene name base. The 'base name' of the gene (no version number).
# Col4:	gene_name_versions					Gene name version. Gene name version for ALL genes. (This is the 'raw mapping')
# Col5: flag_latest_version					Boolean flag indicating if the gene_name_version is the latest version.

### file_out_gene_name_version2ensembl
# This file maps FROM gene_name_optimal TO ensembl_gene_id. 
# All gene_name_optimal (Col1) maps to only ensembl_gene_id. 
# [This is not the case without carefully mapping - sometimes does *NOT* map to only one Ensembl ID (e.g. the human genes 'Y_RNA' or '5S_rRNA')]
# All lines in the file are unique.
# Col1: gene_name_optimal					Gene names 'final mapping'. That is, gene name version for ALL genes EXCEPT the latest version of the gene. 
# 											*IMPORTANTLY*, the gene name with the latest version is listed *TWICE*: BOTH with version number and without.
#											The reason for this, to list all possible gene_name_version in this column - in maximize mapping rate.
#											When a user wants to map a given gene WITHOUT a version number, we want to map the gene to the LATEST ensembl_gene_id.
# Col2: ensembl_gene_id						Ensembl gene ID
# Col3: gene_name_base 						Gene name base. The 'base name' of the gene (no version number).
# Col4: gene_name_version					Gene name version. Gene name version for ALL genes. (This is the 'raw mapping')
# Col5: flag_latest_version					Boolean flag indicating if the gene_name_version is the latest version.


### file_out_unique_gene_mapping
# A file with UNIQUE 1-1 mapping of gene names to Ensembl IDs.


### file_out_gene_name_base_info
# Genes that where filtered out because of none-unique mapping
# Col1: gene_name_base 						Gene name base. The 'base name' of the gene (no version number).
# Col2: gene_name_version_latest			The latest version of the gene_name_base
# Col3: number_of_gene_name_versions		The number of different gene_name_version for the gene_name_base
# Col4: gene_name_versions [list]			List of all gene_name_version for the gene_name_base

### file_out_multi_mapping
# Genes that where filtered out because of none-unique mapping
# ... See the file for your self.

### file_out_lines_could_not_map
# Lines that could not be mapped from the GTF file.


###################################### CONSTANTS ######################################
# TODO: change file_mapping to a command line input

#file_mapping = "Homo_sapiens.GRCh38.90.gtf"
file_mapping = "Mus_musculus.GRCm38.90.gtf"
# file_mapping = "gene_mapping.test_files/test_file.Homo_sapiens.GRCh38.90.gtf"

#file_mapping = "../gtf-cellranger/refdata-cellranger-mm10-1.2.0.genes.gtf"


file_out_ensembl2gene_name_version="{}.ensembl2gene_name_version.txt".format(os.path.splitext(file_mapping)[0])
file_out_gene_name_version2ensembl="{}.gene_name_version2ensembl.txt".format(os.path.splitext(file_mapping)[0])
file_out_unique_gene_mapping="{}.unique_gene_mapping.txt".format(os.path.splitext(file_mapping)[0])
file_out_gene_name_base_info="{}.gene_name_base_info.txt".format(os.path.splitext(file_mapping)[0])
file_out_multi_mapping="{}.multi_mapping.txt".format(os.path.splitext(file_mapping)[0])
file_out_lines_could_not_map="{}.lines_could_not_map.txt".format(os.path.splitext(file_mapping)[0])

###################################### Initializations ######################################

dict_lines_could_no_map = {}

ensembl_gene_id2gene_name_version = defaultdict(list)
gene_name_version2ensembl_gene_id = defaultdict(list)

dict_gene_name_base_info = defaultdict(lambda: defaultdict(list)) # a dict of dict of lists
	# dict_gene_name_base_info[GENE_NAME_BASE]: dict. GENE_NAME_BASE is a string of the 'base' gene name (e.g. 'AC233724')
		# dict_gene_name_base_info[GENE_NAME_BASE]['gene_name_version_latest']: list with ONE string element (the gene_name with the latest version, e.g. 'AC233724.19')
		# dict_gene_name_base_info[GENE_NAME_BASE]['gene_name_versions']: list of strings all gene_name versions (e.g. 'AC233724.1', 'AC233724.2', ...)
		# dict_gene_name_base_info[GENE_NAME_BASE]['version_numbers']: list of integers with all gene_name version_numbers (e.g. 1, 2 ,3...). We use these values to determine the latest version of a gene.
		# dict_gene_name_base_info[GENE_NAME_BASE]['ensembl_gene_id']: list of ensembl_gene_id's associated with the GENE_NAME_BASE. This is used to create a 1-1 mapping file


###################################### Read GTF file ######################################
with open(file_mapping, 'r') as fh_in:
	for line in fh_in:
		if line.startswith("#"):
			continue
		line = line.strip()
		fields = line.split("\t")
		attribute_field = fields[8]
		attributes = attribute_field.strip(";").split(";") # we strip the trailing ";" avoid an 'empty' last field
		tmp_entry_attributes_dict = {}
		for attribute in attributes: # PROCESSING: ensembl_gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name_version "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";
			tag_value_pair = attribute.strip(" ").split(" ")
			tag = tag_value_pair[0]
			value = tag_value_pair[1]
			tmp_entry_attributes_dict[tag] = value
		try:
			ensembl_gene_id = tmp_entry_attributes_dict["gene_id"].replace('"', '') # e.g. ENSMUSG00000102693
			gene_name_version = tmp_entry_attributes_dict["gene_name"].replace('"', '') # e.g. 4933401J01Rik
			ensembl_gene_id2gene_name_version[ensembl_gene_id].append(gene_name_version)
			gene_name_version2ensembl_gene_id[gene_name_version].append(ensembl_gene_id)
		except Exception as e:
			dict_lines_could_no_map[line] = "EXCEPTION MESSAGE: {} | EXCEPTION ARGS: {}".format(e.message, e.args) # saving string presentation of Exception.
		

		
###################################### Remove redundant information from GTF file entries ######################################
### Taking unique set of genes AND alphanumeric sorting 
# This is because we will encounter the same information many times as we read the gtf file.
for ensembl_gene_id in ensembl_gene_id2gene_name_version:
	ensembl_gene_id2gene_name_version[ensembl_gene_id] = sorted(list(set(ensembl_gene_id2gene_name_version[ensembl_gene_id])))

for gene_name_version in gene_name_version2ensembl_gene_id:
	gene_name_version2ensembl_gene_id[gene_name_version] = sorted(list(set(gene_name_version2ensembl_gene_id[gene_name_version])))

###################################### Populating dict_gene_name_base_info ######################################

for gene_name_version in gene_name_version2ensembl_gene_id: # the dict 'gene_name_version2ensembl_gene_id' contains ALL seen gene_name_versions
	tmp_split = gene_name_version.split(".") 
	gene_name_base = tmp_split[0] # this will never fail because split() will always return at least one element
	if len(tmp_split) == 1:
		version_number = -1 # no version information attached to gene name. -1 is an appropriate 'default' value, because we want it to be a string.
	elif len(tmp_split) == 2:
		try:
			version_number = int(tmp_split[1])  # save version number as integer
			# mm10 grf contains a gene "H2-M10.5-ps1", which int() conversion will fail for.
			# hg38 gtf does not have issues
		except: # we do nothing but warn the user (and hope these cases are rare).
			print "ERROR in converting version_number to integer: gene_name_version={}".format(gene_name_version)
			print "Will skip this version of the gene. If there are many such problems, consider correcting for this problem in the code."
		
	else:
		raise Exception("gene_name_version split resulted >2 splits: {}".format(tmp_split))
	
	dict_gene_name_base_info[gene_name_base]['gene_name_versions'].append(gene_name_version) 
	dict_gene_name_base_info[gene_name_base]['version_numbers'].append(version_number)
	
	ensembl_gene_id = gene_name_version2ensembl_gene_id[gene_name_version] # this returns a list with one or more elements
	dict_gene_name_base_info[gene_name_base]['ensembl_gene_id'].extend(ensembl_gene_id) # call .extend() and not .append() to avoid a list of lists.
	if version_number >= max(dict_gene_name_base_info[gene_name_base]['version_numbers']): # Update gene_name_version_latest if we encounter a (new) highest version. This also works in the case where it is the first time we see the gene_name_base
			dict_gene_name_base_info[gene_name_base]['gene_name_version_latest'] = [gene_name_version] # we want 'gene_name_version_latest' to always be a one element list



########################################################################################################
###################################### WRITE PRIMARY OUTPUT FILES ######################################
########################################################################################################

### ensembl2gene_name_version
with open(file_out_ensembl2gene_name_version, 'w') as fh_out:

	### Write header
	fh_out.write("{}\t{}\t{}\t{}\t{}\n".format("ensembl_gene_id", "gene_name_optimal", "gene_name_base", "gene_name_version", "flag_latest_version"))

	for ensembl_gene_id in ensembl_gene_id2gene_name_version: # no sorting nessesary
		if len(ensembl_gene_id2gene_name_version[ensembl_gene_id]) == 1: # unique mapping 'filtering step'
			gene_name_version = ensembl_gene_id2gene_name_version[ensembl_gene_id][0] # we know we only have one element
			gene_name_base = gene_name_version.split(".")[0] # this will never fail because split() will always return at least one element
			
			# get the latest version of the gene:
			gene_name_version_latest = dict_gene_name_base_info[gene_name_base]['gene_name_version_latest'][0] # this list always has only one element list 
			
			if gene_name_version == gene_name_version_latest:
				gene_name_optimal = gene_name_base # use base name if the gene is the most recent version
				flag_latest_version = True
			else:
				gene_name_optimal = gene_name_version
				flag_latest_version = False
			
			# write output
			fh_out.write("{}\t{}\t{}\t{}\t{}\n".format(ensembl_gene_id, gene_name_optimal, gene_name_base, gene_name_version, flag_latest_version))


### gene_name_version2ensembl
# The writing of this file is similar to above, except for:
# 1) The dict we loop over (here we loop over gene_name_version2ensembl_gene_id)
# 2) The order of the two first columns
# 3) IMPORANT: That the gene name with the latest version is WRITTEN *TWICE* to the output file: BOTH with version number and without.
with open(file_out_gene_name_version2ensembl, 'w') as fh_out:

	### Write header
	fh_out.write("{}\t{}\t{}\t{}\t{}\n".format("gene_name_optimal", "ensembl_gene_id", "gene_name_base", "gene_name_version", "flag_latest_version"))

	for gene_name_version in gene_name_version2ensembl_gene_id: # no sorting nessesary
		if len(gene_name_version2ensembl_gene_id[gene_name_version]) == 1: # unique mapping 'filtering step'
			ensembl_gene_id = gene_name_version2ensembl_gene_id[gene_name_version][0] # we know we only have one element
			
			gene_name_base = gene_name_version.split(".")[0] # this will never fail because split() will always return at least one element
			
			# get the latest version of the gene:
			gene_name_version_latest = dict_gene_name_base_info[gene_name_base]['gene_name_version_latest'][0] # this list always has only one element list 
			
			if gene_name_version == gene_name_version_latest:
				gene_name_optimal = gene_name_base # use base name if the gene is the most recent version
				flag_latest_version = True
			else:
				gene_name_optimal = gene_name_version
				flag_latest_version = False
			
			# write output - with VERSION number
			fh_out.write("{}\t{}\t{}\t{}\t{}\n".format(gene_name_version, ensembl_gene_id, gene_name_base, gene_name_version, flag_latest_version))
			
			# write output WITHOUT version number
			# but ONLY for the lastest gene_name_version AND if the gene acually has multiple versions:
			# if we did not check for the last condition, we would end up with non-unique values in column1
			if flag_latest_version and (len(dict_gene_name_base_info[gene_name_base]['gene_name_versions']) > 1): 
				fh_out.write("{}\t{}\t{}\t{}\t{}\n".format(gene_name_optimal, ensembl_gene_id, gene_name_base, gene_name_version, flag_latest_version))
		

### unique_gene_mapping
with open(file_out_unique_gene_mapping, 'w') as fh_out:

	### Write header
	fh_out.write("{}\t{}\n".format("ensembl_gene_id", "gene_name_base"))

	for gene_name_base in dict_gene_name_base_info:
		if len(dict_gene_name_base_info[gene_name_base]['ensembl_gene_id']) == 1: # *** unique mapping of gene_name_base ***
			ensembl_gene_id = dict_gene_name_base_info[gene_name_base]['ensembl_gene_id'][0] # we know we only have one element
			fh_out.write("{}\t{}\n".format(ensembl_gene_id, gene_name_base))


########################################################################################################
###################################### WRITE SECONDARY OUTPUT FILES ######################################
########################################################################################################



with open(file_out_gene_name_base_info, 'w') as fh_out:
	### Write header
	fh_out.write("{}\t{}\t{}\t{}\n".format("gene_name_base", "gene_name_version_latest", "number_of_gene_name_versions", "gene_name_versions"))
	for gene_name_base in sorted(dict_gene_name_base_info, key=lambda k:len(dict_gene_name_base_info[k]['gene_name_versions']), reverse=True):
		gene_name_version_latest = dict_gene_name_base_info[gene_name_base]['gene_name_version_latest'][0] # we know it only contains one element
		number_of_gene_name_versions = len(dict_gene_name_base_info[gene_name_base]['gene_name_versions'])
		gene_name_versions = dict_gene_name_base_info[gene_name_base]['gene_name_versions']

		fh_out.write("{}\t{}\t{}\t{}\n".format(gene_name_base, gene_name_version_latest, number_of_gene_name_versions, ",".join(gene_name_versions)))
		
		


with open(file_out_multi_mapping, 'w') as fh_out:
	fh_out.write("################### Ensembl --> Gene Name ###################\n")
	for ensembl_gene_id, dummy in sorted(ensembl_gene_id2gene_name_version.items(), key=lambda (k, v): int(len(v)), reverse=True): # THIS WORKS
		if len(ensembl_gene_id2gene_name_version[ensembl_gene_id]) > 1: # IMPORTANT
			fh_out.write(str(len(ensembl_gene_id2gene_name_version[ensembl_gene_id])) + "\t" + ensembl_gene_id + "\t" + "\t".join(ensembl_gene_id2gene_name_version[ensembl_gene_id]) + "\n" )
	fh_out.write("################### Gene Name --> Ensembl ###################\n")
	for gene_name_version, dummy in sorted(gene_name_version2ensembl_gene_id.items(), key=lambda (k, v): int(len(v)), reverse=True): # THIS WORKS
		if len(gene_name_version2ensembl_gene_id[gene_name_version]) > 1: # IMPORTANT
			fh_out.write(str(len(gene_name_version2ensembl_gene_id[gene_name_version])) + "\t" + gene_name_version + "\t" + "\t".join(gene_name_version2ensembl_gene_id[gene_name_version]) + "\n" )


with open(file_out_lines_could_not_map, 'w') as fh_out:
	for line in dict_lines_could_no_map:
		fh_out.write(line + "\tEXCEPTION:" + dict_lines_could_no_map[line] + "\n" )

print "Number of lines that could not be mapped:", len(dict_lines_could_no_map)

print "DONE"


