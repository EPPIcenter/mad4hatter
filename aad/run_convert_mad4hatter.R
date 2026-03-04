
# WRITE PATH TO DIRECTORY TO BE TRANSFORMED
input_dir   = "PATH/TO/DIR"

# WRITE PATH TO NEW DIRECTORY TO BE CREATED
output_dir  = "PATH/TO/NEW_DIR"

# SET PATH TO PIPELINE
pipeline_path = "PATH/TO/mad4hatter"
# target_id_conversion_table.tsv WITH ALL TARGETS (PROVIDED)
locus_lookup = file.path(pipeline_path,"aad/target_id_conversion_table.tsv")
# references.fasta WITH ALL TARGETS (PROVIDED)
references = file.path(pipeline_path,"aad/references.fasta")
# amplicon_info.tsv WITH ALL TARGETS (PROVIDED)
amplicon_info = file.path(pipeline_path,"aad/amplicon_info.tsv")
# principal_resistance_marker_info_table.tsv, USE THE ONE THAT CONTAINS ALL MARKERS AS IN fix_principal_list BRANCH
resmarker_info = file.path(pipeline_path,"panel_information/principal_resistance_marker_info_table.tsv")


source(file.path(pipeline_path,"aad/convert_mad4hatter_functions.R"))

# RUN
convert_mad4hatter(
  input_dir       = input_dir,
  output_dir      = output_dir,
  locus_lookup    = locus_lookup,
  amplicon_info   = amplicon_info,
  references      = references,
  resmarker_info  = resmarker_info,
  #release_version = OPTIONAL
  )
