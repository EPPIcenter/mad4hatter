#!/usr/bin/env python3
"""
Build a Portable Microhaplotype Object (PMO) from the pipeline outputs and validate it
against the PMO schema using pmotools-python.
"""
import argparse
import re
import sys

import pandas as pd
from pmotools.pmo_builder.merge_to_pmo import merge_to_pmo
from pmotools.pmo_builder.mhap_table_to_pmo import (
    create_minimum_library_specimen_dict_from_mhap_table,
    mhap_table_to_pmo,
)
from pmotools.pmo_builder.panel_information_to_pmo import panel_info_table_to_pmo
from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_writer import PMOWriter

# Must match standardise_sample_name() in workflows/demultiplex_amplicons.nf
SEQUENCER_SUFFIX = re.compile(r"_S\d+(_L\d{3})?$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allele_data", required=True,
                        help="Allele table produced by the pipeline (allele_data.txt)")
    parser.add_argument("--amplicon_info", required=True,
                        help="Amplicon info table used by the pipeline")
    parser.add_argument("--panel_name", required=True,
                        help="Panel name to use if the amplicon info has no pool column")
    parser.add_argument("--specimen_key", default=None,
                        help="Optional TSV with library_sample_name and specimen_name columns")
    parser.add_argument("--bioinformatics_run_name", default=None,
                        help="Name of the bioinformatics run recorded in the PMO")
    parser.add_argument("--output", default="pmo.json.gz",
                        help="Output PMO path (.json or .json.gz)")
    return parser.parse_args()


def fail(message):
    sys.exit(f"ERROR: {message}")


def read_specimen_key(path, library_sample_names):
    """
    Read and check the library_sample_name -> specimen_name key against the samples in the
    allele table. pmotools only reports the first missing sample and silently keeps the last
    entry for duplicated library sample names, so check everything here first.
    """
    key = pd.read_csv(path, sep="\t", dtype=str)
    key.columns = key.columns.str.strip()

    required_columns = ["library_sample_name", "specimen_name"]
    missing_columns = [col for col in required_columns if col not in key.columns]
    if missing_columns:
        fail(f"Specimen key '{path}' is missing required columns: {missing_columns}. "
             f"Columns found: {list(key.columns)}")

    key = key[required_columns].apply(lambda col: col.str.strip())
    empty_rows = key[key.isna().any(axis=1) | (key == "").any(axis=1)]
    if not empty_rows.empty:
        fail(f"Specimen key '{path}' has empty values on rows: "
             f"{[i + 2 for i in empty_rows.index]}")

    duplicated = sorted(key.loc[key["library_sample_name"].duplicated(), "library_sample_name"].unique())
    if duplicated:
        fail(f"Specimen key '{path}' has duplicated library_sample_name values: {duplicated}")

    key_names = set(key["library_sample_name"])
    missing_from_key = sorted(set(library_sample_names) - key_names)
    if missing_from_key:
        message = (f"{len(missing_from_key)} library sample(s) in the allele table are not in the "
                   f"specimen key '{path}': {missing_from_key}")
        if any(SEQUENCER_SUFFIX.sub("", name) in key_names for name in missing_from_key):
            message += ("\nThe specimen key appears to use sample names without sequencer suffixes "
                        "(e.g. _S1_L001). Rerun with `--standardise_sample_names true`.")
        fail(message)

    not_in_data = sorted(key_names - set(library_sample_names))
    if not_in_data:
        print(f"WARNING: {len(not_in_data)} library sample(s) in the specimen key have no "
              f"microhaplotypes in the allele table and will not be in the PMO: {not_in_data}",
              file=sys.stderr)

    return dict(zip(key["library_sample_name"], key["specimen_name"]))


def main():
    args = parse_args()

    allele_data = pd.read_csv(args.allele_data, sep="\t",
                              dtype={"sample_name": str, "target_name": str})
    amplicon_info = pd.read_csv(args.amplicon_info, sep="\t", dtype={"target_name": str})

    # Pools become the reactions of a single panel
    reaction_name_col = None
    panel_name = args.panel_name
    if "pool" in amplicon_info.columns:
        reaction_name_col = "pool"
        pools = sorted({pool.strip() for pools in amplicon_info["pool"].dropna()
                        for pool in str(pools).split(",")})
        panel_name = "_".join(pools)

    panel_target_info = panel_info_table_to_pmo(
        amplicon_info,
        panel_name=panel_name,
        target_name_col="target_name",
        forward_primers_seq_col="fwd_primer",
        reverse_primers_seq_col="rev_primer",
        reaction_name_col=reaction_name_col,
    )

    mhap_info = mhap_table_to_pmo(
        allele_data,
        bioinformatics_run_name=args.bioinformatics_run_name,
        library_sample_name_col="sample_name",
        target_name_col="target_name",
        seq_col="asv",
        reads_col="reads",
        pseudocigar_col="pseudocigar_unmasked",
    )

    specimen_key = None
    if args.specimen_key:
        specimen_key = read_specimen_key(args.specimen_key, allele_data["sample_name"].unique())

    library_specimen_info = create_minimum_library_specimen_dict_from_mhap_table(
        mhap_info["detected_microhaplotypes"],
        panel_name=panel_name,
        library_sample_specimen_key=specimen_key,
    )

    pmo = merge_to_pmo(
        mhap_info=mhap_info,
        panel_target_info=panel_target_info,
        specimen_info=library_specimen_info["specimen_info"],
        library_sample_info=library_specimen_info["library_sample_info"],
    )

    PMOChecker().validate_pmo_json(pmo)
    PMOWriter.write_out_pmo(pmo, args.output, overwrite=True)


if __name__ == "__main__":
    main()
