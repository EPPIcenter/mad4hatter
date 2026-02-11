# Using Bespoke Pools/Panels

This guide explains how to configure and use custom pools/panels that aren't included in the default pipeline configuration. If you think that the panel you are working with would be useful for others, then please [raise an issue](https://github.com/EPPIcenter/mad4hatter/issues) and label it as a feature request, or submit a pull request if you have already configured it for long-term use. 

## Overview

The pipeline uses **pools** to organize targets. Each pool requires two specific configuration files:

- `amplicon_info.tsv` - Defines amplicon locations, primers, and target information
- `targeted_reference.fasta` - Reference sequences for each amplicon. This can be generated within the pipeline by using the `--genome` flag with a full genome reference that covers all targets.

By default, the pipeline includes pre-configured pools (e.g., D1, R1, R2). See [Pre-Configured Pools](execution.md#pre-configured-pools) for a full list. This page explains how to add your own custom pools. 

!!! note "Pools vs. Panels"
    The term `pools` is used because MAD4HATTER is a modular panel that can be run by combining different pools depending on the use case. If your assay uses a single pool design, you can think of a pool as a panel and use one single pool to define the complete panel. 

---

## When to Use Custom Pools

Use custom pools when:

- You have a custom amplicon panel not included in the default configuration
- You want to add new amplicons to an existing panel
- You're developing a new panel design

---

## Setting Up Custom Pools

There are two options for setting up a custom pool. **Option 1 (Using Command-Line Parameters)** is recommended for one-off or short-term analysis, such as testing a new panel design. **Option 2 (Adding to the Pipeline Configuration)** is recommended for long-term use and sharing functionality between users. If you think this panel would be useful for others, then please follow Option 2 or reach out for help. 
### Option 1: Using Command-Line Parameters

#### Step 1: Prepare Your Amplicon Info File

Create a tab-separated file describing the targets and primers within your panel/pool with the following columns.

!!! tip "Example File"
    See an example amplicon info file: [D1.1_amplicon_info.tsv](https://github.com/EPPIcenter/mad4hatter/blob/develop/panel_information/D1.1/D1.1_amplicon_info.tsv) 

| Column | Description | Example |
|:------:|:------------|:--------|
| `target_id` | Unique identifier for the target | `Pf3D7_01_v3-insertStart-insertEnd` |
| `chrom` | Chromosome name (must match genome/reference) | `Pf3D7_01_v3` |
| `insert_start` | Start position of the amplicon insert | `1000` |
| `insert_end` | End position of the amplicon insert | `1200` |
| `fwd_primer` | Forward primer sequence | `ATCGATCGATCG` |
| `rev_primer` | Reverse primer sequence | `GCTAGCTAGCTA` |
| `pool` | Pool (or panel) name | `MyCustomPool` |

!!! warning "Coordinates"
    - All coordinates should be 0-based (e.g., `insert_start`)
    - If using the resmarker module, make sure the coordinates are relative to the full 3D7 reference.

#### Step 2: Generate a Reference

You can either create a targeted reference that includes only the reference sequences for each target, or supply a whole genome reference that covers all targets in your amplicon info file. If using a genome, the pipeline will use the amplicon info to create a targeted reference automatically.

**Option A: Using a Genome Reference**

- Ensure the genome covers all targets defined in the amplicon info file. If your panel includes targets for multiple species, you may need to combine genomes.
- Ensure the chromosome names in the genome match the `chrom` column in your amplicon info file.
- Supply to the pipeline using the `--genome` flag.

**Option B: Using a Targeted Reference File**

- FASTA file with reference sequences for each amplicon
- Header format: `>target_id` (must match `target_id` in amplicon_info.tsv)
- Supply to the pipeline using the `--refseq_fasta` flag.

!!! tip "Example File"
    See an example targeted reference file: [D1.1_reference.fasta](https://github.com/EPPIcenter/mad4hatter/blob/develop/panel_information/D1.1/D1.1_reference.fasta) 

#### Step 3: Run the Pipeline

**Example: Using a genome reference**
```bash
nextflow run main.nf \
  --readDIR /path/to/data \
  --pools MyCustomPool \
  --amplicon_info /path/to/MyCustomPool_amplicon_info.tsv \
  --genome /path/to/full_genome.fasta \
  -profile docker
```

**Example: Using a targeted reference**
```bash
nextflow run main.nf \
  --readDIR /path/to/data \
  --pools MyCustomPool \
  --amplicon_info /path/to/MyCustomPool_amplicon_info.tsv \
  --refseq_fasta /path/to/MyCustomPool_reference.fasta \
  -profile docker
```

### Option 2: Adding to Pipeline Configuration

#### Step 1: Prepare Your Files

Two files are required to configure the pipeline: an amplicon info file and a targeted reference.

Create a tab-separated file called `<your_pool_name>_amplicon_info.tsv` describing the targets and primers within your panel/pool with the following columns.

!!! tip "Example Files"
    See example files from the D1.1 pool:
    - [D1.1_amplicon_info.tsv](https://github.com/EPPIcenter/mad4hatter/blob/develop/panel_information/D1.1/D1.1_amplicon_info.tsv) - Example amplicon info file
    - [D1.1_reference.fasta](https://github.com/EPPIcenter/mad4hatter/blob/develop/panel_information/D1.1/D1.1_reference.fasta) - Example targeted reference file 

| Column | Description | Example |
|:------:|:------------|:--------|
| `target_id` | Unique identifier for the target | `Pf3D7_01_v3-insertStart-insertEnd` |
| `chrom` | Chromosome name (must match genome/reference) | `Pf3D7_01_v3` |
| `insert_start` | Start position of the amplicon insert | `1000` |
| `insert_end` | End position of the amplicon insert | `1200` |
| `fwd_primer` | Forward primer sequence | `ATCGATCGATCG` |
| `rev_primer` | Reverse primer sequence | `GCTAGCTAGCTA` |

!!! warning "Coordinates"
    - All coordinates should be 0-based (e.g., `insert_start`)
    - If using the resmarker module, make sure the coordinates are relative to the full 3D7 reference.

Create a targeted reference file called `<your_pool_name>_reference.fasta`:

- FASTA file with reference sequences for each amplicon
- Header format: `>target_id` (must match `target_id` in amplicon_info.tsv)

#### Step 2: Add Files to Configuration

Create a folder named after your pool/panel under `panel_information/` and add the files you created in Step 1. For example:

```
panel_information/
└── MyCustomPool/
    ├── MyCustomPool_amplicon_info.tsv
    └── MyCustomPool_reference.fasta
```

Add your new pool to `conf/panel.config`:

```groovy
params {
    pool_options = [
        'MyCustomPool': [
            amplicon_info_path: 'panel_information/MyCustomPool/MyCustomPool_amplicon_info.tsv',
            targeted_reference_path: 'panel_information/MyCustomPool/MyCustomPool_reference.fasta'
        ]
    ]
}
```

Then run with:

```bash
nextflow run main.nf \
  --readDIR /path/to/data \
  --pools MyCustomPool \
  -c conf/custom.config \
  -profile docker
```

---

## Troubleshooting

### Pool Not Found Error

If you see an error like "The following pools were not found in configuration":

1. **Check pool name spelling** - Pool names are case-sensitive
2. **Verify config file** - Ensure your config file is loaded with `-c conf/custom.config`
3. **Check file paths** - Verify paths in your config are correct (relative to project root or absolute)

### File Path Issues

- Use **relative paths** from the project root, or
- Use **absolute paths** if files are outside the project directory
- Ensure file paths in config match actual file locations

### Amplicon Info Format Errors

- Verify all required columns are present
- Check for tab-separated format (not spaces)
- Ensure `target_id` values match between amplicon_info.tsv and reference.fasta headers


## Best Practices

1. **Use descriptive pool names** - Avoid spaces and special characters. Use the format `chrom-insert_start-insert_end` for target IDs to conform with MAD4HATTER naming conventions.
2. **Keep files organized** - Use consistent naming: `PoolName_amplicon_info.tsv` and `PoolName_reference.fasta`
3. **Validate files first** - Check your `amplicon_info.tsv` format before running
4. **Test with QC workflow** - Run `--workflow_name qc` first to verify pool configuration
5. **Document your pools** - Keep notes on pool design and amplicon targets

---

## Getting Help

- **Validation errors**: Check file formats and paths
- **Pool configuration**: Verify your config file syntax
- **File format questions**: See example files in the repository's `panel_information/` directory
- **Issues or questions**: Report issues on [GitHub](https://github.com/EPPIcenter/mad4hatter/issues)


