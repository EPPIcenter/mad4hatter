# Pipeline Resources

This page explains how to configure computational resources (CPU, memory, time) for pipeline modules. This is only necessary if you are having trouble with the default settings. 

In many cases, if the pipeline is failing, it is due to memory issues on the DADA2 step. In most of these cases, using the config file with the extra memory configuration that is already supplied in the repository is sufficient. For example:

```bash
nextflow run main.nf \
  --readDIR tests/example_data/example_fastq \
  --pools D1,R1,R2 \
  -profile sge,apptainer \
  -c conf/custom.config
```

Read the information below if more advanced optimization is required. 

## Understanding Resource Tiers

The pipeline organizes modules into resource tiers defined in `conf/base.config`:

| Tier | Typical Use | Example Modules |
|:----:|:------------|:----------------|
| **Single** | Lightweight modules (file creation, simple processing) | CREATE_PRIMER_FILES, BUILD_ALLELETABLE |
| **Low** | Moderate processing (filtering, quality reports) | CUTADAPT, FILTER_ASVS, QUALITY_REPORT |
| **Medium** | Moderate-heavy processing (masking, reference creation) | MASK_SEQUENCES, CREATE_REFERENCE_FROM_GENOMES |
| **High** | Heavy processing (DADA2, alignment) | DADA2_ANALYSIS, ALIGN_TO_REFERENCE |


## Customizing Resources

### Creating a Custom Config File

Create or edit `conf/custom.config` to override resource settings:

```groovy
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     EPPIcenter / mad4hatter Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Custom config options for institutions
----------------------------------------------------------------------------------------
*/

process {
    // Example: Customize DADA2_ANALYSIS resources
    withName: 'DADA2_ANALYSIS' {
      time = '600m'      // 10 hours
      cpus = 4           // 4 CPU cores
      penv = 'smp'       // Parallel environment (SGE-specific)
      memory = '16 GB'   // 16 GB RAM
    }
    
    // Example: Customize CUTADAPT resources
    withName: 'CUTADAPT' {
      time = '120m'      // 2 hours
      cpus = 2
      memory = '8 GB'
    }
}
```

### Resource Parameters

| Parameter | Description | Example Values |
|:---------:|:------------|:---------------|
| `time` | Maximum execution time | `'120m'` (2 hours), `'600m'` (10 hours), `'24h'` |
| `cpus` | Number of CPU cores | `2`, `4`, `8` |
| `memory` | RAM allocation | `'8 GB'`, `'16 GB'`, `'32 GB'` |
| `penv` | Parallel environment (SGE-specific) | `'smp'` (usually don't change) |

---

## Threading

Only certain modules benefit from multiple CPUs:

### Multithreaded Modules

- **DADA2_ANALYSIS** - Benefits significantly from multiple cores
- **DADA2_POSTPROC** - Can use multiple cores

### Single-threaded Modules

Most other modules are single-threaded. Increasing `cpus` for these won't speed them up, but may help if multiple samples are processed in parallel.

!!! note "Memory vs. CPUs"
    Increasing `cpus` for multithreaded modules also increases memory usage. Monitor your system to find the right balance.

---

## Executor Settings

### Queue Size

Control how many jobs can run simultaneously:

```groovy
executor {
  $sge {
      queueSize = 1000  // Maximum jobs in queue
  }
}

executor {
  $local {
      queueSize = 1000  // For local execution
  }
}
```

**When to adjust:**
- **Decrease** (`500`) if your scheduler becomes overloaded
- **Increase** (`2000`) for large datasets with many parallel samples

!!! info "Queue Size vs. Threading"
    `queueSize` controls the number of jobs submitted, not threading within a module. This is useful for managing cluster load.

---

## Example: Complete Custom Config

```groovy
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Custom config for large dataset (20+ GB)
----------------------------------------------------------------------------------------
*/

process {
    // High-resource modules
    withName: 'DADA2_ANALYSIS' {
      time = '1200m'    // 20 hours
      cpus = 8
      memory = '32 GB'
    }
    
    withName: 'ALIGN_TO_REFERENCE' {
      time = '600m'
      cpus = 4
      memory = '16 GB'
    }
    
    // Medium-resource modules
    withName: 'MASK_SEQUENCES' {
      time = '300m'
      cpus = 2
      memory = '8 GB'
    }
    
    // Low-resource modules
    withName: 'CUTADAPT' {
      time = '180m'
      cpus = 2
      memory = '4 GB'
    }
}

executor {
  $sge {
      queueSize = 2000  // Allow more parallel jobs
  }
}
```

---

## Tips for Resource Management

1. **Start conservative** - Begin with default values, then adjust based on actual usage
2. **Monitor job logs** - Check if jobs are timing out or running out of memory
3. **Adjust incrementally** - Change one parameter at a time to understand its impact
4. **Consider your cluster** - Different HPC systems have different resource limits

---

## Troubleshooting

### Jobs Timing Out

**Symptom**: Jobs fail with timeout errors

**Solution**: Increase `time` parameter for the failing module

### Out of Memory Errors

**Symptom**: Jobs fail with memory errors

**Solution**: Increase `memory` parameter for the failing module

### Slow Processing

**Symptom**: Pipeline runs very slowly

**Solutions**:

- Increase `cpus` for multithreaded modules (DADA2_ANALYSIS, DADA2_POSTPROC)
- Increase `queueSize` to allow more parallel jobs
- Check if other users are using cluster resources
