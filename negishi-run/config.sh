# negishi-run/config.sh: settings shared by every kit script. Sourced, not executed.
# shellcheck shell=bash disable=SC2034  # variables are used by the scripts that source this file
# Edit here, not in the scripts; any value can also be overridden from the environment. Values marked "from ~/rnaseq_env_versions.txt" were
# recorded on Negishi on 2026-10-02; preflight warns if the cluster has changed.

# SLURM settings used by the episodes (keep identical to the lesson)
KIT_ACCOUNT=${KIT_ACCOUNT:-rcac-rnaseq}
KIT_PARTITION=${KIT_PARTITION:-cpu}
KIT_QOS_DEFAULT=${KIT_QOS_DEFAULT:-standby}

# Staged workshop data on Depot (read only for the kit; written only by restage.sh).
# The same path appears in learners/setup.md and .claude/CLAUDE.md; keep them in sync.
STAGED=${STAGED:-/depot/workshop/data/rnaseq-workshop}
# Permission model of the staged copy: group (Depot group-readable) or other (world-readable)
STAGED_PERM=${STAGED_PERM:-group}
# Previous staging location, the default source for restage.sh
STAGED_OLD=${STAGED_OLD:-/scratch/negishi/aseethar/rnaseq-workshop}
STAGED_RESULTS=${STAGED_RESULTS:-/scratch/negishi/aseethar/rnaseq-workshop_results}

# Open OnDemand RStudio (bioconductor) app, R 4.4.0
OOD_SIF=${OOD_SIF:-/depot/itap/aseethar/images/rstudio_bioc_ood_rocky8_r4.4.0_s2025.05.0-496_tex.sif}
OOD_HOST_LIB=${OOD_HOST_LIB:-/apps/biocontainers/extras/r-ood/r4.4.0_s2025.05.0-496}
OOD_HOST_LIB_MOUNT=${OOD_HOST_LIB_MOUNT:-/opt/R/host-site-library}

# Where test runs live. Must be under /scratch/negishi/ (the R episodes build their
# working directory as /scratch/negishi/$USER/rnaseq-workshop; see docs/ADAPTATIONS.md).
KIT_SCRATCH_BASE=${KIT_SCRATCH_BASE:-/scratch/negishi/$(id -un)/rnaseq-kit}

# Expected module defaults (from ~/rnaseq_env_versions.txt). The episodes load the
# default version (module load star), so preflight checks that the default is this one.
EXPECTED_MODULES="fastqc/0.12.1 multiqc/1.23 fastp/0.23.2 star/2.7.11b subread/2.0.1 kallisto/0.48.0 salmon/1.10.1 sra-tools/2.11.0-pl5262 seqtk/1.4"
# r-rnaseq (episode 04b tximport) was not recorded; preflight reports whatever loads.
EXPECTED_R_VERSION=4.4.0
EXPECTED_BIOC_VERSION=3.20

# Output of build_kit.py (do not edit by hand)
KIT_GEN_DIRNAME=generated
