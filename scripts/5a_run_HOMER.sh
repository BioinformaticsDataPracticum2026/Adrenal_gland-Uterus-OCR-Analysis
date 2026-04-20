#!/usr/bin/env bash
#SBATCH -J homer_adrenal
#SBATCH -p RM-shared
#SBATCH -A bio230007p
#SBATCH -t 08:00:00
#SBATCH -c 8
#SBATCH --mem=15G
#SBATCH -o logs/homer_%j.out
#SBATCH -e logs/homer_%j.err

set -euo pipefail

# Run HOMER motif enrichment for adrenal promoter/enhancer peak subsets.
#
# Before running this script:
# 1. Install HOMER yourself and make sure findMotifsGenome.pl is available.
# 2. Update the PATH line below to point to your HOMER bin directory.
# 3. Update HUMAN_FA and MOUSE_FA to your local genome FASTA files.
#
# This script uses local genome FASTA files as the HOMER genome argument
# instead of HOMER's built-in genome shortcuts such as hg38 or mm10.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Repository-local input/output directories.
ROOT="${REPO_ROOT}/data/Promoters_and_Enhancers"
SPECIFIC="${ROOT}/Specific"
CONSERVED="${ROOT}/Conserved"
OUTROOT="${REPO_ROOT}/results/HOMER"

# User-configurable local genome FASTA files.
# Set these to the genome FASTA paths available in your environment.
HUMAN_FA="/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/hg38.fa"
MOUSE_FA="/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/mm10.fa"

mkdir -p "${OUTROOT}"
mkdir -p logs

# User-configurable HOMER installation path.
# You need to install HOMER yourself and update this PATH if needed.
export PATH="/jet/home/xli51/repos/HOMER/bin:$PATH"

echo "Host: $(hostname)"
echo "Start: $(date)"
echo "OUTROOT: ${OUTROOT}"

run_homer () {
    local foreground="$1"
    local genome="$2"
    local outdir="$3"
    local background="$4"

    echo
    echo "Running:"
    echo "  FG: ${foreground}"
    echo "  BG: ${background}"
    echo "  GENOME: ${genome}"
    echo "  OUT: ${outdir}"

    mkdir -p "${outdir}"

    # Here the second argument is a local genome FASTA path, not a short
    # HOMER genome keyword.
    findMotifsGenome.pl \
        "${foreground}" \
        "${genome}" \
        "${outdir}" \
        -size given \
        -mask \
        -bg "${background}" \
        -p "${SLURM_CPUS_PER_TASK:-8}"
}

run_homer \
    "${SPECIFIC}/human_adrenal_enhancer_candidate_specific.bed" \
    "${HUMAN_FA}" \
    "${OUTROOT}/human_enhancer_specific_vs_conserved" \
    "${CONSERVED}/human_adrenal_enhancer_with_mouse_ortholog.bed"

run_homer \
    "${SPECIFIC}/human_adrenal_promoter_candidate_specific.bed" \
    "${HUMAN_FA}" \
    "${OUTROOT}/human_promoter_specific_vs_conserved" \
    "${CONSERVED}/human_adrenal_promoter_with_mouse_ortholog.bed"

run_homer \
    "${SPECIFIC}/mouse_adrenal_enhancer_specific.bed" \
    "${MOUSE_FA}" \
    "${OUTROOT}/mouse_enhancer_specific_vs_conserved" \
    "${CONSERVED}/mouse_adrenal_enhancer_conserved.bed"

run_homer \
    "${SPECIFIC}/mouse_adrenal_promoter_specific.bed" \
    "${MOUSE_FA}" \
    "${OUTROOT}/mouse_promoter_specific_vs_conserved" \
    "${CONSERVED}/mouse_adrenal_promoter_conserved.bed"

echo "Finished at $(date)"
