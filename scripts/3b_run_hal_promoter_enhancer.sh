#!/usr/bin/env bash

set -euo pipefail

# Submit HALPER liftover jobs for adrenal promoter/enhancer peak sets.
#
# Before running this script:
# 1. Install HAL / halLiftover on your system or cluster environment.
# 2. Install HALPER (for example, the halLiftover-postprocessing repository).
# 3. Update HALPER_SCRIPT below to point to your local
#    halper_map_peak_orthologs.sh script.
# 4. Update HAL below to point to your local multi-species .hal alignment file.
#
# This script uses repository-relative input/output paths for peak files, but
# HAL / HALPER themselves are external dependencies and must be configured by
# the user.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

mkdir -p logs

# User-configurable external dependency paths.
# You need to install HALPER yourself and set this to the correct script path.
HALPER_SCRIPT=/jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh

# You need to provide your own HAL alignment file and set this path manually.
HAL=/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal

# Repository-local input/output directories.
DATADIR="$REPO_ROOT/data/Promoters_and_Enhancers"
OUTDIR=$DATADIR

########################################
# HUMAN → MOUSE
########################################

# enhancer
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/human_adrenal_enhancer_peaks.bed \
  -o ${OUTDIR}/human_adrenal_enhancer_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

# promoter
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/human_adrenal_promoter_peaks.bed \
  -o ${OUTDIR}/human_adrenal_promoter_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

########################################
# MOUSE → HUMAN
########################################

# enhancer
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/mouse_adrenal_enhancer_peaks.bed \
  -o ${OUTDIR}/mouse_adrenal_enhancer_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

# promoter
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/mouse_adrenal_promoter_peaks.bed \
  -o ${OUTDIR}/mouse_adrenal_promoter_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

echo "All jobs submitted."
