#!/bin/bash
#SBATCH --job-name=pilotsimu2
#SBATCH --output=logs/pilotsimu2_%x_%A_%a.out
#SBATCH --error=logs/pilotsimu2_%x_%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=12
#SBATCH --array=1-12

# Usage:
#   sbatch submit_simulation.sbatch
#   MODE=cv_only BANDWIDTH=all H_FIXED_VALUES=0.10,0.15,0.20 sbatch --array=1-6 submit_simulation.sbatch
#
# MODE=pilot_fixed runs 3 settings x 2 g_types x 2 n.tr values.
# MODE=cv_only runs 3 settings x 2 g_types with n.tr fixed at 2500.

module --ignore_cache load R

MODE=${MODE:-pilot_fixed}
N=${N:-2500}
NSIM=${NSIM:-500}
NSIM_PILOT=${NSIM_PILOT:-50}
BANDWIDTH=${BANDWIDTH:-all}
H_FIXED=${H_FIXED:-}
H_FIXED_VALUES=${H_FIXED_VALUES:-}

mkdir -p logs

GTYPE_VALUES=("null" "linear")
SETTINGS=(1 2 3)
TASK_ID=${SLURM_ARRAY_TASK_ID}
IDX=$((TASK_ID - 1))

if [ "${MODE}" = "pilot_fixed" ]; then
  NTR_VALUES=(2500 5000)
  OUTDIR=${OUTDIR:-results}

  SETTING_IDX=$((IDX % 3))
  NTR_IDX=$(((IDX / 3) % 2))
  GTYPE_IDX=$((IDX / 6))

  if [ "${GTYPE_IDX}" -ge 2 ]; then
    echo "Invalid task ${TASK_ID} for MODE=pilot_fixed. Use --array=1-12."
    exit 1
  fi

  SETTING=${SETTINGS[$SETTING_IDX]}
  NTR=${NTR_VALUES[$NTR_IDX]}
  GTYPE=${GTYPE_VALUES[$GTYPE_IDX]}

  mkdir -p "${OUTDIR}"

  echo "mode=${MODE}"
  echo "task=${TASK_ID}"
  echo "setting=${SETTING}"
  echo "g_type=${GTYPE}"
  echo "n=${N}"
  echo "n.tr=${NTR}"
  echo "nsim_pilot=${NSIM_PILOT}"
  echo "nsim=${NSIM}"
  echo "outdir=${OUTDIR}"
  echo "bandwidth=pilot_fixed"

  Rscript run_simulation.R \
    --mode pilot_fixed \
    --setting "${SETTING}" \
    --g_type "${GTYPE}" \
    --n "${N}" \
    --n_tr "${NTR}" \
    --nsim_pilot "${NSIM_PILOT}" \
    --nsim "${NSIM}" \
    --bandwidth pilot_fixed \
    --outdir "${OUTDIR}"

elif [ "${MODE}" = "cv_only" ]; then
  OUTDIR=${OUTDIR:-results_cv_ntr2500}
  NTR=${NTR:-2500}

  GTYPE_IDX=$((IDX / 3))
  SETTING_IDX=$((IDX % 3))

  if [ "${GTYPE_IDX}" -ge 2 ]; then
    echo "Invalid task ${TASK_ID} for MODE=cv_only. Use --array=1-6."
    exit 1
  fi

  SETTING=${SETTINGS[$SETTING_IDX]}
  GTYPE=${GTYPE_VALUES[$GTYPE_IDX]}

  mkdir -p "${OUTDIR}"

  echo "mode=${MODE}"
  echo "task=${TASK_ID}"
  echo "setting=${SETTING}"
  echo "g_type=${GTYPE}"
  echo "n=${N}"
  echo "n.tr=${NTR}"
  echo "nsim=${NSIM}"
  echo "outdir=${OUTDIR}"
  echo "bandwidth=${BANDWIDTH}"
  echo "h_fixed=${H_FIXED}"
  echo "h_fixed_values=${H_FIXED_VALUES}"

  CMD=(Rscript run_simulation.R
    --mode cv_only \
    --setting "${SETTING}" \
    --g_type "${GTYPE}" \
    --n "${N}" \
    --n_tr "${NTR}" \
    --nsim "${NSIM}" \
    --bandwidth "${BANDWIDTH}" \
    --outdir "${OUTDIR}")

  if [ -n "${H_FIXED}" ]; then
    CMD+=(--h_fixed "${H_FIXED}")
  fi
  if [ -n "${H_FIXED_VALUES}" ]; then
    CMD+=(--h_fixed_values "${H_FIXED_VALUES}")
  fi

  "${CMD[@]}"

else
  echo "Unknown MODE=${MODE}. Use MODE=pilot_fixed or MODE=cv_only."
  exit 1
fi
