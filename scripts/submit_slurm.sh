#!/usr/bin/env bash
# Submit simulations to an existing Slurm installation; cluster settings are explicit.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${1:-${SCRIPT_DIR}/..}" && pwd)"
RESULTS_ROOT="${2:-${ROOT}/results/mppi-paper-v1}"
CPUS="${MPPI_CPUS:-4}"
B="${MPPI_B:-1000}"
for count in "$CPUS" "$B"; do
  [[ "$count" =~ ^[1-9][0-9]*$ ]] && (( ${#count} <= 9 )) || { echo 'CPUS and B must be positive integers.' >&2; exit 1; }
done
command -v sbatch >/dev/null || { echo 'sbatch is unavailable.' >&2; exit 1; }
command -v Rscript >/dev/null || { echo 'Put Rscript on PATH before submitting.' >&2; exit 1; }
mkdir -p "$ROOT/slurm-logs"
options=(--parsable --nodes=1 --ntasks=1 --cpus-per-task="$CPUS"
         --time="${MPPI_TIME_LIMIT:-2-00:00:00}" --mem="${MPPI_MEMORY:-8G}"
         --job-name="mppi-B${B}" --output="$ROOT/slurm-logs/%x-%j.out" --export=ALL)
[[ -z "${MPPI_PARTITION:-}" ]] || options+=(--partition="$MPPI_PARTITION")
[[ -z "${MPPI_NODE:-}" ]] || options+=(--nodelist="$MPPI_NODE")
printf -v wrap 'export MPPI_AVAILABLE_CORES=%q MPPI_RESERVE_WITHIN_ALLOCATION=0 MPPI_B=%q; exec bash %q %q %q' \
  "$CPUS" "$B" "$ROOT/scripts/run_all_server.sh" "$ROOT" "$RESULTS_ROOT"
sbatch "${options[@]}" --wrap="$wrap"
