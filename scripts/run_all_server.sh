#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${1:-${SCRIPT_DIR}/..}" && pwd)"
RESULTS_ROOT="${2:-${ROOT}/results/mppi-paper-v1}"
command -v Rscript >/dev/null || { echo 'Rscript is unavailable.' >&2; exit 1; }
command -v nproc >/dev/null || { echo 'This launcher requires Linux nproc.' >&2; exit 1; }
VISIBLE_CORES="$(nproc)"
TOTAL_CORES="${MPPI_AVAILABLE_CORES:-${SLURM_CPUS_PER_TASK:-}}"
[[ -n "$TOTAL_CORES" ]] || { echo "Set MPPI_AVAILABLE_CORES to the cores allowed for this job, or run in a Slurm allocation. nproc alone does not establish an allocation." >&2; exit 1; }
for count in "$VISIBLE_CORES" "$TOTAL_CORES"; do
  [[ "$count" =~ ^[1-9][0-9]*$ ]] && (( ${#count} <= 6 )) || { echo 'Invalid CPU allocation.' >&2; exit 1; }
done
(( TOTAL_CORES <= VISIBLE_CORES )) || { echo 'Requested allocation exceeds process-visible CPUs.' >&2; exit 1; }
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  [[ "$SLURM_CPUS_PER_TASK" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid Slurm CPU allocation.' >&2; exit 1; }
  (( TOTAL_CORES <= SLURM_CPUS_PER_TASK )) || { echo 'Requested allocation exceeds Slurm allocation.' >&2; exit 1; }
fi
RESERVE_WITHIN_ALLOCATION="${MPPI_RESERVE_WITHIN_ALLOCATION:-5}"
[[ "$RESERVE_WITHIN_ALLOCATION" =~ ^[0-9]+$ ]] && (( ${#RESERVE_WITHIN_ALLOCATION} <= 6 )) || {
  echo 'MPPI_RESERVE_WITHIN_ALLOCATION must be a nonnegative integer.' >&2
  exit 1
}
(( TOTAL_CORES > RESERVE_WITHIN_ALLOCATION )) || {
  echo 'Allocated cores must exceed MPPI_RESERVE_WITHIN_ALLOCATION.' >&2
  exit 1
}
B="${MPPI_B:-1000}"
[[ "$B" =~ ^[1-9][0-9]*$ ]] && (( ${#B} <= 9 )) || { echo 'MPPI_B must be a positive integer.' >&2; exit 1; }
WORKERS=$((TOTAL_CORES - RESERVE_WITHIN_ALLOCATION))
PARTIAL=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

mkdir -p "$RESULTS_ROOT"
RESULTS_ROOT="$(cd -- "$RESULTS_ROOT" && pwd)"
LOCK="${RESULTS_ROOT}/.launcher_lock"
mkdir "$LOCK" 2>/dev/null || { echo "Another launcher owns $LOCK; inspect its owner.txt before removing a stale lock." >&2; exit 1; }
printf 'host=%s\npid=%s\nstarted=%s\n' "$(hostname)" "$$" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" > "$LOCK/owner.txt"
RUN_ID="$(date -u '+%Y%m%dT%H%M%SZ')-$$"
LOG_DIR="${RESULTS_ROOT}/logs/${RUN_ID}"
cleanup() {
  RC=$?
  trap - EXIT
  if [[ -d "$LOG_DIR" ]]; then printf 'exit_code=%s finished=%s\n' "$RC" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$LOG_DIR/run_environment.log"; fi
  rm -f "$LOCK/owner.txt"
  rmdir "$LOCK" 2>/dev/null || true
  exit "$RC"
}
trap cleanup EXIT
mkdir -p "$LOG_DIR"
printf 'host=%s visible_cores=%s allocated_cores=%s reserve_within_allocation=%s workers=%s B=%s\n' \
  "$(hostname)" "$VISIBLE_CORES" "$TOTAL_CORES" "$RESERVE_WITHIN_ALLOCATION" "$WORKERS" "$B" |
  tee "$LOG_DIR/run_environment.log"
printf 'Logs: %s\n' "$LOG_DIR"
# Resource observations are diagnostics, not a promise of exclusive idle CPUs or memory.
if [[ -r /proc/meminfo ]]; then awk '/MemTotal|MemAvailable/ {print}' /proc/meminfo >> "$LOG_DIR/run_environment.log"; fi
if [[ -r /sys/fs/cgroup/memory.max ]]; then printf 'cgroup_memory_max=' >> "$LOG_DIR/run_environment.log"; cat /sys/fs/cgroup/memory.max >> "$LOG_DIR/run_environment.log"; fi
df -Pk "$RESULTS_ROOT" >> "$LOG_DIR/run_environment.log"
run_step() {
  NAME="$1"; shift
  printf '%s starting %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$NAME"
  if Rscript --vanilla "$@" > "$LOG_DIR/${NAME}.log" 2>&1; then
    return 0
  else
    RC=$?
    printf '%s exited=%s log=%s\n' "$NAME" "$RC" "$LOG_DIR/${NAME}.log" >&2
    if (( RC == 2 )); then PARTIAL=1; return 0; fi
    exit "$RC"
  fi
}
# Preflight failures are fatal even if an external Rscript wrapper uses exit status 2.
if ! Rscript --vanilla "$ROOT/scripts/preflight.R" "--root=$ROOT" "--results=$RESULTS_ROOT" "--B=$B" > "$LOG_DIR/preflight.log" 2>&1; then
  echo "Preflight failed; see $LOG_DIR/preflight.log" >&2; exit 1
fi
for scenario in $(seq 1 12); do
  run_step "scenario_${scenario}" "$ROOT/scripts/run_simulations.R" \
    "--root=$ROOT" "--scenario=$scenario" "--from=1" "--to=$B" "--workers=$WORKERS" "--B=$B" "--results=$RESULTS_ROOT"
done
run_step summary "$ROOT/scripts/summarize_results.R" "--root=$ROOT" "--results=$RESULTS_ROOT"
run_step plots "$ROOT/scripts/plot_results.R" "--root=$ROOT" "--results=$RESULTS_ROOT"
printf 'Finished partial_failures=%s\n' "$PARTIAL" | tee -a "$LOG_DIR/run_environment.log"
if (( PARTIAL )); then exit 2; fi
