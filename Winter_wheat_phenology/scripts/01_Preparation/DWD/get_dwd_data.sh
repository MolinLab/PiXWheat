#!/usr/bin/env bash
# Download the HYRAS variables used for German phenology predictors.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PHENOLOGY_PROJECT_DIR:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"

BASE_URL="https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de"

OUTDIR="${PROJECT_DIR}/Data/DWD"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"

cd "${OUTDIR}"

LOGFILE="${LOGDIR}/download_$(date +%Y%m%d_%H%M%S).log"

# Download years and concurrency (as used in the analysis).

START_YEAR=1978
END_YEAR=2026

PARALLEL=5

VARIABLES=(
    "tasmax"
    "tasmin"
    "pr"
    "rsds"
)

log() {
    echo "[$(date '+%F %T')] $*" | tee -a "${LOGFILE}"
}

download_file() {

    local VAR="$1"
    local YEAR="$2"

    local SUB=""
    local PRE=""
    local VER=""

    case "${VAR}" in

        tasmax)
            SUB="air_temperature_max"
            PRE="tasmax_hyras_1"
            VER="v6-1"
            ;;

        tasmin)
            SUB="air_temperature_min"
            PRE="tasmin_hyras_1"
            VER="v6-1"
            ;;

        pr)
            SUB="precipitation"
            PRE="pr_hyras_1"
            VER="v6-1"
            ;;

        rsds)
            SUB="radiation_global"
            PRE="rsds_hyras_5"
            VER="v4-0"
            ;;

        *)
            log "Unknown variable: ${VAR}"
            return 1
            ;;

    esac

    local FILE="${PRE}_${YEAR}_${VER}_de.nc"

    local URL="${BASE_URL}/${SUB}/${FILE}"

    local TARGET_DIR="${OUTDIR}/${VAR}"

    mkdir -p "${TARGET_DIR}"

    log "Downloading ${VAR} ${YEAR}"
    log "URL: ${URL}"

    if wget \
        --continue \
        --tries=20 \
        --retry-connrefused \
        --waitretry=10 \
        --read-timeout=60 \
        --timeout=60 \
        --random-wait \
        --timestamping \
        --no-verbose \
        --directory-prefix="${TARGET_DIR}" \
        "${URL}" \
        >> "${LOGFILE}" 2>&1
    then
        log "Finished ${VAR} ${YEAR}"
    else
        log "FAILED ${VAR} ${YEAR}"
        return 1
    fi
}

export -f download_file
export -f log

export BASE_URL
export OUTDIR
export LOGFILE


TASKS=$(mktemp)

trap 'rm -f "${TASKS}"' EXIT

for VAR in "${VARIABLES[@]}"; do
    for YEAR in $(seq "${START_YEAR}" "${END_YEAR}"); do
        echo "${VAR} ${YEAR}" >> "${TASKS}"
    done
done

log "Starting HYRAS downloads"

xargs -n 2 -P "${PARALLEL}" bash -c '
download_file "$1" "$2"
' _ < "${TASKS}"

log "All downloads complete"
