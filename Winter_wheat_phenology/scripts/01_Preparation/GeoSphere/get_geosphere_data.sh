#!/usr/bin/env bash
# Download SPARTACUS daily weather for Austrian phenology predictors.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PHENOLOGY_PROJECT_DIR:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"

BASE_URL="https://public.hub.geosphere.at/datahub/resources/spartacus-v2-1d-1km/filelisting"

OUTDIR="${PROJECT_DIR}/Data/GeoSphere"

LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"

LOGFILE="${LOGDIR}/download_$(date +%Y%m%d_%H%M%S).log"

# Variables:
# TN = minimum temperature
# TX = maximum temperature
# RR = precipitation
# SA = sunshine duration

VARIABLES=("TN" "TX" "RR" "SA")

START_YEAR=1988
END_YEAR=2026

PARALLEL=5

log() {
    echo "[$(date '+%F %T')] $*" | tee -a "${LOGFILE}"
}

download_file() {

    local VAR="$1"
    local YEAR="$2"

    local URL="${BASE_URL}/${VAR}/SPARTACUS2-DAILY_${VAR}_${YEAR}.nc"

    local TARGET_DIR="${OUTDIR}/${VAR}"

    mkdir -p "${TARGET_DIR}"

    log "Downloading ${VAR} ${YEAR}"

    wget \
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

    log "Finished ${VAR} ${YEAR}"
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

xargs -n 2 -P "${PARALLEL}" bash -c '
download_file "$0" "$1"
' < "${TASKS}"

log "All downloads complete."
