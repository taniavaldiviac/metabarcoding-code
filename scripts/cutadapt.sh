#!/usr/bin/env bash
set -euo pipefail

# Recorte de primers paired-end con cutadapt

# Configuración de rutas
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RAW_DIR="${PROJECT_ROOT}/raw_fastqs"
OUT_DIR="${PROJECT_ROOT}/for_dada2"
REPORT_DIR="${PROJECT_ROOT}/cutadapt_reports"
mkdir -p "${OUT_DIR}" "${REPORT_DIR}"

# Primers (override con FWD=... REV=...)
FWD="${FWD:-ACACCGCCCGTCACTCT}"
REV="${REV:-CTTCCGGTACACTTACCATG}"
revcomp(){ echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev; }
RC="$(revcomp "${REV}")"
ADAPTER_R2="${REV}"   # R2 comienza con REV (no RC)

# Parámetros (override con THREADS= MIN_LEN= ERROR_RATE= DISCARD_UNTRIMMED=0/1)
THREADS="${THREADS:-4}"
ERROR_RATE="${ERROR_RATE:-0.10}"
MIN_LEN="${MIN_LEN:-40}"
DISCARD_UNTRIMMED="${DISCARD_UNTRIMMED:-0}"

# Localizar cutadapt
if command -v cutadapt >/dev/null 2>&1; then
  CUTADAPT_BIN="cutadapt"
elif python3 -m cutadapt --version >/dev/null 2>&1; then
  CUTADAPT_BIN="python3 -m cutadapt"
else
  echo "ERROR: cutadapt no encontrado"; exit 1
fi

# Verificaciones
[[ -d "${RAW_DIR}" ]] || { echo "ERROR: falta ${RAW_DIR}"; exit 1; }
shopt -s nullglob
R1_LIST=( "${RAW_DIR}"/*_R1_001.fastq )
(( ${#R1_LIST[@]} > 0 )) || { echo "ERROR: no hay *_R1_001.fastq"; exit 1; }

# Resumen (cutadapt 5.x)
SUMMARY_TXT="${REPORT_DIR}/overall_report.txt"
SUMMARY_TSV="${REPORT_DIR}/overall_summary.tsv"
: > "${SUMMARY_TXT}"
echo -e "sample\tpairs_total\tread1_adapter\tread2_adapter\tpairs_written\tpct_written" > "${SUMMARY_TSV}"

echo "[INFO] FWD=${FWD} REV=${REV} RC=${RC} THREADS=${THREADS} MIN_LEN=${MIN_LEN} DISCARD_UNTRIMMED=${DISCARD_UNTRIMMED}"

for R1 in "${R1_LIST[@]}"; do
  base_R1="$(basename "${R1}")"
  R2="${R1/_R1_001.fastq/_R2_001.fastq}"
  [[ -f "${R2}" ]] || { echo "[WARN] Falta par R2 para ${base_R1}, omite"; continue; }
  sample="${base_R1/_R1_001.fastq/}"

  out_R1="${OUT_DIR}/${base_R1}"
  out_R2="${OUT_DIR}/$(basename "${R2}")"
  report="${REPORT_DIR}/${sample}_cutadapt.txt"

  # Construir argumentos (seguro con set -u)
  args=( -j "${THREADS}" -g "^${FWD}" -G "^${ADAPTER_R2}" -e "${ERROR_RATE}" -m "${MIN_LEN}" )
  if [[ "${DISCARD_UNTRIMMED}" == "1" ]]; then
    args+=( --discard-untrimmed )
  fi
  args+=( -o "${out_R1}" -p "${out_R2}" "${R1}" "${R2}" )

  echo "[INFO] Procesando ${sample}"
  "${CUTADAPT_BIN}" "${args[@]}" > "${report}"

  {
    echo "=== ${sample} ==="
    grep -E "^Total read pairs processed:|^  Read 1 with adapter:|^  Read 2 with adapter:|^Pairs written \(passing filters\):" "${report}" || true
    echo
  } >> "${SUMMARY_TXT}"

    # Parseo robusto (cutadapt 5.x)
  total=$(awk -F'[[:space:]]+' '/^Total read pairs processed:/ {gsub(",","",$5); print $5}' "${report}")
  r1adp=$(awk -F'[[:space:]]+' '/^ *Read 1 with adapter:/ {gsub(",","",$5); print $5}' "${report}")
  r2adp=$(awk -F'[[:space:]]+' '/^ *Read 2 with adapter:/ {gsub(",","",$5); print $5}' "${report}")
  written=$(awk -F'[[:space:]]+' '/^Pairs written \(passing filters\):/ {gsub(",","",$5); print $5}' "${report}")
  # Extraer porcentaje con sed (compatible BSD)
  pct=$(grep -E "^Pairs written \(passing filters\):" "${report}" | sed -E 's/.*\(([^%]+)%\).*/\1/')
  [[ -z "${pct}" ]] && pct="NA"

  echo -e "${sample}\t${total}\t${r1adp}\t${r2adp}\t${written}\t${pct}" >> "${SUMMARY_TSV}"
  done

echo "[OK] Terminado"
echo "[OK] FASTQ recortados: ${OUT_DIR}"
echo "[OK] Resumen TXT: ${SUMMARY_TXT}"
echo "[OK] Resumen TSV: ${SUMMARY_TSV}"