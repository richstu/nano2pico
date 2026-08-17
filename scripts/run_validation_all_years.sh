#!/bin/bash
# Run the Type 1 (correctness) and Type 2 (distributional impact) photon-SF
# validation scripts against each year's processed pico file, i.e. the
# output of scripts/run_all_years_mChi200.sh.
#
# 2025 is deliberately excluded -- it's not processed yet (deferred, see
# run_all_years_mChi200.sh's FORCE_YEAR mechanism).
#
# Run this from the nano2pico repo root (where scripts/validate_photon_sf.py
# etc. live).
#
# Usage:
#   bash scripts/run_validation_all_years.sh
#   BASE_OUT=test_sf_higgsino_allyears bash scripts/run_validation_all_years.sh
#   WEIGHT_BRANCHES=w_lumi,w_photon,w_phshape bash scripts/run_validation_all_years.sh

set -u

MASS_CHI=200
MASS_LSP=0
BASE_OUT="${BASE_OUT:-test_sf_higgsino_allyears}"
VALID_OUT="${VALID_OUT:-${BASE_OUT}/validation_reports}"
# Restrict Type 2's "weight" to just the validated photon SF chain --
# other SF components (b-tag etc.) aren't validated yet.
WEIGHT_BRANCHES="${WEIGHT_BRANCHES:-w_lumi,w_photon,w_phshape}"

YEARS=(2016 2016APV 2017 2018 2022 2022EE 2023 2023BPix 2024)

mkdir -p "$VALID_OUT"

n_ok=0
n_fail=0

for year in "${YEARS[@]}"; do
  unskimmed_dir="${BASE_OUT}/${year}/unskimmed"
  log_prefix="[$year]"

  if [ ! -d "$unskimmed_dir" ]; then
    echo "${log_prefix} ERROR: $unskimmed_dir does not exist -- skipping" >&2
    n_fail=$((n_fail+1))
    continue
  fi

  # apply_corrections.cxx writes unskimmed/pico_<original_nanoaod_filename>
  pico_file=$(ls "$unskimmed_dir" 2>/dev/null | grep -E "^pico_.*mChi-${MASS_CHI}_mLSP-${MASS_LSP}_" | head -1)
  if [ -z "$pico_file" ]; then
    echo "${log_prefix} ERROR: no pico_*mChi-${MASS_CHI}_mLSP-${MASS_LSP}* file found in $unskimmed_dir -- skipping" >&2
    n_fail=$((n_fail+1))
    continue
  fi

  in_path="${unskimmed_dir}/${pico_file}"
  echo "${log_prefix} validating $in_path"

  ok=1

  echo "${log_prefix}   Type 1 (photon SF correctness)..."
  if ! python3 scripts/validate_photon_sf.py "$in_path" \
      --out "${VALID_OUT}/${year}_photon_sf_type1.pdf" \
      > "${VALID_OUT}/${year}_type1.log" 2>&1; then
    echo "${log_prefix}   Type 1 FAILED -- see ${VALID_OUT}/${year}_type1.log" >&2
    ok=0
  fi

  echo "${log_prefix}   Type 2 (distributional impact, inclusive)..."
  if ! python3 scripts/validate_photon_sf_type2.py "$in_path" \
      --out "${VALID_OUT}/${year}_photon_sf_type2_inclusive.pdf" \
      --weight-branches "$WEIGHT_BRANCHES" \
      > "${VALID_OUT}/${year}_type2_inclusive.log" 2>&1; then
    echo "${log_prefix}   Type 2 (inclusive) FAILED -- see ${VALID_OUT}/${year}_type2_inclusive.log" >&2
    ok=0
  fi

  echo "${log_prefix}   Type 2 (distributional impact, photon SR cuts)..."
  if ! python3 scripts/validate_photon_sf_type2.py "$in_path" \
      --out "${VALID_OUT}/${year}_photon_sf_type2_srcuts.pdf" \
      --weight-branches "$WEIGHT_BRANCHES" \
      --photon-sr-cuts \
      > "${VALID_OUT}/${year}_type2_srcuts.log" 2>&1; then
    echo "${log_prefix}   Type 2 (SR cuts) FAILED -- see ${VALID_OUT}/${year}_type2_srcuts.log" >&2
    ok=0
  fi

  if [ "$ok" == "1" ]; then
    echo "${log_prefix} DONE -- 3 PDFs in ${VALID_OUT}/"
    n_ok=$((n_ok+1))
  else
    n_fail=$((n_fail+1))
  fi
done

echo ""
echo "Finished: ${n_ok} years OK, ${n_fail} years with errors."
echo "Reports under: ${VALID_OUT}/"

# Lightweight sanity check: confirm each PDF actually has pages (catches
# silent empty-report failures that still exit 0).
python3 - "$VALID_OUT" <<'PYEOF'
import sys, glob, os
out_dir = sys.argv[1]
try:
    from pypdf import PdfReader
except ImportError:
    print("(pypdf not installed -- skipping PDF page-count sanity check)")
    sys.exit(0)

pdfs = sorted(glob.glob(os.path.join(out_dir, "*.pdf")))
if not pdfs:
    print("No PDFs found to check.")
    sys.exit(0)

print("\nPDF page-count sanity check:")
for p in pdfs:
    try:
        n = len(PdfReader(p).pages)
        flag = "" if n > 0 else "  <-- EMPTY, investigate"
        print(f"  {os.path.basename(p)}: {n} pages{flag}")
    except Exception as exc:
        print(f"  {os.path.basename(p)}: ERROR reading PDF ({exc})")
PYEOF
