#!/bin/bash
# Run the full 3-stage nano2pico pipeline (process_nano -> merge_corrections ->
# apply_corrections) on the mChi-200/mLSP-0 higgsino signal point, for every
# year, in parallel. Mirrors the single-file 2016 commands used earlier in
# this validation, generalized across years.
#
# "2025" reuses the 2024 physical MC file but is processed with 2025
# corrections via --is25mc (see IS25MC map below). This requires
# process_nano.exe to have the --is25mc flag
# added (handoff diffs given separately) -- it's a general mechanism, not
# specific to this script, so the same pattern can be reused for background
# samples / other signal points once those get their own batch script.
#
# Run this from the nano2pico repo root (where run/process_nano.exe etc. live).
#
# Usage:
#   bash scripts/run_all_years_mChi200.sh            # process all years
#   DRY_RUN=1 bash scripts/run_all_years_mChi200.sh   # just show what file
#                                                      # would be used per
#                                                      # year, don't run
#   MAX_PARALLEL=3 bash scripts/run_all_years_mChi200.sh   # cap concurrency
#
# ALWAYS do a DRY_RUN=1 pass first to confirm the file-discovery glob finds
# exactly one file per year before launching 9 parallel processing jobs.

set -o pipefail

MASS_CHI=200
MASS_LSP=0
MAX_PARALLEL="${MAX_PARALLEL:-3}"
DRY_RUN="${DRY_RUN:-0}"
BASE_OUT="test_sf_higgsino_allyears"

declare -A YEAR_DIRS=(
  ["2016"]="/net/cms11/cms11r0/pico/NanoAODv9/nano/2016/SMS-TChiHH-Hto2G-FullSIM-Central-Run2_split"
  ["2016APV"]="/net/cms11/cms11r0/pico/NanoAODv9/nano/2016APV/SMS-TChiHH-Hto2G-FullSIM-Central-Run2_split"
  ["2017"]="/net/cms11/cms11r0/pico/NanoAODv9/nano/2017/SMS-TChiHH-Hto2G-FullSIM-Central-Run2_split"
  ["2018"]="/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/SMS-TChiHH-Hto2G-FullSIM-Central-Run2_split"
  ["2022"]="/net/cms11/cms11r0/pico/NanoAODv12/nano/2022/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
  ["2022EE"]="/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
  ["2023"]="/net/cms11/cms11r0/pico/NanoAODv12/nano/2023/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
  ["2023BPix"]="/net/cms11/cms11r0/pico/NanoAODv12/nano/2023BPix/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
  ["2024"]="/net/cms11/cms11r0/pico/NanoAODv15/nano/2024/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
  ["2025"]="/net/cms11/cms11r0/pico/NanoAODv15/nano/2024/SMS-TChiHH-Hto2G-FullSIM-Central-Run3_split"
)


# Years listed here get --is25mc passed to process_nano.exe, which flips the
# filename-detected year 2024 -> 2025 and keeps only odd-numbered events. The
# 2024 job keeps the even ones, so the two eras use statistically independent
# halves of the same physical MC file.
declare -A IS25MC=( ["2025"]="1" )

process_year() {
  local year="$1"
  local in_dir="${YEAR_DIRS[$year]}"
  local out_dir="${BASE_OUT}/${year}"
  local log_prefix="[$year]"

  if [ ! -d "$in_dir" ]; then
    echo "${log_prefix} ERROR: input directory does not exist: $in_dir" >&2
    return 1
  fi

  # exact match on mChi-200_mLSP-0_ (trailing underscore avoids accidentally
  # matching mChi-2000_mLSP-0_ etc.)
  local in_file
  in_file=$(ls "$in_dir" 2>/dev/null | grep -E "mChi-${MASS_CHI}_mLSP-${MASS_LSP}_" | head -1)
  local n_matches
  n_matches=$(ls "$in_dir" 2>/dev/null | grep -cE "mChi-${MASS_CHI}_mLSP-${MASS_LSP}_")

  if [ -z "$in_file" ]; then
    echo "${log_prefix} ERROR: no mChi-${MASS_CHI}_mLSP-${MASS_LSP} file found in $in_dir" >&2
    return 1
  fi
  if [ "$n_matches" -gt 1 ]; then
    echo "${log_prefix} WARNING: ${n_matches} files matched mChi-${MASS_CHI}_mLSP-${MASS_LSP}, using first: $in_file" >&2
  fi

  echo "${log_prefix} input file: $in_file"

  local is25_args=()
  if [ -n "${IS25MC[$year]:-}" ]; then
    echo "${log_prefix} 2024 MC processed with 2025 corrections (odd events only)"
    is25_args=(--is25mc)
  fi

  if [ "$DRY_RUN" == "1" ]; then
    return 0
  fi

  mkdir -p "${out_dir}/raw_pico" "${out_dir}/wgt_sums" "${out_dir}/corrections" "${out_dir}/unskimmed"

  echo "${log_prefix} Step 1: process_nano.exe"
  if ! run/process_nano.exe --in_dir "$in_dir" --in_file "$in_file" --out_dir "$out_dir" \
      ${is25_args[@]+"${is25_args[@]}"} \
      > "${out_dir}/step1.log" 2>&1; then
    echo "${log_prefix} FAILED at step 1 -- see ${out_dir}/step1.log" >&2
    return 1
  fi

  local corr_name
  corr_name="corr_$(echo "$in_file" | sed 's/\.root$//').root"
  echo "${log_prefix} Step 2: merge_corrections.exe"
  # merge_corrections.exe takes no flags -- its long_options array is empty
  if ! run/merge_corrections.exe "${out_dir}/corrections/${corr_name}" \
      "${out_dir}/wgt_sums/wgt_sums_${in_file}" \
      > "${out_dir}/step2.log" 2>&1; then
    echo "${log_prefix} FAILED at step 2 -- see ${out_dir}/step2.log" >&2
    return 1
  fi

  echo "${log_prefix} Step 3: apply_corrections.exe"
  if ! run/apply_corrections.exe --in_dir "${out_dir}/raw_pico/" \
      --in_file "raw_pico_${in_file}" --corr_file "$corr_name" \
      > "${out_dir}/step3.log" 2>&1; then
    echo "${log_prefix} FAILED at step 3 -- see ${out_dir}/step3.log" >&2
    return 1
  fi

  echo "${log_prefix} DONE -> ${out_dir}/unskimmed/pico_${in_file}"
}

years=("${!YEAR_DIRS[@]}")
n=${#years[@]}
echo "Processing ${n} years, up to ${MAX_PARALLEL} in parallel at a time. DRY_RUN=${DRY_RUN}"
echo ""

for ((i=0; i<n; i+=MAX_PARALLEL)); do
  batch=("${years[@]:i:MAX_PARALLEL}")
  for year in "${batch[@]}"; do
    process_year "$year" &
  done
  wait
done

echo ""
echo "All years finished (check output above for any ERROR/FAILED lines)."
echo "Outputs under: ${BASE_OUT}/<year>/unskimmed/"
