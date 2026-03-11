#!/usr/bin/env bash
set -euo pipefail

# Single entrypoint for benchmark workflow.
# Usage examples:
#   bash benchmark/run_benchmark.sh all
#   bash benchmark/run_benchmark.sh prepare flare rfmix convert model sweep format

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$ROOT_DIR"

PERF_DIR="benchmark/results/performance"
PERF_CSV="${PERF_DIR}/stage_performance.csv"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"

mkdir -p "$PERF_DIR"
if [[ ! -f "$PERF_CSV" ]]; then
  echo "run_id,timestamp_utc,stage,task,elapsed_seconds,max_rss_kb" > "$PERF_CSV"
fi

run_profiled() {
  local stage="$1"
  local task="$2"
  shift 2

  local tmp_file
  tmp_file="$(mktemp)"
  local start_ns end_ns elapsed_seconds exit_code max_rss_kb

  start_ns="$(date +%s%N)"
  if [[ -x /usr/bin/time ]]; then
    set +e
    /usr/bin/time -f "max_rss_kb=%M" -o "$tmp_file" "$@"
    exit_code=$?
    set -e
  else
    set +e
    "$@"
    exit_code=$?
    set -e
    echo "max_rss_kb=NA" > "$tmp_file"
  fi
  end_ns="$(date +%s%N)"

  elapsed_seconds="$(awk "BEGIN { printf \"%.3f\", (${end_ns}-${start_ns})/1000000000 }")"
  max_rss_kb="$(awk -F= '/^max_rss_kb=/{print $2}' "$tmp_file" | tail -n1)"
  if [[ -z "$max_rss_kb" ]]; then
    max_rss_kb="NA"
  fi

  printf "%s,%s,%s,%s,%s,%s\n" \
    "$RUN_ID" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$stage" \
    "$task" \
    "$elapsed_seconds" \
    "$max_rss_kb" >> "$PERF_CSV"

  rm -f "$tmp_file"

  if [[ $exit_code -ne 0 ]]; then
    echo "[ERROR] ${stage}/${task} failed with exit code ${exit_code}" >&2
    return $exit_code
  fi
}

run_prepare() {
  run_profiled prepare prepare_benchmark_data python3 benchmark/scripts/prepare_benchmark_data.py
  run_profiled prepare fix_genetic_map python3 benchmark/scripts/fix_genetic_map.py
  run_profiled prepare make_flare_ref_panel python3 benchmark/scripts/make_flare_ref_panel.py
  run_profiled prepare make_flare_map python3 benchmark/scripts/make_flare_map.py
}

run_flare() {
  run_profiled flare flare_inference java -Xmx8g -jar benchmark/flare.jar \
    ref=benchmark/data/reference_yri_ceu_chr22.vcf.gz \
    ref-panel=benchmark/data/flare_ref_panel_yri_ceu.txt \
    gt=benchmark/data/query_asw_chr22.vcf.gz \
    map=benchmark/data/genetic_map_GRCh37_chr22.flare.map \
    out=benchmark/results/flare_output \
    nthreads=4

  # Optional but recommended for faster downstream reads.
  run_profiled flare flare_tabix tabix -p vcf benchmark/results/flare_output.anc.vcf.gz || true
}

run_rfmix() {
  run_profiled rfmix rfmix_inference rfmix \
    -f benchmark/data/query_asw_chr22.vcf.gz \
    -r benchmark/data/reference_yri_ceu_chr22.vcf.gz \
    -g benchmark/data/genetic_map_GRCh37_chr22.strict.txt \
    -m data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
    --chromosome=22 \
    -o benchmark/results/rfmix_output
}

run_convert() {
  run_profiled convert convert_flare_to_snp_csv python3 benchmark/scripts/convert_flare_to_snp_csv.py \
    --flare benchmark/results/flare_output.anc.vcf.gz \
    --tie-policy unknown \
    --out benchmark/predictions/flare/flare_predictions.csv \
    --chromosome 22

  run_profiled convert convert_rfmix_to_snp_csv python3 benchmark/scripts/convert_rfmix_to_snp_csv.py \
    --msp benchmark/results/rfmix_output.msp.tsv \
    --vcf benchmark/data/query_asw_chr22.vcf.gz \
    --tie-policy unknown \
    --out benchmark/predictions/rfmix/rfmix_predictions.csv \
    --chromosome 22
}

run_model() {
  run_profiled model model_export python3 benchmark/scripts/export_model_predictions.py \
    --vcf data/processed/chr22_slice.vcf.gz \
    --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
    --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
    --query-pop ASW \
    --out benchmark/predictions/model/model_predictions.csv
}

run_sweep() {
  run_profiled sweep run_sample_size_sweep python3 benchmark/scripts/run_sample_size_sweep.py \
    --sample-sizes 5,20,50 \
    --seeds 0,1,2,3,4 \
    --sample-strategy random \
    --rfmix benchmark/predictions/rfmix/rfmix_predictions.csv \
    --flare benchmark/predictions/flare/flare_predictions.csv \
    --out-prefix sample_sweep
}

run_format() {
  run_profiled format format_sweep_table python3 benchmark/scripts/format_sweep_table.py \
    --input benchmark/results/sweep/sample_sweep_summary.csv \
    --out-csv benchmark/plots/sample_sweep_summary_wide.csv \
    --out-dir benchmark/plots

  run_profiled format format_performance_table python3 benchmark/scripts/format_performance_table.py \
    --input benchmark/results/performance/stage_performance.csv \
    --out-dir benchmark/plots
}

run_all() {
  run_prepare
  run_flare
  run_rfmix
  run_convert
  run_model
  run_sweep
  run_format
}

if [[ $# -eq 0 ]]; then
  echo "No step specified. Running: all"
  run_all
  exit 0
fi

for step in "$@"; do
  case "$step" in
    all) run_all ;;
    prepare) run_prepare ;;
    flare) run_flare ;;
    rfmix) run_rfmix ;;
    convert) run_convert ;;
    model) run_model ;;
    sweep) run_sweep ;;
    format) run_format ;;
    *) echo "Unknown step: $step"; exit 1 ;;
  esac
done

echo "Saved performance profile rows to: ${PERF_CSV}"
