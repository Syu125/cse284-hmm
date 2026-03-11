#!/usr/bin/env bash
set -euo pipefail

# Single entrypoint for benchmark workflow.
# Usage examples:
#   bash benchmark/run_benchmark.sh all
#   bash benchmark/run_benchmark.sh prepare flare rfmix convert model sweep format

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$ROOT_DIR"

run_prepare() {
  python3 benchmark/scripts/prepare_benchmark_data.py
  python3 benchmark/scripts/fix_genetic_map.py
  python3 benchmark/scripts/make_flare_ref_panel.py
  python3 benchmark/scripts/make_flare_map.py
}

run_flare() {
  java -Xmx8g -jar benchmark/flare.jar \
    ref=benchmark/data/reference_yri_ceu_chr22.vcf.gz \
    ref-panel=benchmark/data/flare_ref_panel_yri_ceu.txt \
    gt=benchmark/data/query_asw_chr22.vcf.gz \
    map=benchmark/data/genetic_map_GRCh37_chr22.flare.map \
    out=benchmark/results/flare_output \
    nthreads=4

  # Optional but recommended for faster downstream reads.
  tabix -p vcf benchmark/results/flare_output.anc.vcf.gz || true
}

run_rfmix() {
  rfmix \
    -f benchmark/data/query_asw_chr22.vcf.gz \
    -r benchmark/data/reference_yri_ceu_chr22.vcf.gz \
    -g benchmark/data/genetic_map_GRCh37_chr22.strict.txt \
    -m data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
    --chromosome=22 \
    -o benchmark/results/rfmix_output
}

run_convert() {
  python3 benchmark/scripts/convert_flare_to_snp_csv.py \
    --flare benchmark/results/flare_output.anc.vcf.gz \
    --tie-policy unknown \
    --out benchmark/predictions/flare/flare_predictions.csv \
    --chromosome 22

  python3 benchmark/scripts/convert_rfmix_to_snp_csv.py \
    --msp benchmark/results/rfmix_output.msp.tsv \
    --vcf benchmark/data/query_asw_chr22.vcf.gz \
    --tie-policy unknown \
    --out benchmark/predictions/rfmix/rfmix_predictions.csv \
    --chromosome 22
}

run_model() {
  python3 benchmark/scripts/export_model_predictions.py \
    --vcf data/processed/chr22_slice.vcf.gz \
    --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
    --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
    --query-pop ASW \
    --out benchmark/predictions/model/model_predictions.csv
}

run_sweep() {
  python3 benchmark/scripts/run_sample_size_sweep.py \
    --sample-sizes 5,20,50 \
    --seeds 0,1,2,3,4 \
    --sample-strategy random \
    --rfmix benchmark/predictions/rfmix/rfmix_predictions.csv \
    --flare benchmark/predictions/flare/flare_predictions.csv \
    --out-prefix sample_sweep
}

run_format() {
  python3 benchmark/scripts/format_sweep_table.py \
    --input benchmark/results/sweep/sample_sweep_summary.csv \
    --out-csv benchmark/plots/sample_sweep_summary_wide.csv \
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
