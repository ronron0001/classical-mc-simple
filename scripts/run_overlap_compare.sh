#!/usr/bin/env bash
# Run MC_simple with EXMC and without EXMC, save overlap.dat for each.
# Usage: from repo root, run:
#   ./scripts/run_overlap_compare.sh
# or from samples/square_L16_Ising:
#   ../../scripts/run_overlap_compare.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MC_SIMPLE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC="$MC_SIMPLE_DIR/src"
SAMPLE="$MC_SIMPLE_DIR/samples/square_L16_Ising"
BIN="$SRC/build/MC_simple"

cd "$SAMPLE"

echo "=== Run with EXMC (use_exmc=1) ==="
"$BIN" param.def lattice.def interaction.def
mv -f overlap.dat overlap_exmc.dat 2>/dev/null || true

echo ""
echo "=== Run without EXMC (use_exmc=0) ==="
"$BIN" param_no_exmc.def lattice.def interaction.def
mv -f overlap.dat overlap_no_exmc.dat 2>/dev/null || true

echo ""
python3 "$SCRIPT_DIR/plot_overlap.py" "$SAMPLE"
cd "$SAMPLE" && gnuplot "$SCRIPT_DIR/plot_overlap.gp" 2>/dev/null && echo "Saved $SAMPLE/overlap_compare.png" || true
echo "Done. See overlap_compare.tsv and overlap_compare.png in $SAMPLE"
