#!/bin/bash
#
# Efficiency test: Normal MC vs Exchange MC
# Single run per method, stderr estimated from internal Sample bins.
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MC_BIN="$PROJECT_ROOT/src/build/MC_simple"
SAMPLE_DIR="$PROJECT_ROOT/samples/square_L16_Ising"
BASE_SEED=10000

if [[ ! -x "$MC_BIN" ]]; then
    echo "Error: MC_simple not found. Build first:"
    echo "  cd src && cmake -S . -B build && cmake --build build -j"
    exit 1
fi

cd "$SCRIPT_DIR"
rm -f normal_sample_bins.dat exchange_sample_bins.dat \
    normal_single_run.log exchange_single_run.log \
    _tmp_param_normal.def _tmp_param_exchange.def

BASE_PARAM="$SAMPLE_DIR/param.def"
LATTICE="$SAMPLE_DIR/lattice.def"
INTERACTION="$SAMPLE_DIR/interaction.def"

if [[ ! -f "$BASE_PARAM" || ! -f "$LATTICE" || ! -f "$INTERACTION" ]]; then
    echo "Error: param.def / lattice.def / interaction.def not found in $SAMPLE_DIR"
    exit 1
fi

cp "$BASE_PARAM" _tmp_param_normal.def
cp "$BASE_PARAM" _tmp_param_exchange.def

if grep -q "^[[:space:]]*exchange_interval[[:space:]]*=" _tmp_param_normal.def; then
    sed -E -i.bak 's/^[[:space:]]*exchange_interval[[:space:]]*=.*/exchange_interval = 0/' _tmp_param_normal.def
else
    echo "exchange_interval = 0" >> _tmp_param_normal.def
fi

if grep -q "^[[:space:]]*exchange_interval[[:space:]]*=" _tmp_param_exchange.def; then
    sed -E -i.bak 's/^[[:space:]]*exchange_interval[[:space:]]*=.*/exchange_interval = 1/' _tmp_param_exchange.def
else
    echo "exchange_interval = 1" >> _tmp_param_exchange.def
fi

rm -f _tmp_param_normal.def.bak _tmp_param_exchange.def.bak

echo "Running efficiency test (stderr from Sample bins in param.def)"
echo ""

echo "Running Normal MC (exchange_interval=0)..."
MC_SIMPLE_SEED=$BASE_SEED "$MC_BIN" _tmp_param_normal.def "$LATTICE" "$INTERACTION" \
    > normal_single_run.log 2>&1
cp MC_simple_sample_bins.dat normal_sample_bins.dat

echo "Running Exchange MC (exchange_interval=1)..."
MC_SIMPLE_SEED=$((BASE_SEED + 1)) "$MC_BIN" _tmp_param_exchange.def "$LATTICE" "$INTERACTION" \
    > exchange_single_run.log 2>&1
cp MC_simple_sample_bins.dat exchange_sample_bins.dat

rm -f _tmp_param_normal.def _tmp_param_exchange.def

echo ""
echo "Running analysis..."
python3 analyze_sample_bins.py

echo ""
echo "Done. See efficiency_comparison.png"
