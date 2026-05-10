#!/usr/bin/env bash
# test_full_pipeline.sh — End-to-end full pipeline test for startracker
# Usage: bash test_full_pipeline.sh (from any directory)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=== Full Pipeline Test ==="

# 1. Generate synthetic star image (identity quaternion → boresight = +Z)
echo "[1/4] Generating synthetic star image..."
cd "$REPO_ROOT/tools"
python3 generate_synthetic_data.py --quat 0 0 0 1 --out "$REPO_ROOT/data/test_phase1"

# 2. Generate binary star catalog
echo "[2/4] Generating star catalog..."
python3 generate_catalog.py

# 3. Run the C++ startracker engine
echo "[3/4] Running startracker engine..."
cd "$REPO_ROOT"
./build/startracker \
  data/test_phase1/synthetic_starfield.png \
  data/catalog_stars.bin \
  data/catalog_pairs.bin \
  20 | tee /tmp/startracker_result.txt

# 4. Verify accuracy
echo "[4/4] Verifying accuracy..."
cd "$REPO_ROOT/tools"
python3 verify_accuracy.py \
  --truth "$REPO_ROOT/data/test_phase1/truth.json" \
  --result-file /tmp/startracker_result.txt

echo ""
echo "=== Full Pipeline Test Complete ==="
