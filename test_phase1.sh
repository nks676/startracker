#!/usr/bin/env bash
# test_phase1.sh — End-to-end Phase 1 smoke test for startracker
# Runs inside Docker: docker run --rm startracker bash /app/test_phase1.sh
set -euo pipefail

echo "=== Phase 1 Smoke Test ==="

# 1. Generate synthetic star image (identity quaternion → boresight = +Z)
echo "[1/4] Generating synthetic star image..."
cd /app/tools
python3 generate_synthetic_data.py --quat 0 0 0 1 --out /app/data/test_phase1

# 2. Generate binary star catalog
echo "[2/4] Generating star catalog..."
python3 generate_catalog.py

# 3. Run the C++ startracker engine
echo "[3/4] Running startracker engine..."
cd /app
./build/startracker \
  data/test_phase1/synthetic_starfield.png \
  data/catalog_stars.bin \
  data/catalog_pairs.bin \
  20 | tee /tmp/result.txt

# 4. Verify accuracy (use --result-file to avoid argparse issues with
#    negative scientific-notation values)
echo "[4/4] Verifying accuracy..."
cd /app/tools
python3 verify_accuracy.py \
  --truth /app/data/test_phase1/truth.json \
  --result-file /tmp/result.txt

echo ""
echo "=== Phase 1 Smoke Test Complete ==="
