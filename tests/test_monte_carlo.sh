#!/usr/bin/env bash
# test_monte_carlo.sh — Monte Carlo accuracy validation for startracker
# Usage: bash tests/test_monte_carlo.sh (from repo root)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=== Monte Carlo Validation ==="
cd "$REPO_ROOT"
python3 tools/monte_carlo.py --num-trials 50 --fov 20 --noise 5
echo "=== Monte Carlo Complete ==="
