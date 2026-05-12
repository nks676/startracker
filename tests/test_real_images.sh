#!/usr/bin/env bash
# test_real_images.sh — runs the real-image regression against committed
# truth fixtures in tests/data/real_images/. The TIFFs are vendored
# alongside the JSON fixtures (phase 3g.4), so this runs fully offline;
# the Python driver still falls back to a URL download if a new fixture's
# TIFF hasn't been committed yet.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=== Real-image Regression ==="
cd "$REPO_ROOT"
python3 tools/test_real_images.py
echo "=== Real-image Regression Complete ==="
