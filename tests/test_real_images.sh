#!/usr/bin/env bash
# test_real_images.sh — runs the real-image regression against committed
# truth fixtures in tests/data/real_images/. The Python driver downloads
# the source TIFFs on first run (cached in data/real_images/).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=== Real-image Regression ==="
cd "$REPO_ROOT"
python3 tools/test_real_images.py
echo "=== Real-image Regression Complete ==="
