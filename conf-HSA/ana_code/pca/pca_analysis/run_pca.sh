#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

python3 combined_pca.py \
    --input-dir ../com_md/ \
    --output-dir . \
    --top-k-eigs 0 \
    --align-chunk 5000 \
    --cov-chunk 5000 \
    --force
