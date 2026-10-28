#!/bin/bash
# Build the ACES manual: HTML into docs/_build/html (pass "latexpdf" for a PDF).
set -euo pipefail
cd "$(dirname "$0")"
target=${1:-html}
if [ "$target" = html ]; then
  python3 -m sphinx -b html -W --keep-going . _build/html
else
  python3 -m sphinx -M "$target" . _build
fi
