#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$analysis_dir/scripts/01_extract_parikh_pgtd.R"
Rscript "$analysis_dir/scripts/02_make_expression_panels.R"

echo "Expression panels B and C completed"
