#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
wolframscript -file run_tests.wls
