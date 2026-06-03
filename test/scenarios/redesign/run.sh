#!/usr/bin/env bash
# Run every redesign fixture, writing .dat outputs under out/. Usage: bash run.sh
set -e
cd "$(git rev-parse --show-toplevel)"
rm -rf test/scenarios/redesign/out
mkdir -p test/scenarios/redesign/out
for ini in test/scenarios/redesign/*.ini; do
  echo "=== $ini ==="
  ./class "$ini"
done
