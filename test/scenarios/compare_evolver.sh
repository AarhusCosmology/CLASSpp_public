#!/usr/bin/env bash
# Compare rkdp45 (evolver=2) against the ndf15 default across all test scenarios.
# Pass = every scenario's outputs agree within test/scenarios/compare_tol.py (RTOL=1e-3).
# Usage: test/scenarios/compare_evolver.sh [path-to-class-binary]
set -u
CLASS="${1:-./class}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(mktemp -d)"
overall=0

# Scenarios that are intentionally invalid (expected to error under ndf15 too,
# e.g. error-path rejection tests). They are not solver comparisons, so skip
# them rather than letting them mark the whole run as failed.
SKIP=" dncdm_dr_bare_omega "

for ini in "$HERE"/*.ini; do
  name="$(basename "$ini" .ini)"
  case "$SKIP" in
    *" $name "*) echo "SKIP $name (intentionally-invalid scenario)"; continue ;;
  esac
  refd="$WORK/$name/ref"; newd="$WORK/$name/new"
  mkdir -p "$refd" "$newd"
  # strip any existing root=, then append our own + evolver setting
  grep -vE '^\s*root\s*=' "$ini" > "$WORK/$name/base.ini"
  { cat "$WORK/$name/base.ini"; printf "\nroot = %s/out_\n" "$refd"; }                       > "$WORK/$name/ref.ini"
  { cat "$WORK/$name/base.ini"; printf "\nroot = %s/out_\nevolver = 2\n" "$newd"; }          > "$WORK/$name/new.ini"
  if ! "$CLASS" "$WORK/$name/ref.ini" > "$WORK/$name/ref.log" 2>&1; then
    echo "RUNFAIL(ndf15) $name"; overall=1; continue
  fi
  if ! "$CLASS" "$WORK/$name/new.ini" > "$WORK/$name/new.log" 2>&1; then
    echo "RUNFAIL(rkdp45) $name (see $WORK/$name/new.log)"; overall=1; continue
  fi
  if ! ls "$refd"/out_*.dat >/dev/null 2>&1; then
    echo "NOOUT $name (scenario produces no .dat)"; continue
  fi
  echo "== $name =="
  if ! python3 "$HERE/compare_tol.py" "$refd" "$newd" 'out_*.dat'; then
    overall=1
  fi
done
echo "WORKDIR=$WORK"
[ "$overall" -eq 0 ] && echo "ALL SCENARIOS PASS" || echo "SOME SCENARIOS FAILED"
exit $overall
