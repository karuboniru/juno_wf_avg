#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <eos_pattern>"
    echo "  Example: $0 '/eos/juno/juno-rtraw/J25.7.1/global_trigger/00015000/00015000/15010/*.rtraw'"
    exit 1
fi

EOS_PATTERN="$1"
ROOT_PREFIX="root://junoeos01.ihep.ac.cn/"

EOS_DIR="$(dirname "$EOS_PATTERN")"
RUN_ID="$(basename "$EOS_DIR")"
OUTFILE="${RUN_ID}.list"

echo "Pattern : $EOS_PATTERN"
echo "Run ID  : $RUN_ID"

eos ls "$EOS_PATTERN" | while IFS= read -r filename; do
    [[ -z "$filename" ]] && continue
    echo "${ROOT_PREFIX}${EOS_DIR}/${filename}"
done > "$OUTFILE"

count=$(wc -l < "$OUTFILE")
echo "Wrote $count files to $OUTFILE"
