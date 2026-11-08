#!/usr/bin/env bash
# Optional one-time bootstrap for qualification:* labels. Never deletes or
# overwrites existing labels. Requires an authenticated GitHub CLI.
set -euo pipefail

if ! command -v gh >/dev/null 2>&1; then
  echo "GitHub CLI (gh) is required: https://cli.github.com/" >&2
  exit 2
fi

gh auth status >/dev/null
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
existing="$(gh label list --limit 500 --json name --jq '.[].name')"

while IFS=$'\t' read -r label description; do
  if grep -Fqx -- "$label" <<< "$existing"; then
    echo "Already exists: $label"
  else
    gh label create "$label" --color 1D76DB --description "$description"
    echo "Created: $label"
  fi
done < <(python3 "$ROOT/scripts/ci/qualification.py" labels --tsv)