#!/usr/bin/env bash
# Distill an output tree into a served bundle and upload to S3.
#
# Run from your local machine (or any host with `aws` configured and
# the repo checked out).  The input `output_dir` can be a local
# pipeline run or an output tree rsynced down from Bouchet.
#
# Usage:
#     deploy/sync_to_s3.sh <output_dir> <version>
#     deploy/sync_to_s3.sh output v1.0.0
#     BUCKET=corpus-bundles-staging deploy/sync_to_s3.sh output v1.0.0-rc1
#
# Env overrides:
#     BUCKET      — target S3 bucket (default: corpus-bundles)
#     BUNDLE_DIR  — where to stage the distilled bundle on local disk
#                   (default: /tmp/corpus-bundle-<version>).  Kept
#                   after the upload so reruns can reuse mtime-based
#                   skip in package_for_serve.py.
#     AWS_PROFILE — honoured if set; otherwise falls through to the
#                   default aws cli configuration.
#     INCLUDE_PDFS=1 — also copy processed.pdf per document (~3 GB
#                   extra at 2000 papers).  Default: excluded.
#
# Exits non-zero on any failure, so this is safe to chain in a CI
# pipeline.

set -euo pipefail

die() { printf 'sync_to_s3: %s\n' "$*" >&2; exit 1; }

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

OUTPUT_DIR="${1:-}"
VERSION="${2:-}"
[[ -n "$OUTPUT_DIR" && -n "$VERSION" ]] || die "usage: $0 <output_dir> <version>"
[[ -d "$OUTPUT_DIR" ]] || die "output_dir not a directory: $OUTPUT_DIR"
[[ -d "$OUTPUT_DIR/documents" ]] || die "no documents/ under $OUTPUT_DIR — is this really a pipeline output?"

BUCKET="${BUCKET:-corpus-bundles}"
BUNDLE_DIR="${BUNDLE_DIR:-/tmp/corpus-bundle-${VERSION}}"

# Using a plain string instead of an array because macOS ships bash 3.2,
# which errors on empty-array expansion under `set -u` ("unbound
# variable").  Safer to build a flag string and let the shell split on
# whitespace; we only pass fixed flags here, never user input.
EXTRA_FLAGS=""
if [[ "${INCLUDE_PDFS:-0}" == "1" ]]; then
    EXTRA_FLAGS="--include-pdfs"
fi

command -v aws >/dev/null 2>&1 || die "aws cli not on PATH; install awscli"

echo "── Distill ────────────────────────────────────────────────────"
echo "  input:  $OUTPUT_DIR"
echo "  bundle: $BUNDLE_DIR"
echo "  version: $VERSION"
# shellcheck disable=SC2086
python "$REPO_DIR/package_for_serve.py" \
    "$OUTPUT_DIR" "$BUNDLE_DIR" \
    --version "$VERSION" \
    $EXTRA_FLAGS

echo
echo "── Upload ─────────────────────────────────────────────────────"
S3_URL="s3://${BUCKET}/${VERSION}/"
echo "  target: $S3_URL"
# --delete prunes objects in S3 that aren't in the local bundle.
# Safe because each version is a distinct S3 prefix; older versions
# live under s3://$BUCKET/v0.9.0/ and are untouched.
aws s3 sync --delete "$BUNDLE_DIR/" "$S3_URL"

echo
echo "── Done ───────────────────────────────────────────────────────"
echo "  Deploy on EC2:"
echo "      sudo -u corpus /srv/corpus/repo/deploy/update.sh $VERSION"
