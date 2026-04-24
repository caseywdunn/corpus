#!/usr/bin/env bash
# Fetch a bundle version from S3 onto the EC2 host, atomically swap
# the active bundle symlink, and restart the corpus-mcp service.
#
# Intended to be run as the `corpus` user on the EC2 host:
#
#     sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.0.0
#
# Idempotent: re-running with the same version is a no-op apart
# from the restart (rsync --delete inside `aws s3 sync` will prune
# any drift).  Keeps prior version directories in place for
# instant rollback — see DEPLOY.md.
#
# Env overrides:
#     BUCKET — source S3 bucket (default: corpus-bundles)
#     ROOT   — root of the corpus install on this host
#              (default: /srv/corpus)
#     SERVICE — systemd unit name (default: corpus-mcp)

set -euo pipefail

die() { printf 'update: %s\n' "$*" >&2; exit 1; }

VERSION="${1:-}"
[[ -n "$VERSION" ]] || die "usage: $0 <version>   (e.g. v1.0.0)"

BUCKET="${BUCKET:-corpus-bundles}"
ROOT="${ROOT:-/srv/corpus}"
SERVICE="${SERVICE:-corpus-mcp}"

BUNDLES_DIR="$ROOT/bundles"
SYMLINK="$ROOT/bundle"
TARGET_DIR="$BUNDLES_DIR/$VERSION"

[[ -d "$ROOT" ]] || die "$ROOT does not exist — has the host been bootstrapped?"
command -v aws >/dev/null 2>&1 || die "aws cli not on PATH"

mkdir -p "$TARGET_DIR"

echo "── Fetch ──────────────────────────────────────────────────────"
S3_URL="s3://${BUCKET}/${VERSION}/"
echo "  source: $S3_URL"
echo "  target: $TARGET_DIR"
aws s3 sync --delete "$S3_URL" "$TARGET_DIR"

# Sanity check the manifest landed.  If not, the S3 prefix is wrong
# or the upload was incomplete — bail without touching the symlink.
[[ -f "$TARGET_DIR/bundle_manifest.json" ]] || \
    die "no bundle_manifest.json in $TARGET_DIR — aborting before symlink swap"

echo
echo "── Activate ───────────────────────────────────────────────────"
# Atomic symlink swap: create a sibling, rename over.  Nothing sees
# a half-updated state.
TMP_LINK="${SYMLINK}.new"
ln -sfn "$TARGET_DIR" "$TMP_LINK"
mv -T "$TMP_LINK" "$SYMLINK"
echo "  $SYMLINK -> $(readlink "$SYMLINK")"

echo
echo "── Restart ────────────────────────────────────────────────────"
# `restart` rather than `reload` because mcp_server.py builds its
# in-memory index at startup; a SIGHUP handler isn't plumbed in.
# Brief interruption (a few seconds); in-flight requests drain.
if command -v systemctl >/dev/null 2>&1; then
    sudo systemctl restart "$SERVICE"
    sudo systemctl is-active "$SERVICE" --quiet || \
        die "$SERVICE did not come back up — journalctl -u $SERVICE for logs"
    echo "  $SERVICE is active"
else
    echo "  systemctl not available; restart the service manually" >&2
fi

echo
echo "── Done ───────────────────────────────────────────────────────"
echo "  Active bundle: $VERSION"
echo "  Rollback:      $0 <older-version>"
