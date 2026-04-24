# AWS deployment runbook

This walks through the first deploy and the update flow.  The
architecture is spelled out in [PLAN.md §10](../PLAN.md); this
document is the operational instructions.

## Layout on EC2

```
/srv/corpus/
  repo/               git checkout of this repository (code)
  venv/               Python virtualenv with the server's deps
  bundles/
    v1.0.0/           copy of s3://corpus-bundles/v1.0.0/
    v1.1.0/           copy of s3://corpus-bundles/v1.1.0/
    …                 older versions retained for instant rollback
  bundle -> bundles/v1.1.0/   atomic symlink to the active version

/etc/corpus/
  mcp.token           bearer token, mode 600, owned root:corpus

/etc/systemd/system/corpus-mcp.service    the unit from deploy/
```

## One-time setup

### 1. AWS

- **S3 bucket** — e.g. `corpus-bundles`.  Enable versioning.  IAM
  policy giving the EC2 instance role `s3:GetObject` and
  `s3:ListBucket`.  Optional lifecycle rule moving objects in
  unused version prefixes to Glacier after 90 days.
- **EC2 instance** — `t3.small` is plenty for ~20 collaborators.
  Ubuntu 24.04 LTS.  10 GB gp3 EBS.  Security group allows 22 from
  your IP and 8080 from the ALB only.
- **ACM cert** — for your chosen domain.  Attach to the ALB (or
  CloudFront; CloudFront → ALB → EC2 is the pattern from PLAN §10).
- **Route 53 record** — CNAME or alias pointing at CloudFront.
- **IAM role** attached to the EC2 instance giving it read on the S3
  bucket (no keys in files).

### 2. EC2 bootstrap

```bash
# On the EC2 host as root or via sudo:
sudo useradd --system --home-dir /srv/corpus --create-home corpus
sudo mkdir -p /srv/corpus /etc/corpus
sudo chown corpus:corpus /srv/corpus
sudo apt-get update
sudo apt-get install -y python3.12-venv python3-pip awscli

# Clone + venv + deps (still as root, then chown to corpus)
sudo git clone https://github.com/caseywdunn/corpus.git /srv/corpus/repo
sudo python3.12 -m venv /srv/corpus/venv
sudo /srv/corpus/venv/bin/pip install -r /srv/corpus/repo/requirements.txt
sudo chown -R corpus:corpus /srv/corpus

# Bearer token — generate on your laptop, scp over, or:
sudo python3 -c 'import secrets; print(secrets.token_urlsafe(32))' | \
    sudo tee /etc/corpus/mcp.token >/dev/null
sudo chown root:corpus /etc/corpus/mcp.token
sudo chmod 640 /etc/corpus/mcp.token    # readable by the corpus group

# systemd unit
sudo cp /srv/corpus/repo/deploy/corpus-mcp.service \
    /etc/systemd/system/corpus-mcp.service
sudo systemctl daemon-reload
sudo systemctl enable corpus-mcp
# Don't start yet — there's no bundle to serve.
```

Distribute the contents of `/etc/corpus/mcp.token` to your ~20
collaborators via Slack DM, 1Password, signal, whatever.

### 3. First bundle deploy (from your laptop)

Bouchet pipeline produces `output/`.  `rsync` it down or run the
pipeline locally on a smaller corpus.  Then:

```bash
# From your local corpus repo checkout
BUCKET=corpus-bundles deploy/sync_to_s3.sh ~/corpus-output v1.0.0
```

This runs `package_for_serve.py` to distill, then `aws s3 sync` to
upload.

### 4. Pull it onto EC2 and start the service

```bash
ssh ec2-host
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.0.0
sudo systemctl start corpus-mcp
sudo systemctl status corpus-mcp
```

The service should come up in a few seconds and start listening on
127.0.0.1:8080.  ALB/CloudFront terminates TLS on 443 and forwards.

## Updating to a new bundle

Whenever new papers process, new indices are built (e.g., the
geography layer when it lands), or the biblio authority is rebuilt:

```bash
# 1. Generate + upload the new bundle from your laptop
deploy/sync_to_s3.sh ~/corpus-output v1.1.0

# 2. On EC2
ssh ec2-host
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.1.0
```

The update script is idempotent; re-running with the same version
is safe.  Collaborators see the new `bundle_version` via the
`bundle_info` MCP tool the next time they call it.

## Rollback

Every version stays on disk under `/srv/corpus/bundles/`.  To revert:

```bash
ssh ec2-host
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.0.0
```

Seconds.

## Code-only updates

If the change is pure code (bug fix, new MCP tool) without needing
a new bundle:

```bash
ssh ec2-host
cd /srv/corpus/repo
sudo -u corpus git pull
sudo systemctl restart corpus-mcp
```

## Observability

```bash
sudo journalctl -u corpus-mcp -f           # live logs
sudo systemctl status corpus-mcp           # status + last few lines
sudo ss -ltnp | grep 8080                  # confirm the listener
curl -H "Authorization: Bearer $(sudo cat /etc/corpus/mcp.token)" \
    -N http://127.0.0.1:8080/sse           # auth GET from the host
```

## Removing a retired version from S3

```bash
aws s3 rm --recursive s3://corpus-bundles/v0.9.0/
# Versioning keeps a delete marker; purge that too if disk matters.
```

The local copy in `/srv/corpus/bundles/` can be cleaned up with a
simple cron; not urgent since the total disk footprint at 5–10 GB
per version stays manageable for a long time on 10 GB EBS.
