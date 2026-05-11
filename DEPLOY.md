# AWS deployment runbook

End-to-end CLI-only deploy — no AWS console clicks required.  The EC2
stack is defined in [deploy/stack.yaml](deploy/stack.yaml) (CloudFormation);
on-host setup and bundle operations happen via `ssh` + the shell scripts
in [deploy/](deploy/).

## Architecture

```
Client → ALB (TLS, ACM wildcard *.siphonophores.org)
           └─ Host: <organism>.siphonophores.org
              → <organism>-mcp-tg → EC2:80 → nginx → 127.0.0.1:8080
```

- **ALB `siphonophores-mcp-alb`** — shared across all organism MCP servers.
  TLS terminates here (ACM cert); EC2 only listens on :80.  One listener
  rule per organism (host-header match → target group).
- **EC2 `t3.large`** on the default VPC, Ubuntu 24.04 LTS.
- **nginx** on the instance is ALB-aware: no TLS, no certbot.  It serves
  `GET /health` directly (for ALB health checks) and proxies `/sse`,
  `/messages/`, `/healthz` to `mcpsrv.main` on `127.0.0.1:8080`.
- **S3** bucket holds versioned serve bundles; the instance reads from it
  via IAM instance role (no credentials).
- **DNS** (`siphonophores.org` at DreamHost) — one CNAME per organism
  pointing at the ALB DNS name.  DNS never changes after first deploy;
  cutover is done by swapping target group registrations.

**Shared infrastructure** (one-time, already exists for siphonophores.org):

| Resource | Value |
|---|---|
| ALB | `siphonophores-mcp-alb` |
| ALB SG | `mcp-alb-sg` (TCP 443+80 from 0.0.0.0/0) |
| ACM cert ARN | `arn:aws:acm:us-east-1:656748891540:certificate/23a1234d-e908-4bcb-8a49-c8977b5df9a0` |
| ALB DNS name | `siphonophores-mcp-alb-570029510.us-east-1.elb.amazonaws.com` |
| Region | `us-east-1` |
| VPC | `vpc-0428ad714156dbb54` |

Layout on the EC2 host:

```text
/srv/corpus/
  repo/               git checkout (dev branch)
  venv/               Python 3.12 virtualenv
  bundles/
    v0.4.0/           copy of s3://<bucket>/v0.4.0/
    v0.5.0/
  bundle -> bundles/v0.5.0/       atomic symlink to active version

/etc/corpus/
  mcp.token           bearer token, mode 640, root:corpus
  update.conf         BUCKET=<name> for deploy/update.sh

/etc/nginx/sites-available/corpus-mcp   reverse-proxy config
/etc/systemd/system/corpus-mcp.service  unit from deploy/
```

---

## Quick remote deploy (any host)

The minimum to bring up an SSE server on any host you already have — no AWS, no CloudFormation. The full AWS runbook starts in §"Prerequisites" below.

**1. Generate a bearer token.** Strong URL-safe random; mode-600 file; never pass on the CLI (would leak via `ps`):

```bash
python -c 'import secrets; print(secrets.token_urlsafe(32))' > ~/corpus-mcp.token
chmod 600 ~/corpus-mcp.token
```

**2. Start the server over SSE.**

```bash
corpus serve --output-dir <bundle> -- \
    --transport sse --host 127.0.0.1 --port 8080 \
    --auth-token-file ~/corpus-mcp.token
```

Clients send `Authorization: Bearer <token>` on every request. Without `--auth-token-file` or `CORPUS_MCP_TOKEN` the server runs open and logs a loud warning — fine for localhost experiments, never safe for public-facing deploys.

**3. Smoke-test before exposing to anyone else.** [tools/smoke_test_sse.py](tools/smoke_test_sse.py) launches its own server on a free port, generates its own token, and drives the full stack — 401 without auth, 200 with auth, MCP `initialize` → `list_tools` → `bundle_info` → `list_papers`:

```bash
python tools/smoke_test_sse.py <bundle>
```

All seven checks should pass locally before you move to a real server.

**4. Connect a client.** Remote-MCP support exists natively in Claude Desktop, claude.ai web, and Claude Code. Recommended in order:

- **Claude Code (CLI)** — tested and working against our server:

  ```bash
  claude mcp add corpus-remote http://127.0.0.1:18080/sse \
      --transport sse \
      --scope user \
      --header "Authorization: Bearer $(cat ~/corpus-mcp.token)"
  ```

  Use `--scope project` to tie to one repo, `claude mcp remove corpus-remote --scope user` to undo. `/mcp` in any session lists connected servers and their tools.

- **Claude Desktop / claude.ai web (Custom Connectors UI)** — Settings → Connectors → "Add custom connector", paste the URL. Available on Free / Pro / Max / Team / Enterprise plans ([Anthropic docs](https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp)). **Caveat:** Custom Connectors target the **Streamable HTTP** transport with OAuth-style auth, while our server speaks **SSE** with a static bearer token. Whether the UI accepts our URL depends on how strict the client is. If it rejects the SSE endpoint, fall back to the bridge below.

- **Claude Desktop fallback — mcp-remote bridge.** Edit `~/Library/Application Support/Claude/claude_desktop_config.json` to launch [`mcp-remote`](https://www.npmjs.com/package/mcp-remote) as a stdio subprocess that bridges to the SSE endpoint:

  ```json
  {
    "mcpServers": {
      "corpus-remote": {
        "command": "/opt/homebrew/bin/npx",
        "args": ["-y", "mcp-remote",
                 "http://127.0.0.1:18080/sse",
                 "--header", "Authorization: Bearer <token>"]
      }
    }
  }
  ```

A cleaner Custom Connectors experience would add a **Streamable HTTP** transport (`mcp.streamable_http_app()` in FastMCP) with OAuth discovery, but this is deferred indefinitely — SSE + bearer token works for the ~20-collaborator deploy target.

---

## Prerequisites

```bash
# aws cli v2
brew install awscli          # macOS; other OSes similar
aws configure                # or set AWS_PROFILE
aws sts get-caller-identity  # sanity check
```

SSH key `corpus-mcp` already exists in us-east-1 and at
`~/.ssh/corpus-mcp.pem` (on Bouchet: `/nfs/roberts/home/cwd7/.ssh/corpus-mcp.pem`).

---

## Deploying a new organism (first time)

### 1. Pick parameters

```bash
export REGION=us-east-1
export ORGANISM=siphonophores          # lowercase, no spaces
export DOMAIN=${ORGANISM}.siphonophores.org
export BUCKET=${ORGANISM}              # S3 bucket name
export STACK=${ORGANISM}-mcp
export KEYPAIR=corpus-mcp
export VPC_ID=vpc-0428ad714156dbb54
export SUBNET_ID=subnet-06be250021470525f   # any public subnet in the VPC
export MY_IP=$(curl -4 -s https://ifconfig.co)/32
export ALB_ARN=arn:aws:elasticloadbalancing:us-east-1:656748891540:loadbalancer/app/siphonophores-mcp-alb/107a154e420a47f8
export ACM_ARN=arn:aws:acm:us-east-1:656748891540:certificate/23a1234d-e908-4bcb-8a49-c8977b5df9a0
```

### 2. Deploy CloudFormation stack

First deploy creates the S3 bucket (default `CreateBucket=true`):

```bash
aws cloudformation deploy \
    --region $REGION \
    --stack-name $STACK \
    --template-file deploy/stack.yaml \
    --capabilities CAPABILITY_IAM \
    --parameter-overrides \
        BucketName=$BUCKET \
        KeyPairName=$KEYPAIR \
        SSHAllowCIDR=$MY_IP \
        VpcId=$VPC_ID \
        SubnetId=$SUBNET_ID

export EC2_IP=$(aws cloudformation describe-stacks --stack-name $STACK --region $REGION \
    --query 'Stacks[0].Outputs[?OutputKey==`PublicIp`].OutputValue' --output text)
export EC2_SG=$(aws cloudformation describe-stack-resources --stack-name $STACK --region $REGION \
    --query 'StackResources[?ResourceType==`AWS::EC2::SecurityGroup`].PhysicalResourceId' \
    --output text)
echo "EC2 IP: $EC2_IP  SG: $EC2_SG"
```

### 3. Wire up ALB

```bash
# Restrict EC2 port-80 to ALB SG only (default rule is 0.0.0.0/0)
ALB_SG_ID=$(aws ec2 describe-security-groups --region $REGION \
    --filters Name=group-name,Values=mcp-alb-sg \
    --query 'SecurityGroups[0].GroupId' --output text)

aws ec2 revoke-security-group-ingress --region $REGION --group-id $EC2_SG \
    --protocol tcp --port 80 --cidr 0.0.0.0/0
aws ec2 authorize-security-group-ingress --region $REGION --group-id $EC2_SG \
    --ip-permissions "[{\"IpProtocol\":\"tcp\",\"FromPort\":80,\"ToPort\":80,
      \"UserIdGroupPairs\":[{\"GroupId\":\"$ALB_SG_ID\",\"Description\":\"HTTP from mcp-alb-sg only\"}]}]"

# Create target group
TG_ARN=$(aws elbv2 create-target-group \
    --region $REGION \
    --name ${ORGANISM}-mcp-tg \
    --protocol HTTP --port 80 \
    --vpc-id $VPC_ID \
    --target-type instance \
    --health-check-path /health \
    --healthy-threshold-count 2 \
    --health-check-interval-seconds 30 \
    --query 'TargetGroups[0].TargetGroupArn' --output text)
echo "Target group: $TG_ARN"

# Register EC2 in target group
INSTANCE_ID=$(aws cloudformation describe-stacks --stack-name $STACK --region $REGION \
    --query 'Stacks[0].Outputs[?OutputKey==`InstanceId`].OutputValue' --output text)
aws elbv2 register-targets --region $REGION \
    --target-group-arn $TG_ARN --targets Id=$INSTANCE_ID

# Add host-header listener rule to the shared ALB (HTTPS:443)
HTTPS_LISTENER_ARN=$(aws elbv2 describe-listeners --region $REGION \
    --load-balancer-arn $ALB_ARN \
    --query 'Listeners[?Port==`443`].ListenerArn' --output text)

# Pick a priority not already in use:
aws elbv2 describe-rules --listener-arn $HTTPS_LISTENER_ARN --region $REGION \
    --query 'Rules[*].Priority' --output text

aws elbv2 create-rule \
    --region $REGION \
    --listener-arn $HTTPS_LISTENER_ARN \
    --priority 2 \
    --conditions "[{\"Field\":\"host-header\",\"HostHeaderConfig\":{\"Values\":[\"$DOMAIN\"]}}]" \
    --actions "[{\"Type\":\"forward\",\"TargetGroupArn\":\"$TG_ARN\"}]"
```

### 4. DNS (DreamHost)

In the DreamHost panel → `siphonophores.org` DNS → Custom Records → Add:
- **Name**: `<organism>` (e.g. `corpus`)
- **Type**: CNAME
- **Value**: `siphonophores-mcp-alb-570029510.us-east-1.elb.amazonaws.com`

Wait for propagation: `dig $DOMAIN CNAME +short`

### 5. On-host setup

```bash
ssh -i ~/.ssh/${KEYPAIR}.pem ubuntu@$EC2_IP
cloud-init status --wait    # must print 'status: done'

# Repo + Python venv
sudo git clone -b dev https://github.com/caseywdunn/corpus.git /srv/corpus/repo
sudo python3.12 -m venv /srv/corpus/venv
sudo /srv/corpus/venv/bin/pip install -e /srv/corpus/repo
sudo chown -R corpus:corpus /srv/corpus

# Bearer token
sudo mkdir -p /etc/corpus
sudo bash -c 'python3 -c "import secrets; print(secrets.token_urlsafe(32))" > /etc/corpus/mcp.token'
sudo chown root:corpus /etc/corpus/mcp.token && sudo chmod 640 /etc/corpus/mcp.token
sudo cat /etc/corpus/mcp.token    # ← copy this; share with collaborators

# Systemd unit
sudo cp /srv/corpus/repo/deploy/corpus-mcp.service /etc/systemd/system/corpus-mcp.service
sudo systemctl daemon-reload && sudo systemctl enable corpus-mcp

# S3 bucket config
echo "BUCKET=<bucket-name>" | sudo tee /etc/corpus/update.conf
sudo chmod 644 /etc/corpus/update.conf

# Sudoers
sudo tee /etc/sudoers.d/corpus-mcp > /dev/null <<'EOF'
corpus ALL=(root) NOPASSWD: /bin/systemctl restart corpus-mcp, /bin/systemctl is-active corpus-mcp
EOF
sudo chmod 440 /etc/sudoers.d/corpus-mcp
sudo visudo -c -f /etc/sudoers.d/corpus-mcp

# nginx (ALB-aware config — no TLS, no certbot)
sudo cp /srv/corpus/repo/deploy/nginx.conf /etc/nginx/sites-available/corpus-mcp
sudo sed -i "s/REPLACE.example.edu/<your-domain>/g" /etc/nginx/sites-available/corpus-mcp
sudo ln -sf /etc/nginx/sites-available/corpus-mcp /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default
sudo nginx -t && sudo systemctl reload nginx

# Pull bundle from S3 and start
sudo -u corpus /srv/corpus/repo/deploy/update.sh <version>
sudo systemctl start corpus-mcp
sudo systemctl status corpus-mcp
```

### 6. Smoke test

```bash
# Direct to EC2 (bypasses ALB — no TLS)
curl -s http://$EC2_IP/health
# → ok

# Via ALB + DNS (after propagation)
TOKEN=<bearer-token>
curl -s -o /dev/null -w '%{http_code}\n' https://$DOMAIN/sse
# → 401
curl -s -o /dev/null -w '%{http_code} %{content_type}\n' \
    -H "Authorization: Bearer $TOKEN" https://$DOMAIN/sse
# → 200 text/event-stream

# Add to Claude Code
claude mcp add ${ORGANISM}-prod https://$DOMAIN/sse \
    --transport sse --scope user \
    --header "Authorization: Bearer $TOKEN"
```

---

## Replacing an EC2 (zero-DNS-change cutover)

DNS always points at the ALB; cutover is done by swapping target
group registrations.  The old EC2 keeps running until you're satisfied
and can be accessed directly via its IP for verification at any time.

### 1. Deploy new stack

```bash
# CreateBucket=false — bucket already exists in the old stack
aws cloudformation deploy \
    --region $REGION \
    --stack-name ${ORGANISM}-mcp-v2 \
    --template-file deploy/stack.yaml \
    --capabilities CAPABILITY_IAM \
    --parameter-overrides \
        BucketName=$BUCKET \
        CreateBucket=false \
        KeyPairName=$KEYPAIR \
        SSHAllowCIDR=$MY_IP \
        VpcId=$VPC_ID \
        SubnetId=$SUBNET_ID

NEW_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name ${ORGANISM}-mcp-v2 --region $REGION \
    --query 'Stacks[0].Outputs[?OutputKey==`InstanceId`].OutputValue' --output text)
NEW_EC2_IP=$(aws cloudformation describe-stacks \
    --stack-name ${ORGANISM}-mcp-v2 --region $REGION \
    --query 'Stacks[0].Outputs[?OutputKey==`PublicIp`].OutputValue' --output text)
```

Restrict the new EC2 SG port-80 to `mcp-alb-sg` (same as step 3 above).

### 2. On-host setup (same as §5 above)

Use `$NEW_EC2_IP` for SSH.  Smoke-test directly:
```bash
curl -s http://$NEW_EC2_IP/health   # → ok
```

### 3. Cut over

```bash
TG_ARN=$(aws elbv2 describe-target-groups --region $REGION \
    --names ${ORGANISM}-mcp-tg \
    --query 'TargetGroups[0].TargetGroupArn' --output text)

OLD_INSTANCE_ID=<old-instance-id>

# Register new instance — ALB will route to both once new one is healthy
aws elbv2 register-targets --region $REGION \
    --target-group-arn $TG_ARN --targets Id=$NEW_INSTANCE_ID

# Wait for new instance to pass health checks
aws elbv2 wait target-in-service --region $REGION \
    --target-group-arn $TG_ARN --targets Id=$NEW_INSTANCE_ID
echo "New instance healthy — deregistering old"

# Remove old instance from rotation (it keeps running; traffic stops)
aws elbv2 deregister-targets --region $REGION \
    --target-group-arn $TG_ARN --targets Id=$OLD_INSTANCE_ID
```

The old EC2 is still accessible at its public IP for verification.
When satisfied, delete the old stack:

```bash
aws cloudformation delete-stack --stack-name ${ORGANISM}-mcp --region $REGION
aws cloudformation wait stack-delete-complete --stack-name ${ORGANISM}-mcp --region $REGION
# Rename new stack if desired (not directly possible in CF — just document the new name)
```

---

## Updating to a new bundle (no EC2 replacement)

```bash
# On Bouchet — upload new bundle
BUCKET=$BUCKET deploy/sync_to_s3.sh /path/to/pipeline/output v0.5.0

# On EC2
ssh -i ~/.ssh/${KEYPAIR}.pem ubuntu@$EC2_IP
sudo -u corpus /srv/corpus/repo/deploy/update.sh v0.5.0
# (reads BUCKET from /etc/corpus/update.conf)
```

Rollback to any retained version is the same command with the old version string.
Every version stays on disk under `/srv/corpus/bundles/` until manually pruned.

---

## Code-only updates

```bash
ssh -i ~/.ssh/${KEYPAIR}.pem ubuntu@$EC2_IP
sudo git config --global --add safe.directory /srv/corpus/repo
sudo git -C /srv/corpus/repo pull origin dev
sudo systemctl restart corpus-mcp
```

---

## Observability

```bash
sudo journalctl -u corpus-mcp -f        # live server logs
sudo journalctl -u nginx -f             # reverse-proxy logs
sudo systemctl status corpus-mcp        # status + last few lines
curl -s http://localhost/health         # nginx health (no host header needed)
curl -s -H "Host: $DOMAIN" http://localhost/healthz   # MCP liveness
```

---

## Teardown

```bash
# Remove from ALB target group first
aws elbv2 deregister-targets --region $REGION \
    --target-group-arn $TG_ARN --targets Id=$INSTANCE_ID
aws elbv2 delete-target-group --region $REGION --target-group-arn $TG_ARN

# Remove listener rule
aws elbv2 delete-rule --region $REGION --rule-arn <rule-arn>

# Delete CF stack (EC2, IAM, security group — bucket is retained)
aws cloudformation delete-stack --stack-name $STACK --region $REGION
aws cloudformation wait stack-delete-complete --stack-name $STACK --region $REGION

# Only if you want to destroy the bucket and all bundles:
aws s3 rm --recursive s3://$BUCKET/
aws s3api delete-bucket --bucket $BUCKET --region $REGION
```

---

## Gotchas

- **ALB idle timeout must be 3600s** for SSE connections.  The default
  is 60s, which silently kills long-running SSE streams mid-session.
  Already set on `siphonophores-mcp-alb`; if you create a new ALB set it
  immediately: `aws elbv2 modify-load-balancer-attributes --attributes Key=idle_timeout.timeout_seconds,Value=3600`.

- **FastMCP 421 Invalid Host header.**  FastMCP's SSE transport rejects
  any `Host` that isn't `127.0.0.1:*` or `localhost:*` (DNS-rebinding
  protection).  `deploy/nginx.conf` rewrites `Host: localhost:8080` before
  proxying.  If you rebuild nginx config from scratch, keep the rewrite —
  passing `$host` upstream silently breaks every client with a 421.

- **Health check path is `/health`, not `/healthz`.**  nginx serves
  `/health` directly (no auth, returns `ok`); `/healthz` is proxied to
  the MCP server and requires the correct Host header.  ALB health checks
  must use `/health`.

- **`CreateBucket=false` for replacement deploys.**  When the S3 bucket
  already exists (owned by a previous stack that's still running), pass
  `CreateBucket=false` or CloudFormation will fail trying to create a
  duplicate bucket.  The IAM policy always grants access by bucket name
  regardless of this flag.

- **SSH CIDR changes between Bouchet login nodes.**  Bouchet routes
  outbound through different IPs (`192.31.2.x/24`).  The CF stack sets
  `SSHAllowCIDR=192.31.2.0/24` to cover the whole block.  If SSH times
  out, check your current IP (`curl -4 -s https://ifconfig.co`) and
  update the SG rule.

- **git safe.directory on EC2.**  When running `sudo git` on files owned
  by another user, git may refuse with "dubious ownership".  Fix:
  `sudo git config --global --add safe.directory /srv/corpus/repo`

- **Clone the `dev` branch**, not `main`.  The ALB-aware nginx config,
  bearer-auth transport, and taxa path scrubber are on `dev` and not yet
  merged to `main`.  Use `git clone -b dev ...` or `git checkout dev`
  after cloning.

- **EBS full** during `update.sh` shows as `[Errno 28] No space left`.
  Check `df -h /`.  Resize live without reboot:

  ```bash
  VOLUME_ID=$(aws ec2 describe-volumes \
      --filters "Name=attachment.instance-id,Values=$INSTANCE_ID" \
      --region $REGION --query 'Volumes[0].VolumeId' --output text)
  aws ec2 modify-volume --volume-id $VOLUME_ID --size 60 --region $REGION
  until [[ "$(aws ec2 describe-volumes-modifications \
      --volume-ids $VOLUME_ID --region $REGION \
      --query 'VolumesModifications[0].ModificationState' --output text)" \
      =~ ^(optimizing|completed)$ ]]; do sleep 10; done

  # On EC2:
  sudo growpart /dev/nvme0n1 1
  sudo resize2fs /dev/nvme0n1p1
  df -h /
  ```

- **SCP denies `ec2:AllocateAddress`.**  The template uses the instance's
  auto-assigned public IP (no Elastic IP).  If you stop/start the instance
  the IP changes — but the ALB routes by target group registration, not IP,
  so the public URL is unaffected.  Only SSH access needs the new IP.
