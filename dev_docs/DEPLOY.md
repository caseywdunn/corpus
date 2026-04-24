# AWS deployment runbook

End-to-end CLI-only deploy — no AWS console clicks.  The stack is
defined declaratively in [deploy/stack.yaml](../deploy/stack.yaml)
(CloudFormation); on-host setup and bundle operations happen via
`ssh` + the shell scripts in [deploy/](../deploy/).

## Architecture

Deliberately boring, per PLAN.md §10, scaled to ~20 collaborators:

- **EC2 `t3.small`** on the default VPC, Ubuntu 24.04 LTS, Elastic IP.
- **nginx** on the instance terminates TLS via Let's Encrypt, reverse-
  proxies to `mcp_server.py` on `127.0.0.1:8080`.
- **S3** bucket holds versioned serve bundles; the instance reads
  from it via IAM role (no keys).
- **Your DNS** (any provider) has an A record pointing at the
  Elastic IP.

Deliberately not yet: CloudFront, ALB, WAF, Route 53.  Upgrade if you
hit rate-limit problems or scale beyond one instance.

Layout on the EC2 host:

```text
/srv/corpus/
  repo/               git checkout of this repository
  venv/               Python 3.12 virtualenv
  bundles/
    v1.0.0/           copy of s3://<bucket>/v1.0.0/
    v1.1.0/
  bundle -> bundles/v1.1.0/       atomic symlink to active version

/etc/corpus/
  mcp.token           bearer token, mode 640, root:corpus

/etc/nginx/sites-available/corpus-mcp   reverse-proxy config
/etc/systemd/system/corpus-mcp.service  unit from deploy/
```

## Prerequisites on your laptop

```bash
# aws cli v2
brew install awscli                   # macOS; other OSes similar

# Authenticate with an AWS account that can create EC2 / IAM / S3.
aws configure                         # ACCESS_KEY, SECRET_KEY, region
aws sts get-caller-identity           # sanity check
```

You'll also need:

- A domain under your control (any registrar / DNS host works).
- An email address for Let's Encrypt renewal notices.

## 1. Pick parameters

Set these once per shell session so subsequent commands can reference them:

```bash
export REGION=us-east-1                             # or your preferred
export STACK=corpus-mcp
export BUCKET=corpus-bundles-$(whoami)              # globally unique
export KEYPAIR=corpus-mcp
export DOMAIN=corpus.example.edu                    # your hostname
export ACME_EMAIL=you@example.edu
export MY_IP=$(curl -4 -s https://ifconfig.co)/32   # SSH source CIDR
```

## 2. Create the SSH keypair

```bash
aws ec2 create-key-pair \
    --region "$REGION" \
    --key-name "$KEYPAIR" \
    --query KeyMaterial --output text > ~/.ssh/${KEYPAIR}.pem
chmod 600 ~/.ssh/${KEYPAIR}.pem
```

## 3. Deploy the CloudFormation stack

```bash
aws cloudformation deploy \
    --region "$REGION" \
    --stack-name "$STACK" \
    --template-file deploy/stack.yaml \
    --capabilities CAPABILITY_IAM \
    --parameter-overrides \
        BucketName=$BUCKET \
        KeyPairName=$KEYPAIR \
        SSHAllowCIDR=$MY_IP
```

Takes ~2 minutes.  When it's done, grab the outputs:

```bash
aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query 'Stacks[0].Outputs[*].[OutputKey,OutputValue]' --output table

export EIP=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query 'Stacks[0].Outputs[?OutputKey==`PublicIp`].OutputValue' --output text)
echo "Elastic IP: $EIP"
```

## 4. Point DNS at the Elastic IP

**Not automated** — depends on your DNS provider.  Create an A record:

```text
corpus.example.edu.  A  <EIP from step 3>
```

Wait for propagation (`dig +short $DOMAIN` should return `$EIP`).

## 5. On-host setup — one-time

SSH in using the keypair from step 2:

```bash
ssh -i ~/.ssh/${KEYPAIR}.pem ubuntu@$EIP
# (cloud-init from the template has already installed python3.12,
# nginx, certbot, awscli, and created the `corpus` user.)

# Clone the repo + create the venv.
sudo git clone https://github.com/caseywdunn/corpus.git /srv/corpus/repo
sudo python3.12 -m venv /srv/corpus/venv
sudo /srv/corpus/venv/bin/pip install -r /srv/corpus/repo/requirements.txt
sudo chown -R corpus:corpus /srv/corpus

# Generate the bearer token and stash at /etc/corpus/mcp.token.
sudo bash -c 'python3 -c "import secrets; print(secrets.token_urlsafe(32))" \
    > /etc/corpus/mcp.token'
sudo chown root:corpus /etc/corpus/mcp.token
sudo chmod 640 /etc/corpus/mcp.token
sudo cat /etc/corpus/mcp.token    # copy for distribution to collaborators

# Install the systemd unit.  Don't start yet — there's no bundle.
sudo cp /srv/corpus/repo/deploy/corpus-mcp.service \
    /etc/systemd/system/corpus-mcp.service
sudo systemctl daemon-reload
sudo systemctl enable corpus-mcp

# nginx reverse-proxy: install the config + get a TLS cert.
# The DNS A record from step 4 must be resolving before certbot runs.
sudo cp /srv/corpus/repo/deploy/nginx.conf /etc/nginx/sites-available/corpus-mcp
sudo sed -i "s/REPLACE.example.edu/$DOMAIN/g" /etc/nginx/sites-available/corpus-mcp
# (on the instance, $DOMAIN isn't your local shell var — hardcode the hostname)
sudo ln -sf /etc/nginx/sites-available/corpus-mcp /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default
sudo nginx -t
sudo systemctl reload nginx

sudo certbot --nginx -d corpus.example.edu --agree-tos --email you@example.edu \
    --redirect --non-interactive
# certbot edits the nginx config to wire in the cert + schedules renewals.
```

## 6. First bundle — from your laptop

Still on your laptop:

```bash
# The output dir can be a local pipeline run or an rsync'd Bouchet tree.
BUCKET=$BUCKET deploy/sync_to_s3.sh ~/corpus-output v1.0.0
```

## 7. Pull it on EC2 and start the service

```bash
ssh -i ~/.ssh/${KEYPAIR}.pem ubuntu@$EIP
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.0.0
sudo systemctl start corpus-mcp
sudo systemctl status corpus-mcp
```

## 8. Smoke-test end-to-end

From your laptop:

```bash
# Unauth should be 401
curl -s -o /dev/null -w '%{http_code}\n' https://corpus.example.edu/sse
# → 401

# Authed should be 200 + text/event-stream
curl -s -o /dev/null -w '%{http_code} %{content_type}\n' \
    -H "Authorization: Bearer <the-token-you-copied-above>" \
    https://corpus.example.edu/sse
# → 200 text/event-stream; ...
```

Then add the server to Claude Code for a full client check:

```bash
claude mcp add corpus-prod https://corpus.example.edu/sse \
    --transport sse --scope user \
    --header "Authorization: Bearer <token>"
```

## Updating to a new bundle

Whenever new data lands (new papers, geography layer, rebuilt biblio
authority):

```bash
# On your laptop
deploy/sync_to_s3.sh ~/corpus-output v1.1.0

# On EC2
ssh ubuntu@$EIP
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.1.0
```

Collaborators see the new `bundle_version` the next time their
client calls the `bundle_info` tool.

## Rollback

Every version stays on disk under `/srv/corpus/bundles/`:

```bash
ssh ubuntu@$EIP
sudo -u corpus /srv/corpus/repo/deploy/update.sh v1.0.0
```

Seconds.  Rollback from S3 alone (if /srv/corpus was wiped) takes
one extra `aws s3 sync`.

## Code-only updates

```bash
ssh ubuntu@$EIP
cd /srv/corpus/repo
sudo -u corpus git pull
sudo systemctl restart corpus-mcp
```

## Observability

```bash
sudo journalctl -u corpus-mcp -f           # live server logs
sudo journalctl -u nginx -f                # reverse-proxy logs
sudo systemctl status corpus-mcp           # status + last few lines
sudo ss -ltnp | grep -E '8080|443'         # confirm listeners
curl -sS -H "Authorization: Bearer $(sudo cat /etc/corpus/mcp.token)" \
    -N https://localhost/sse --resolve localhost:443:127.0.0.1 -k
```

## Teardown

When you're done or rebuilding from scratch:

```bash
# Empty the bucket (CloudFormation retains it on purpose).
aws s3 rm --recursive s3://$BUCKET/
aws s3api delete-bucket --bucket $BUCKET --region $REGION

# Delete the stack (EC2, EIP, IAM, security group).
aws cloudformation delete-stack --stack-name $STACK --region $REGION
aws cloudformation wait stack-delete-complete --stack-name $STACK --region $REGION
```

## Gotchas

- **Default VPC required** by the template as shipped.  If your
  account has the default VPC deleted, add `VpcId` + `SubnetId`
  parameters to [stack.yaml](../deploy/stack.yaml) and pass them
  through.
- **certbot fails the first time** if DNS hasn't propagated yet
  (step 4).  Re-run it after `dig` confirms the A record.
- **SSH times out** if your home IP changed since stack deploy.
  Update `SSHAllowCIDR` and re-deploy:
  `aws cloudformation deploy --stack-name $STACK --template-file deploy/stack.yaml --parameter-overrides SSHAllowCIDR=$(curl -4 -s https://ifconfig.co)/32 BucketName=$BUCKET KeyPairName=$KEYPAIR --capabilities CAPABILITY_IAM`
- **ACM in us-east-1** isn't relevant here (we use Let's Encrypt on
  the instance) — but if you later add CloudFront, that's the region
  its cert must live in.
