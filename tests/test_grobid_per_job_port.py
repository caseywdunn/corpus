"""Concurrent Grobid jobs get a port of their own (#279).

Grobid's Dropwizard service binds a fixed 8070, and SLURM is free to
co-schedule several of those jobs onto one node — at which point every
instance after the first dies ~10 s in with a Jetty BindException.
Submitting six pipelines put five Grobid jobs on two nodes and three
failed.

Worse than a plain failure: the job has already reached RUNNING, so a
chain waiting on job state alone points Stage 1 at that node and is
served by *another chain's* server. That happened — two documents in a
build whose own Grobid had died recorded `grobid.outcome = extracted`.

The port pair is derived from the job ID by one shared function, so the
client cannot drift from the server.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SLURM = ROOT / "slurm"
PATHS_SH = SLURM / "bouchet_paths.sh"


def _ports_for(job_ids) -> dict:
    """Run the real shell functions, not a Python re-implementation.

    `bouchet_paths.sh` cannot be sourced off-cluster — it cd's into the
    project root — so extract just the two functions and evaluate those.
    Re-implementing the arithmetic here would let this file agree with
    itself while disagreeing with the scripts, which is the whole thing
    being guarded. One subprocess for the whole batch, so checking every
    slot stays cheap enough to actually do.
    """
    text = PATHS_SH.read_text()
    funcs = re.findall(
        r"^(corpus_grobid_(?:admin_)?port\(\)\s*\{.*?^\})",
        text, re.MULTILINE | re.DOTALL,
    )
    assert len(funcs) == 2, f"expected both port functions, found {len(funcs)}"
    ids = list(job_ids)
    script = "\n".join(funcs) + "\nfor j in " + " ".join(str(i) for i in ids) + (
        '; do echo "$j $(corpus_grobid_port "$j") '
        '$(corpus_grobid_admin_port "$j")"; done\n'
    )
    out = subprocess.run(["bash", "-c", script], capture_output=True,
                         text=True, check=True).stdout.split()
    assert len(out) == 3 * len(ids)
    return {int(out[i]): (int(out[i + 1]), int(out[i + 2]))
            for i in range(0, len(out), 3)}


def _ports(job_id: int) -> tuple:
    return _ports_for([job_id])[job_id]


def test_the_five_job_ids_that_actually_collided_no_longer_do():
    """The real IDs from the failed 2026-09-06 run: two Grobid jobs
    survived and three died on 8070."""
    collided = [25125616, 25125651, 25125701, 25125722, 25125743]
    pairs = [_ports(j) for j in collided]
    used = [p for pair in pairs for p in pair]
    assert len(set(used)) == len(used), f"overlapping ports: {pairs}"


def test_the_admin_connector_moves_with_the_application_one():
    """The trap the issue's own suggested fix would have hit. Dropwizard
    binds an admin connector too — default 8071 — so overriding only the
    application port still dies with the same BindException. Verified
    against lfoppiano/grobid:0.8.1: app-only exits 1 on
    `java.net.BindException: Address already in use`, both-overridden runs
    three instances side by side returning byte-identical TEI."""
    app, admin = _ports(12345)
    assert admin == app + 1


def test_no_two_job_ids_can_produce_overlapping_pairs():
    """A stride of 2 with admin = app + 1 is what makes this true; a
    stride of 1 would put one job's admin port on the next job's
    application port."""
    seen: dict = {}
    derived = _ports_for(range(400))
    for job_id, (app, admin) in sorted(derived.items()):
        for port in (app, admin):
            assert port not in seen, (
                f"job {job_id} reuses port {port} from job {seen[port]}"
            )
            seen[port] = job_id


def test_the_default_grobid_port_is_left_free():
    """8070/8071 stay clear so a hand-started Grobid — the documented
    debugging path, and the local docker-compose one — never collides with
    a batch job."""
    for job_id in (0, 1, 399, 25125616):
        assert 8070 not in _ports(job_id)
        assert 8071 not in _ports(job_id)


def test_the_port_wraps_rather_than_running_away():
    assert _ports(400) == _ports(0)
    assert _ports(0)[0] == 8100


# --- the scripts use the shared function rather than a literal ----------


def test_the_grobid_job_overrides_both_dropwizard_connectors():
    text = (SLURM / "batch_grobid.sh").read_text()
    assert "dw.server.applicationConnectors[0].port" in text
    assert "dw.server.adminConnectors[0].port" in text
    # SINGULARITYENV_ is honoured by both Singularity 3.x and Apptainer,
    # unlike `--env`, whose support depends on the runtime version.
    assert "SINGULARITYENV_JAVA_OPTS" in text


def test_the_grobid_job_derives_its_port():
    text = (SLURM / "batch_grobid.sh").read_text()
    assert 'corpus_grobid_port "${SLURM_JOB_ID:-0}"' in text
    assert 'corpus_grobid_admin_port "${SLURM_JOB_ID:-0}"' in text


def test_the_pipeline_derives_the_same_port_from_the_same_job_id():
    """Not a copy of the arithmetic — the same function, so they cannot
    drift."""
    text = (SLURM / "batch_pipeline.sh").read_text()
    assert 'corpus_grobid_port "$GROBID_JOB"' in text
    assert ':${GROBID_PORT}"' in text


def test_no_script_hardcodes_the_grobid_url_port():
    """A literal `:8070` in a URL is the defect. Prose and the
    localhost/docker-compose default are fine."""
    for script in sorted(SLURM.glob("*.sh")):
        for i, line in enumerate(script.read_text().splitlines(), 1):
            if line.lstrip().startswith("#"):
                continue
            assert not re.search(r'\$\{?GROBID_NODE\}?:8070', line), (
                f"{script.name}:{i} hardcodes the port: {line.strip()}"
            )


def test_the_preflight_check_tests_the_derived_port_not_8070():
    """The `ss` backstop has to look at the port this job will actually
    bind, or it guards nothing."""
    text = (SLURM / "batch_grobid.sh").read_text()
    assert '$GROBID_PORT|$GROBID_ADMIN_PORT' in text


def test_the_bouchet_runbook_does_not_teach_the_fixed_port():
    """BOUCHET.md's manual path is what an operator copies; leaving :8070
    in it would reintroduce the collision by hand."""
    text = (ROOT / "dev_docs" / "BOUCHET.md").read_text()
    assert 'squeue -j "$GROBID_JOB" -h -o %N):8070' not in text
    assert "corpus_grobid_port" in text
