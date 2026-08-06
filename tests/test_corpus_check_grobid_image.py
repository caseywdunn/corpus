"""Tests for the Grobid image-drift warning in `corpus check`.

``/api/isalive`` answering only proves *something* is listening on the
port. Because the compose service carries ``restart: unless-stopped``, a
container started from an earlier ``docker-compose.yml`` survives pulls
and keeps serving that port — so the image can silently disagree with
the compose file, and the image is what decides whether references are
parsed by the CRF or the DeLFT model.

``pipeline.cli._grobid_image_warning`` holds the branch table as a pure
function (running, expected, apple_silicon) -> Optional[str], so these
tests pin the contract without a Docker daemon. The IO around it
(``_grobid_container_image`` / ``_compose_grobid_image``) is tested via
mocks for the cases that must stay silent.
"""
from __future__ import annotations

from unittest import mock

from pipeline.cli import (
    _check_grobid_image,
    _grobid_container_image,
    _grobid_image_warning,
)

LFOPPIANO = "lfoppiano/grobid:0.8.1"
DELFT = "grobid/grobid:0.8.1"


def test_silent_when_images_agree():
    assert _grobid_image_warning(LFOPPIANO, LFOPPIANO, apple_silicon=True) is None


def test_silent_when_running_image_unknown():
    """Grobid under Apptainer or a bare JVM has no container to inspect."""
    assert _grobid_image_warning(None, LFOPPIANO, apple_silicon=False) is None


def test_silent_when_no_compose_expectation():
    """A non-clone install has no docker-compose.yml to compare against."""
    assert _grobid_image_warning(DELFT, None, apple_silicon=True) is None


def test_warns_on_drift():
    msg = _grobid_image_warning(DELFT, LFOPPIANO, apple_silicon=False)
    assert msg is not None
    assert DELFT in msg and LFOPPIANO in msg
    assert "docker compose up -d grobid" in msg


def test_drift_warning_is_generic_off_apple_silicon():
    """The AVX/Rosetta note is Apple-Silicon-specific; DeLFT is a
    legitimate opt-in on an AVX-capable Linux x86_64 host."""
    msg = _grobid_image_warning(DELFT, LFOPPIANO, apple_silicon=False)
    assert "AVX" not in msg


def test_delft_on_apple_silicon_gets_the_platform_note():
    msg = _grobid_image_warning(DELFT, LFOPPIANO, apple_silicon=True)
    assert "AVX" in msg
    assert "Apple Silicon" in msg


def test_unrelated_drift_on_apple_silicon_stays_generic():
    msg = _grobid_image_warning("someone/grobid:0.7", LFOPPIANO, apple_silicon=True)
    assert msg is not None
    assert "AVX" not in msg


def test_container_lookup_skips_remote_grobid():
    """On a cluster grobid.url names a compute node — not ours to inspect."""
    with mock.patch("shutil.which") as which:
        assert _grobid_container_image("http://c03n05.cluster:8070") is None
        which.assert_not_called()


def test_container_lookup_skips_without_docker():
    with mock.patch("shutil.which", return_value=None):
        assert _grobid_container_image("http://localhost:8070") is None


def test_container_lookup_returns_single_match():
    proc = mock.Mock(returncode=0, stdout=f"{DELFT}\n")
    with mock.patch("shutil.which", return_value="/usr/bin/docker"), \
         mock.patch("subprocess.run", return_value=proc):
        assert _grobid_container_image("http://127.0.0.1:8070") == DELFT


def test_container_lookup_is_silent_when_ambiguous():
    """Two containers publishing the port — we can't say which answered."""
    proc = mock.Mock(returncode=0, stdout=f"{DELFT}\n{LFOPPIANO}\n")
    with mock.patch("shutil.which", return_value="/usr/bin/docker"), \
         mock.patch("subprocess.run", return_value=proc):
        assert _grobid_container_image("http://localhost:8070") is None


def test_container_lookup_survives_a_broken_docker():
    """A docker CLI present but daemon down must not fail the check."""
    with mock.patch("shutil.which", return_value="/usr/bin/docker"), \
         mock.patch("subprocess.run", side_effect=OSError("boom")):
        assert _grobid_container_image("http://localhost:8070") is None


def test_check_wires_the_pieces_together():
    with mock.patch("pipeline.cli._grobid_container_image", return_value=DELFT), \
         mock.patch("pipeline.cli._compose_grobid_image", return_value=LFOPPIANO):
        msg = _check_grobid_image("http://localhost:8070")
    assert msg is not None and DELFT in msg


def test_compose_image_matches_the_shipped_compose_file():
    """Guards the parse against a compose-file restructure."""
    from pipeline.cli import _compose_grobid_image
    image = _compose_grobid_image()
    assert image is not None, "docker-compose.yml should be readable from a clone"
    assert "grobid" in image
