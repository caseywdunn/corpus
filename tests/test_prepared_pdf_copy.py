"""Preparation resume over frozen/read-only PDF copies."""
import shutil
import stat

import pytest

from pipeline import scan


def test_prepare_replaces_readonly_pdf_without_changing_source(tmp_path, monkeypatch):
    from pipeline import native_text_recovery
    monkeypatch.setattr(native_text_recovery, 'inspect_native_text_regions',
                        lambda *_: {'candidate_count': 0})
    source = tmp_path / 'source.pdf'
    output = tmp_path / 'processed.pdf'
    source.write_bytes(b'new source bytes')
    source.chmod(0o444)
    output.write_bytes(b'previous prepared bytes')
    output.chmod(0o444)
    prior_inode = output.stat().st_ino
    scan.prepare_pdf(source, {'needs_ocr': False}, output)
    assert output.read_bytes() == source.read_bytes() == b'new source bytes'
    assert output.stat().st_ino != prior_inode
    assert stat.S_IMODE(source.stat().st_mode) == 0o444
    # The copied mode remains read-only; subsequent preparation still works.
    scan.prepare_pdf(source, {'needs_ocr': False}, output)
    assert output.read_bytes() == b'new source bytes'
    assert not list(tmp_path.glob('.prepared-*'))


def test_failed_preparation_copy_preserves_previous_output(tmp_path, monkeypatch):
    source = tmp_path / 'source.pdf'
    output = tmp_path / 'processed.pdf'
    source.write_bytes(b'new')
    output.write_bytes(b'previous complete PDF')
    def fail(_source, temporary):
        temporary.write_bytes(b'partial')
        raise OSError('disk full')
    monkeypatch.setattr(scan.shutil, 'copy2', fail)
    with pytest.raises(OSError, match='disk full'):
        scan._copy_prepared_pdf(source, output)
    assert output.read_bytes() == b'previous complete PDF'
    assert not list(tmp_path.glob('.prepared-*'))


def test_preparation_replaces_symlink_without_overwriting_its_target(tmp_path):
    source = tmp_path / 'source.pdf'
    other = tmp_path / 'other.pdf'
    output = tmp_path / 'processed.pdf'
    source.write_bytes(b'new')
    other.write_bytes(b'untouched')
    output.symlink_to(other)
    scan._copy_prepared_pdf(source, output)
    assert not output.is_symlink()
    assert output.read_bytes() == b'new'
    assert other.read_bytes() == b'untouched'
    with pytest.raises(shutil.SameFileError):
        scan._copy_prepared_pdf(source, source)
    assert source.read_bytes() == b'new'
