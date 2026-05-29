"""HF implicit-token warning suppressed at pipeline entry (#97).

Importing the pipeline package must set HF_HUB_DISABLE_IMPLICIT_TOKEN
before any submodule pulls in huggingface_hub, so the three model-load
sites (embeddings.py, vision.py) don't emit the implicit-token warning
for the public models we fetch. setdefault must also leave an explicit
operator override intact.
"""
from __future__ import annotations

import importlib
import os
import subprocess
import sys


def test_importing_pipeline_sets_disable_implicit_token():
    importlib.import_module("pipeline")  # import side effect under test
    assert os.environ.get("HF_HUB_DISABLE_IMPLICIT_TOKEN") == "1"


def test_explicit_override_is_preserved():
    # In a clean subprocess with the var pre-set to 0, the setdefault in
    # pipeline/__init__ must not clobber the operator's choice.
    env = dict(os.environ, HF_HUB_DISABLE_IMPLICIT_TOKEN="0")
    out = subprocess.run(
        [sys.executable, "-c",
         "import pipeline, os; print(os.environ['HF_HUB_DISABLE_IMPLICIT_TOKEN'])"],
        env=env, capture_output=True, text=True, check=True,
    )
    assert out.stdout.strip() == "0"
