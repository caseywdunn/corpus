"""Shared `rich` console layer for the `corpus` CLI (#63).

A single :data:`console` instance backs every operator-facing
subcommand. ``console.is_terminal`` gates the rich rendering: a TTY
operator gets emoji + color + bars + spinners; SLURM ``.out``
capture, CI logs, and ``journalctl`` readers see clean ASCII —
uniform call sites, no per-site ``if tty:`` branches.

Critically, the diagnostic ``logging`` stream stays plain text. Rich
is the operator-facing layer (interactive run, success summary,
``corpus check`` report, ``corpus status --report``); the structured
event log that lands in journals + grep pipelines stays untouched.

Pattern matches sharkmer's ``indicatif`` usage: one ``show_progress``
bool, hidden Progress when not a TTY.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator, Literal, Optional

from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

# Single Console shared across every `corpus` subcommand.
# Auto-detects TTY; force-disable via env (``NO_COLOR``) is honored
# by rich natively.
console: Console = Console()


_Status = Literal["ok", "fail", "warn", "info"]


def _symbol(status: _Status) -> str:
    """Emoji status with ASCII fallback when not on a TTY."""
    if console.is_terminal:
        return {"ok": "✓", "fail": "✗", "warn": "⚠", "info": "•"}[status]
    return {"ok": "[ok]", "fail": "[fail]", "warn": "[warn]", "info": "[info]"}[status]


_STATUS_STYLE = {
    "ok": "bold green",
    "fail": "bold red",
    "warn": "bold yellow",
    "info": "bold cyan",
}


def print_status(message: str, status: _Status = "info") -> None:
    """Print ``message`` prefixed with the configured status symbol.

    On TTY: colored emoji + styled message. Off TTY: plain ASCII tag.

    ``message`` is treated as literal text in both modes. Square
    brackets in the body (e.g. ``[build bundle] ...``) are escaped
    before being handed to rich so they aren't interpreted as markup
    tags and silently dropped from terminal output.
    """
    from rich.markup import escape as _rich_escape

    sym = _symbol(status)
    if console.is_terminal:
        console.print(
            f"[{_STATUS_STYLE[status]}]{sym}[/] {_rich_escape(message)}"
        )
    else:
        # `[ok]` etc. would be interpreted as rich markup tags and
        # stripped; print plain.
        print(f"{sym} {message}")


@contextmanager
def progress(
    description: str,
    total: Optional[int] = None,
    *,
    transient: bool = False,
) -> Iterator["Progress"]:
    """Context manager that yields a :class:`Progress` configured for TTY/non-TTY.

    Caller adds tasks via the yielded instance. When not a TTY, the
    bars are hidden (``disable=True``) so SLURM ``.out`` files don't
    fill with carriage returns; the underlying tasks still complete.

    ``total=None`` produces a spinner; ``total=N`` produces a bar.
    """
    columns: list = [
        TextColumn("[progress.description]{task.description}"),
    ]
    if total is None:
        columns.append(SpinnerColumn(spinner_name="dots"))
    else:
        columns += [BarColumn(), TaskProgressColumn()]
    columns.append(TimeElapsedColumn())

    with Progress(
        *columns,
        console=console,
        transient=transient,
        disable=not console.is_terminal,
    ) as p:
        p.add_task(description, total=total)
        yield p


__all__ = [
    "console",
    "print_status",
    "progress",
]
