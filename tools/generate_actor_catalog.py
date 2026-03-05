#!/usr/bin/env python3
"""
Generate LaTeX actor catalog tables for the overview paper.

Primary data source: `knowledge/fuse_knowledge_base.json` (structured actor descriptions + I/O paths).
Completeness guard: scans `src/actors/**/*.jl` for `struct Actor*` definitions and merges in any actors
missing from the knowledge base. Actors present in the knowledge base but not in source are omitted.

Output: `actor_catalog.tex` at repo root.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
KB_PATH = ROOT / "knowledge" / "fuse_knowledge_base.json"
OUT_PATH = ROOT / "actor_catalog.tex"


def latex_escape(text: str) -> str:
    # Minimal escaping for table content.
    return (
        text.replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("$", r"\$")
        .replace("#", r"\#")
        .replace("_", r"\_")
        .replace("{", r"\{")
        .replace("}", r"\}")
        .replace("~", r"\textasciitilde{}")
        .replace("^", r"\textasciicircum{}")
    )


def first_sentence(text: str, max_len: int = 140) -> str:
    text = " ".join((text or "").split())
    if not text:
        return ""
    m = re.search(r"\.(\s|$)", text)
    sent = text[: m.start() + 1] if m else text
    if len(sent) > max_len:
        sent = sent[: max_len - 3].rstrip() + "..."
    return sent


def top_level_ids(paths: Iterable[str]) -> list[str]:
    ids: set[str] = set()
    for p in paths or []:
        p = (p or "").strip()
        if not p:
            continue
        if p.startswith("dd."):
            p = p[3:]
        if p.startswith("dd") and (len(p) == 2 or p[2] in ".[["):
            p = p[2:]
            if p.startswith("."):
                p = p[1:]
        token = re.split(r"[.\[]", p, 1)[0].strip()
        if token:
            ids.add(token)
    return sorted(ids)


def tt_join(items: Iterable[str]) -> str:
    items = list(items)
    if not items:
        return ""
    return ", ".join([rf"\texttt{{{latex_escape(i)}}}" for i in items])


def is_generic_model_desc(desc: str) -> bool:
    d = (desc or "").strip()
    if not d:
        return False
    dl = d.lower()
    # Strong signals that this switch selects between multiple implementations.
    if "options are" in dl or "switch" in dl or "selection" in dl:
        return True
    if "solver to use" in dl or "actor model" in dl:
        return True
    if "replay" in dl or "none" in dl:
        return True
    # Parenthetical list of multiple options (avoid default-only hints).
    m = re.search(r"\(([^)]*)\)", d)
    if m:
        inside = m.group(1)
        if "," in inside and "default" not in inside.lower():
            return True
    return False


@dataclass(frozen=True)
class ActorRow:
    name: str
    category: str
    purpose: str
    reads: list[str]
    writes: list[str]
    notes: str
    kind: str
    model_opts: str


def scan_actor_definitions() -> dict[str, Path]:
    actor_to_file: dict[str, Path] = {}
    for path in (ROOT / "src" / "actors").rglob("*.jl"):
        text = path.read_text(errors="ignore")
        for m in re.finditer(r"\b(?:mutable\s+)?struct\s+(Actor[A-Za-z0-9_]+)\b", text):
            name = m.group(1)
            actor_to_file.setdefault(name, path)
    return actor_to_file


def scan_extension_actor_names() -> set[str]:
    # Identify actors defined in files only included by Julia package extensions.
    extension_files: set[Path] = set()
    for ext_path in (ROOT / "ext").glob("*.jl"):
        text = ext_path.read_text(errors="ignore")
        for m in re.finditer(r'include\(\s*joinpath\(([^)]*)\)\s*\)', text):
            parts = re.findall(r'"([^"]+)"', m.group(1))
            if len(parts) >= 4 and parts[:3] == ["..", "src", "actors"]:
                extension_files.add(ROOT / "src" / "actors" / Path(*parts[3:]))
    extension_actors: set[str] = set()
    for path in extension_files:
        if not path.exists():
            continue
        text = path.read_text(errors="ignore")
        for m in re.finditer(r"\b(?:mutable\s+)?struct\s+(Actor[A-Za-z0-9_]+)\b", text):
            extension_actors.add(m.group(1))
    return extension_actors


def infer_category(name: str, path: Path) -> str:
    # Mirror the folder organization under src/actors where possible.
    try:
        rel = path.relative_to(ROOT / "src" / "actors")
    except ValueError:
        return "uncategorized"
    if len(rel.parts) >= 2:
        return rel.parts[0]
    return "uncategorized"


def find_docstring_purpose(name: str, path: Path) -> str:
    text = path.read_text(errors="ignore")
    # Try to find the docstring that documents the main `ActorX(dd, act; ...)` entrypoint.
    # This pattern is common across the repo and yields good one-line purposes.
    pat = re.compile(
        r'"""(.*?)"""\s*function\s+'
        + re.escape(name)
        + r"\s*\(\s*dd::IMAS\.dd\s*,\s*act::ParametersAllActors",
        re.S,
    )
    m = pat.search(text)
    if m:
        raw = m.group(1).strip()
        # Many docstrings start with a signature line like `ActorX(...)`; drop it.
        lines = raw.splitlines()
        if lines and name in lines[0] and "(" in lines[0]:
            raw = "\n".join(lines[1:]).lstrip()
        return first_sentence(raw)
    return ""


def heuristic_rw_from_source(path: Path) -> tuple[list[str], list[str]]:
    text = path.read_text(errors="ignore")
    touched = set(re.findall(r"\bdd\.([A-Za-z0-9_]+)\b", text))
    if not touched:
        return [], []

    writes: set[str] = set()
    # Common mutating patterns.
    for fn in ["empty", "resize", "push", "append", "fill", "setindex"]:
        if fn == "setindex":
            rx = re.compile(r"\bsetindex!\(\s*dd\.([A-Za-z0-9_]+)\b")
        else:
            rx = re.compile(rf"\b{fn}!\(\s*dd\.([A-Za-z0-9_]+)\b")
        writes.update(rx.findall(text))

    # Mutating calls where the first arg is a dd field, e.g. IMAS.magnetics!(dd.magnetics, ...).
    writes.update(re.findall(r"\b\w+!\(\s*dd\.([A-Za-z0-9_]+)\b", text))

    reads = touched - writes
    return sorted(reads), sorted(writes)


def actor_calls_other_actors(name: str, path: Path) -> bool:
    text = path.read_text(errors="ignore")
    called = set(re.findall(r"\b(Actor[A-Za-z0-9_]+)\s*\(", text))
    called.discard(name)
    # Ignore mentions in docstrings.
    return bool(called)


def main() -> None:
    kb = json.loads(KB_PATH.read_text())
    kb_actors: dict[str, dict] = kb["actors"]

    actor_to_file = scan_actor_definitions()
    extension_actors = scan_extension_actor_names()

    # Authoritative inventory: actors that exist in source.
    actor_names = sorted(actor_to_file.keys())

    # Manual corrections for actors missing from the knowledge base.
    # (The KB was extracted from `*_actor.jl` files; some actors live in other files.)
    manual = {
        "ActorControllerIp": {
            "category": "control",
            "purpose": "Controls the loop voltage to track the plasma current setpoint.",
            "reads": ["global_time", "pulse_schedule", "core_profiles", "controllers"],
            "writes": ["controllers"],
        },
        "ActorMagnetics": {
            "category": "diagnostics",
            "purpose": "Calculates synthetic magnetic diagnostics from equilibrium data.",
            "reads": ["equilibrium"],
            "writes": ["magnetics"],
        },
        "ActorInterferometer": {
            "category": "diagnostics",
            "purpose": "Calculates synthetic interferometer measurements from equilibrium and core profiles.",
            "reads": ["equilibrium", "core_profiles", "interferometer"],
            "writes": ["interferometer"],
        },
        "ActorSawteethSource": {
            "category": "hcd",
            "purpose": "Applies sawtooth reconnection effects to core sources when q < 1.",
            "reads": ["global_time", "sawteeth", "equilibrium", "core_profiles"],
            "writes": ["core_sources"],
        },
        "ActorThermalSystemModels": {
            "category": "balance_plant",
            "purpose": "Detailed thermal plant model using ThermalSystemModels.jl (extension).",
            "reads": ["balance_of_plant"],
            "writes": ["balance_of_plant"],
        },
    }

    rows: list[ActorRow] = []
    categories: dict[str, list[str]] = {}
    model_map: list[tuple[str, str]] = []

    for name in actor_names:
        path = actor_to_file[name]
        kb_entry = kb_actors.get(name)

        if kb_entry is not None:
            category = kb_entry.get("category") or infer_category(name, path)
            purpose = first_sentence(kb_entry.get("description", ""))
            reads = top_level_ids(kb_entry.get("data_inputs", []))
            writes = top_level_ids(kb_entry.get("data_outputs", []))
            model_desc = str((kb_entry.get("key_parameters") or {}).get("model") or "")
        else:
            category = infer_category(name, path)
            purpose = find_docstring_purpose(name, path)
            reads, writes = heuristic_rw_from_source(path)
            model_desc = ""

        if name in manual:
            category = manual[name]["category"]
            reads = manual[name]["reads"]
            writes = manual[name]["writes"]
            purpose = manual[name].get("purpose", purpose)

        if not purpose:
            purpose = first_sentence(f"{name} actor.")

        # Kind: Generic if it selects among multiple models; Compound if it instantiates other actors; else Specific.
        kind = "Specific"
        model_opts = ""
        if model_desc and is_generic_model_desc(model_desc):
            kind = "Generic"
            m = re.search(r"\(([^)]{3,})\)", model_desc)
            if m:
                model_opts = m.group(1).strip().strip(".")
            else:
                m = re.search(r"options\s+are\s+(.*)$", model_desc, re.I)
                if m:
                    model_opts = m.group(1).strip().strip(".")
            if not model_opts:
                syms = sorted(set(re.findall(r":[A-Za-z0-9_]+", model_desc)))
                model_opts = ", ".join(syms) if syms else model_desc.strip().strip(".")
            model_map.append((name, model_opts))
        elif actor_calls_other_actors(name, path):
            kind = "Compound"

        notes = r"\textit{extension}" if name in extension_actors else ""

        rows.append(
            ActorRow(
                name=name,
                category=category,
                purpose=purpose,
                reads=reads,
                writes=writes,
                notes=notes,
                kind=kind,
                model_opts=model_opts,
            )
        )
        categories.setdefault(category, []).append(name)

    # Preferred category order (paper-friendly).
    preferred_order = [
        "equilibrium",
        "stability",
        "pedestal",
        "transport",
        "current",
        "hcd",
        "ec",
        "nbi",
        "control",
        "pf",
        "build",
        "nuclear",
        "wall_loading",
        "divertors",
        "sol",
        "balance_plant",
        "costing",
        "diagnostics",
        "compound",
        "uncategorized",
    ]
    cat_order = {c: i for i, c in enumerate(preferred_order)}

    title_map = {
        "hcd": "HCD",
        "ec": "EC",
        "nbi": "NBI",
        "pf": "PF",
        "sol": "SOL",
        "balance_plant": "Balance of Plant",
        "wall_loading": "Wall Loading",
    }

    out: list[str] = []
    out.append("% Auto-generated by tools/generate_actor_catalog.py")
    out.append(r"% Requires: \usepackage{longtable} and (recommended) \usepackage{booktabs}.")
    out.append("")
    out.append(r"\section{Actor catalog}\label{appendix:actors}")
    out.append(
        r"\noindent This catalog lists FUSE actors grouped by domain. "
        r"The \texttt{Reads} and \texttt{Writes} columns show the first-level \texttt{dd} entries "
        r"(IMAS IDSs) touched by each actor."
    )
    out.append("")

    out.append(r"\subsection{Generic actor model selection}")
    out.append(r"\begin{longtable}{p{0.22\linewidth}p{0.68\linewidth}}")
    out.append(r"\toprule")
    out.append(r"\textbf{Generic actor} & \textbf{Available models (\texttt{act.\ldots.model})}\\")
    out.append(r"\midrule")
    out.append(r"\endhead")
    for name, opts in sorted(model_map, key=lambda t: t[0]):
        out.append(rf"\texttt{{{latex_escape(name)}}} & {latex_escape(opts)}\\")
    out.append(r"\bottomrule")
    out.append(r"\end{longtable}")
    out.append("")

    rows_by_name = {r.name: r for r in rows}
    for cat in sorted(categories.keys(), key=lambda c: (cat_order.get(c, 999), c)):
        names = sorted(set(categories[cat]))
        title = title_map.get(cat, cat.replace("_", " ").title())
        out.append(rf"\subsection{{{latex_escape(title)}}}")
        out.append(
            r"\begin{longtable}{p{0.18\linewidth}p{0.09\linewidth}p{0.33\linewidth}p{0.16\linewidth}p{0.16\linewidth}p{0.06\linewidth}}"
        )
        out.append(r"\toprule")
        out.append(r"\textbf{Actor} & \textbf{Kind} & \textbf{Purpose} & \textbf{Reads} & \textbf{Writes} & \textbf{Notes}\\")
        out.append(r"\midrule")
        out.append(r"\endhead")
        for name in names:
            r = rows_by_name[name]
            out.append(
                rf"\texttt{{{latex_escape(r.name)}}} & {latex_escape(r.kind)} & {latex_escape(r.purpose)} & "
                rf"{tt_join(r.reads)} & {tt_join(r.writes)} & {r.notes}\\"
            )
        out.append(r"\bottomrule")
        out.append(r"\end{longtable}")
        out.append("")

    OUT_PATH.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
