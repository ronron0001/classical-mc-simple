#!/usr/bin/env python3
"""
Create overlap_comparison.png from per-Sample overlap bins directly.

No external N_REPEAT loop:
  - one Normal run (exchange_interval=0)
  - one Exchange run (exchange_interval=1)
stderr is estimated from variance across internal Sample bins.
"""

import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
MC_BIN = PROJECT_ROOT / "src" / "build" / "MC_simple"
SAMPLE_DIR = PROJECT_ROOT / "samples" / "square_L16_Ising"

PARAM_BASE = SAMPLE_DIR / "param.def"
LATTICE = SAMPLE_DIR / "lattice.def"
INTERACTION = SAMPLE_DIR / "interaction.def"

NORMAL_FILE = SCRIPT_DIR / "overlap_samples_normal.dat"
EXCHANGE_FILE = SCRIPT_DIR / "overlap_samples_exchange.dat"
FIG_FILE = SCRIPT_DIR / "overlap_comparison.png"
TXT_FILE = SCRIPT_DIR / "overlap_comparison_summary.txt"

SEED_NORMAL = int(os.environ.get("SEED_NORMAL", "90000"))
SEED_EXCHANGE = int(os.environ.get("SEED_EXCHANGE", "190000"))


def make_temp_param(base_param, out_param, exchange_interval):
    lines = base_param.read_text().splitlines()
    out_lines = []
    found = False
    for raw in lines:
        line = raw.strip()
        if line.startswith("#") or "=" not in line:
            out_lines.append(raw)
            continue
        key, _ = [x.strip() for x in line.split("=", 1)]
        if key == "exchange_interval":
            out_lines.append(f"exchange_interval = {exchange_interval}")
            found = True
        else:
            out_lines.append(raw)
    if not found:
        out_lines.append(f"exchange_interval = {exchange_interval}")
    out_param.write_text("\n".join(out_lines) + "\n")


def run_once(label, param_file, seed, out_file):
    env = os.environ.copy()
    env["MC_SIMPLE_SEED"] = str(seed)
    log_path = SCRIPT_DIR / f"{label}_overlap.log"
    cmd = [str(MC_BIN), str(param_file), str(LATTICE), str(INTERACTION)]
    with open(log_path, "w") as logf:
        proc = subprocess.run(
            cmd, cwd=SCRIPT_DIR, env=env, stdout=logf, stderr=subprocess.STDOUT
        )
    if proc.returncode != 0:
        raise RuntimeError(f"{label} run failed. See {log_path}")
    src = SCRIPT_DIR / "MC_simple_overlap_samples.dat"
    if not src.exists():
        raise RuntimeError(f"{label}: MC_simple_overlap_samples.dat not found")
    shutil.copy2(src, out_file)


def parse_overlap_samples(path):
    temps = None
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# sample  step"):
                toks = line.split()[3:]
                temps = [float(tok.replace("q_T", "")) for tok in toks]
            continue
        p = line.split()
        s = int(p[0])
        step = int(p[1])
        qvals = [float(x) for x in p[2:]]
        if temps is None or len(qvals) != len(temps):
            raise RuntimeError(f"invalid overlap row in {path}")
        rows.append((s, step, qvals))
    if temps is None or not rows:
        raise RuntimeError(f"failed to parse {path}")
    return temps, rows


def mean_and_stderr(vals):
    n = len(vals)
    if n == 0:
        return 0.0, 0.0
    m = sum(vals) / float(n)
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in vals) / float(n - 1)
    return m, math.sqrt(var / float(n))


def summarize_grid(temps, rows):
    sample_ids = sorted({s for s, _, _ in rows})
    steps = sorted({st for _, st, _ in rows})
    s_to_i = {s: i for i, s in enumerate(sample_ids)}
    st_to_i = {st: i for i, st in enumerate(steps)}

    # grid[step][temp] = list over samples
    grid = [[[] for _ in temps] for _ in steps]
    for s, st, q in rows:
        si = s_to_i[s]
        _ = si  # for clarity; sample index implied by append order
        sti = st_to_i[st]
        for t_idx, val in enumerate(q):
            grid[sti][t_idx].append(val)

    mean = [[0.0 for _ in temps] for _ in steps]
    se = [[0.0 for _ in temps] for _ in steps]
    for sti in range(len(steps)):
        for t_idx in range(len(temps)):
            m, e = mean_and_stderr(grid[sti][t_idx])
            mean[sti][t_idx] = m
            se[sti][t_idx] = e
    return steps, mean, se, len(sample_ids)


def integrated_overlap(series):
    return sum(series)


def main():
    if not MC_BIN.exists():
        print("Error: MC_simple not found. Build first.")
        return 1
    for fp in (PARAM_BASE, LATTICE, INTERACTION):
        if not fp.exists():
            print(f"Error: missing file {fp}")
            return 1

    tmp_normal = SCRIPT_DIR / "_tmp_param_normal.def"
    tmp_exchange = SCRIPT_DIR / "_tmp_param_exchange.def"
    make_temp_param(PARAM_BASE, tmp_normal, exchange_interval=0)
    make_temp_param(PARAM_BASE, tmp_exchange, exchange_interval=1)

    try:
        run_once("normal", tmp_normal, SEED_NORMAL, NORMAL_FILE)
        run_once("exchange", tmp_exchange, SEED_EXCHANGE, EXCHANGE_FILE)
    finally:
        if tmp_normal.exists():
            tmp_normal.unlink()
        if tmp_exchange.exists():
            tmp_exchange.unlink()

    t_n, rows_n = parse_overlap_samples(NORMAL_FILE)
    t_e, rows_e = parse_overlap_samples(EXCHANGE_FILE)
    if t_n != t_e:
        print("Error: temperature grids differ between normal and exchange")
        return 1
    temps = t_n

    steps_n, mean_n, se_n, n_n = summarize_grid(temps, rows_n)
    steps_e, mean_e, se_e, n_e = summarize_grid(temps, rows_e)
    if steps_n != steps_e:
        print("Error: step grids differ between normal and exchange")
        return 1
    steps = steps_n

    selected = sorted(set([0, len(temps) // 2, len(temps) - 1]))
    fig, axes = plt.subplots(
        len(selected), 1, figsize=(7.0, 2.6 * len(selected)), sharex=True
    )
    if len(selected) == 1:
        axes = [axes]

    lines = []
    lines.append("Overlap Comparison: Normal MC vs Exchange MC")
    lines.append("stderr estimated from per-Sample variance (within one run)")
    lines.append(f"n_normal={n_n}, n_exchange={n_e}")
    lines.append("")

    for ax, t_idx in zip(axes, selected):
        temp = temps[t_idx]
        qn = [mean_n[s][t_idx] for s in range(len(steps))]
        qe = [mean_e[s][t_idx] for s in range(len(steps))]
        sen = [se_n[s][t_idx] for s in range(len(steps))]
        see = [se_e[s][t_idx] for s in range(len(steps))]

        n_lo = [max(-1.0, qn[i] - sen[i]) for i in range(len(qn))]
        n_hi = [min(1.0, qn[i] + sen[i]) for i in range(len(qn))]
        e_lo = [max(-1.0, qe[i] - see[i]) for i in range(len(qe))]
        e_hi = [min(1.0, qe[i] + see[i]) for i in range(len(qe))]

        ax.plot(steps, qn, color="#1f77b4", label=f"Normal (n={n_n})", linewidth=1.8)
        ax.plot(steps, qe, color="#d62728", label=f"Exchange (n={n_e})", linewidth=1.8)
        ax.fill_between(steps, n_lo, n_hi, color="#1f77b4", alpha=0.2)
        ax.fill_between(steps, e_lo, e_hi, color="#d62728", alpha=0.2)
        ax.set_ylabel(f"q(t)\nT={temp:.2f}")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=9)

        area_n = integrated_overlap(qn)
        area_e = integrated_overlap(qe)
        ratio = area_e / area_n if area_n != 0 else 0.0
        lines.append(
            f"T={temp:.4f}: sum_t q_normal={area_n:.6f}, sum_t q_exchange={area_e:.6f}, ratio(exchange/normal)={ratio:.6f}"
        )

    axes[-1].set_xlabel("MC step t (measurement phase)")
    fig.suptitle("Spin-Configuration Overlap Decay (Sample-bin stderr)")
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(FIG_FILE, dpi=170)
    plt.close(fig)

    TXT_FILE.write_text("\n".join(lines) + "\n")
    print(f"Wrote {FIG_FILE}")
    print(f"Wrote {TXT_FILE}")
    print(f"Saved raw sample bins: {NORMAL_FILE}, {EXCHANGE_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
