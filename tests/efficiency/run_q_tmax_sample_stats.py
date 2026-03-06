#!/usr/bin/env python3
"""
Plot <overlap> vs T from MC_simple_overlap_samples.dat.

<overlap> = 時間平均 (1/N) S(t)·S(0). 平均・分散は Sample ビンで計算.
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

FIG_FILE = SCRIPT_DIR / "q_tmax_vs_T.png"
TXT_FILE = SCRIPT_DIR / "q_tmax_vs_T_summary.txt"

SEED_NORMAL = int(os.environ.get("SEED_NORMAL", "500000"))
SEED_EXCHANGE = int(os.environ.get("SEED_EXCHANGE", "600000"))


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


def parse_overlap_samples(path):
    """Parse MC_simple_overlap_samples.dat. Returns (temps, rows) where rows = (sample, step, qvals)."""
    temps = None
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            if line.startswith("# sample  step"):
                toks = line.split()[3:]
                temps = [float(tok.replace("q_T", "")) for tok in toks]
            continue
        p = line.split()
        s, step = int(p[0]), int(p[1])
        qvals = [float(x) for x in p[2:]]
        if temps is None or len(qvals) != len(temps):
            raise RuntimeError(f"invalid overlap row in {path}")
        rows.append((s, step, qvals))
    if temps is None or not rows:
        raise RuntimeError(f"failed to parse {path}")
    return temps, rows


def per_sample_time_avg_overlap(temps, rows):
    """Compute per-sample time-averaged overlap for each T. step 1..Total_Step (measurement phase)."""
    sample_ids = sorted({s for s, _, _ in rows})
    steps = sorted({st for _, st, _ in rows})
    total_step = max(steps)
    # measurement phase: step 1 .. Total_Step (step 0 = initial)
    meas_steps = [st for st in steps if 1 <= st <= total_step]
    if not meas_steps:
        meas_steps = steps[1:] if len(steps) > 1 else []

    # grid[step][temp] = list of (sample_idx, q)
    st_to_i = {st: i for i, st in enumerate(steps)}
    s_to_i = {s: i for i, s in enumerate(sample_ids)}
    # sample_overlap[sample_idx][temp] = sum over meas steps of q
    n_sample = len(sample_ids)
    n_temp = len(temps)
    sums = [[0.0 for _ in range(n_temp)] for _ in range(n_sample)]

    for s, st, qvals in rows:
        if st < 1 or st > total_step:
            continue
        si = s_to_i[s]
        for t_idx, q in enumerate(qvals):
            sums[si][t_idx] += q

    n_meas = len(meas_steps)
    overlap = [[sums[si][t] / n_meas if n_meas > 0 else 0.0
                for t in range(n_temp)] for si in range(n_sample)]
    return temps, overlap, n_sample


def mean_stderr(vals):
    n = len(vals)
    if n == 0:
        return 0.0, 0.0
    m = sum(vals) / float(n)
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in vals) / float(n - 1)
    return m, math.sqrt(var / float(n))


def run_once(label, param_file, seed, out_file):
    env = os.environ.copy()
    env["MC_SIMPLE_SEED"] = str(seed)
    proc = subprocess.run(
        [str(MC_BIN), str(param_file), str(LATTICE), str(INTERACTION)],
        cwd=SCRIPT_DIR,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{label} run failed: {proc.stderr}")
    src = SCRIPT_DIR / "MC_simple_overlap_samples.dat"
    if not src.exists():
        raise RuntimeError(f"{label}: MC_simple_overlap_samples.dat not found")
    shutil.copy2(src, out_file)


def main():
    if not MC_BIN.exists():
        print("Error: MC_simple not found. Build first.")
        return 1
    for fp in (PARAM_BASE, LATTICE, INTERACTION):
        if not fp.exists():
            print(f"Error: missing {fp}")
            return 1

    tmp_normal = SCRIPT_DIR / "_tmp_param_normal.def"
    tmp_exchange = SCRIPT_DIR / "_tmp_param_exchange.def"
    normal_file = SCRIPT_DIR / "q_tmax_overlap_normal.dat"
    exchange_file = SCRIPT_DIR / "q_tmax_overlap_exchange.dat"
    make_temp_param(PARAM_BASE, tmp_normal, exchange_interval=0)
    make_temp_param(PARAM_BASE, tmp_exchange, exchange_interval=1)

    try:
        print("Running Normal MC (1 run)...")
        run_once("normal", tmp_normal, SEED_NORMAL, normal_file)
        print("Running Exchange MC (1 run)...")
        run_once("exchange", tmp_exchange, SEED_EXCHANGE, exchange_file)
    finally:
        for p in (tmp_normal, tmp_exchange):
            if p.exists():
                p.unlink()

    temps_n, overlap_n, n_n = per_sample_time_avg_overlap(*parse_overlap_samples(normal_file))
    temps_e, overlap_e, n_e = per_sample_time_avg_overlap(*parse_overlap_samples(exchange_file))

    if temps_n != temps_e:
        print("Error: temperature grids differ")
        return 1

    temps = temps_n
    mean_n = [mean_stderr([overlap_n[s][t] for s in range(n_n)])[0] for t in range(len(temps))]
    se_n = [mean_stderr([overlap_n[s][t] for s in range(n_n)])[1] for t in range(len(temps))]
    mean_e = [mean_stderr([overlap_e[s][t] for s in range(n_e)])[0] for t in range(len(temps))]
    se_e = [mean_stderr([overlap_e[s][t] for s in range(n_e)])[1] for t in range(len(temps))]

    plt.figure(figsize=(7.0, 4.6))
    plt.errorbar(
        temps, mean_n, yerr=se_n, fmt="o-", capsize=3, color="#1f77b4",
        label=f"Normal (Sample={n_n})",
    )
    plt.errorbar(
        temps, mean_e, yerr=se_e, fmt="s-", capsize=3, color="#d62728",
        label=f"Exchange (Sample={n_e})",
    )
    plt.axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.xlabel("Temperature T")
    plt.ylabel("<overlap>")
    plt.title("<overlap> vs T (stderr from Sample bins)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_FILE, dpi=170)
    plt.close()

    lines = [
        "<overlap> vs T (time-averaged S(t)·S(0))",
        "stderr from per-Sample variance (1 run each)",
        f"n_normal={n_n}, n_exchange={n_e}",
        f"base_param={PARAM_BASE}",
        "T  overlap_normal_mean  overlap_normal_stderr  overlap_exchange_mean  overlap_exchange_stderr",
    ]
    for i, T in enumerate(temps):
        lines.append(
            f"{T:.6f}\t{mean_n[i]:.10f}\t{se_n[i]:.10f}\t{mean_e[i]:.10f}\t{se_e[i]:.10f}"
        )
    TXT_FILE.write_text("\n".join(lines) + "\n")

    for f in (normal_file, exchange_file):
        if f.exists():
            f.unlink()

    print(f"Wrote {FIG_FILE}")
    print(f"Wrote {TXT_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
