#!/usr/bin/env python3
"""
Plot <overlap> vs T from MC_simple_result.dat (yyoshiyy 互換).

<overlap> = 時間平均 (1/N) S(t)·S(0). エラーバーは n_runs 回実行の標準誤差.
"""

import math
import os
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

N_RUNS = int(os.environ.get("N_RUNS", "20"))
BASE_SEED = int(os.environ.get("BASE_SEED", "10000"))


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


def parse_result_dat(path):
    """Extract T and overlap from MC_simple_result.dat."""
    temps, overlaps = [], []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 6:
            temps.append(float(parts[0]))
            overlaps.append(float(parts[5]))
    return temps, overlaps


def run_and_collect(label, param_file, n_runs, seed_offset):
    all_temps, all_overlaps = [], []
    for run in range(n_runs):
        seed = BASE_SEED + seed_offset + run * 1007
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
            raise RuntimeError(f"{label} run {run} failed: {proc.stderr}")
        t, ov = parse_result_dat(SCRIPT_DIR / "MC_simple_result.dat")
        all_temps.append(t)
        all_overlaps.append(ov)
    return all_temps, all_overlaps


def mean_stderr(vals):
    n = len(vals)
    if n == 0:
        return 0.0, 0.0
    m = sum(vals) / float(n)
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in vals) / float(n - 1)
    return m, math.sqrt(var / float(n))


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
    make_temp_param(PARAM_BASE, tmp_normal, exchange_interval=0)
    make_temp_param(PARAM_BASE, tmp_exchange, exchange_interval=1)

    try:
        print(f"Running Normal MC {N_RUNS} times...")
        t_n, ov_n = run_and_collect("normal", tmp_normal, N_RUNS, 0)
        print(f"Running Exchange MC {N_RUNS} times...")
        t_e, ov_e = run_and_collect("exchange", tmp_exchange, N_RUNS, 100000)
    finally:
        for p in (tmp_normal, tmp_exchange):
            if p.exists():
                p.unlink()

    temps = t_n[0]
    if t_e[0] != temps:
        print("Error: temperature grids differ")
        return 1

    mean_n = [mean_stderr([ov_n[r][i] for r in range(N_RUNS)])[0]
              for i in range(len(temps))]
    se_n = [mean_stderr([ov_n[r][i] for r in range(N_RUNS)])[1]
            for i in range(len(temps))]
    mean_e = [mean_stderr([ov_e[r][i] for r in range(N_RUNS)])[0]
              for i in range(len(temps))]
    se_e = [mean_stderr([ov_e[r][i] for r in range(N_RUNS)])[1]
            for i in range(len(temps))]

    plt.figure(figsize=(7.0, 4.6))
    plt.errorbar(
        temps, mean_n, yerr=se_n, fmt="o-", capsize=3, color="#1f77b4",
        label=f"Normal (n_runs={N_RUNS})",
    )
    plt.errorbar(
        temps, mean_e, yerr=se_e, fmt="s-", capsize=3, color="#d62728",
        label=f"Exchange (n_runs={N_RUNS})",
    )
    plt.axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.xlabel("Temperature T")
    plt.ylabel("<overlap>")
    plt.title("<overlap> vs T (stderr from n_runs)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_FILE, dpi=170)
    plt.close()

    lines = [
        "<overlap> vs T (time-averaged S(t)·S(0))",
        f"n_runs={N_RUNS}",
        f"base_param={PARAM_BASE}",
        "T  overlap_normal_mean  overlap_normal_stderr  overlap_exchange_mean  overlap_exchange_stderr",
    ]
    for i, T in enumerate(temps):
        lines.append(
            f"{T:.6f}\t{mean_n[i]:.10f}\t{se_n[i]:.10f}\t{mean_e[i]:.10f}\t{se_e[i]:.10f}"
        )
    TXT_FILE.write_text("\n".join(lines) + "\n")

    print(f"Wrote {FIG_FILE}")
    print(f"Wrote {TXT_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
