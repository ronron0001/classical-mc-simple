#!/usr/bin/env python3
"""
Compute q(t_max) vs T using variance across internal Sample bins directly.

Flow:
  1) Read one base param file: samples/square_L16_Ising/param.def.
  2) Create temporary params by changing only exchange_interval:
       - Normal: exchange_interval=0
       - Exchange: exchange_interval=1
  3) Run each method once.
  3) Read MC_simple_q_tmax_samples.dat from each run.
  4) Estimate mean and stderr over Sample bins.
  5) Output q_tmax_vs_T.png and q_tmax_vs_T_summary.txt.
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

NORMAL_FILE = SCRIPT_DIR / "q_tmax_samples_normal.dat"
EXCHANGE_FILE = SCRIPT_DIR / "q_tmax_samples_exchange.dat"
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


def run_once(label, param_file, seed, out_file):
    env = os.environ.copy()
    env["MC_SIMPLE_SEED"] = str(seed)
    log_path = SCRIPT_DIR / f"{label}_q_tmax.log"
    cmd = [str(MC_BIN), str(param_file), str(LATTICE), str(INTERACTION)]
    with open(log_path, "w") as logf:
        proc = subprocess.run(
            cmd,
            cwd=SCRIPT_DIR,
            env=env,
            stdout=logf,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"{label} run failed. See {log_path}")

    src = SCRIPT_DIR / "MC_simple_q_tmax_samples.dat"
    if not src.exists():
        raise RuntimeError(f"{label}: MC_simple_q_tmax_samples.dat not found")
    shutil.copy2(src, out_file)


def parse_sample_q_file(path):
    temps = None
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# sample"):
                tokens = line.split()[2:]
                temps = [float(tok.replace("q_T", "")) for tok in tokens]
            continue
        parts = line.split()
        rows.append([float(x) for x in parts[1:]])
    if temps is None:
        raise RuntimeError(f"failed to parse header in {path}")
    if not rows:
        raise RuntimeError(f"no sample rows in {path}")
    for r in rows:
        if len(r) != len(temps):
            raise RuntimeError(f"column mismatch in {path}")
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


def summarize(rows):
    n_sample = len(rows)
    n_temp = len(rows[0])
    means = []
    ses = []
    for t_idx in range(n_temp):
        vals = [rows[s][t_idx] for s in range(n_sample)]
        m, se = mean_and_stderr(vals)
        means.append(m)
        ses.append(se)
    return means, ses, n_sample


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

    t_n, rows_n = parse_sample_q_file(NORMAL_FILE)
    t_e, rows_e = parse_sample_q_file(EXCHANGE_FILE)
    if t_n != t_e:
        print("Error: temperature grids differ between normal and exchange")
        return 1
    temps = t_n

    mean_n, se_n, n_n = summarize(rows_n)
    mean_e, se_e, n_e = summarize(rows_e)

    plt.figure(figsize=(7.0, 4.6))
    plt.errorbar(
        temps,
        mean_n,
        yerr=se_n,
        fmt="o-",
        capsize=3,
        color="#1f77b4",
        label=f"Normal (Sample bins: n={n_n})",
    )
    plt.errorbar(
        temps,
        mean_e,
        yerr=se_e,
        fmt="s-",
        capsize=3,
        color="#d62728",
        label=f"Exchange (Sample bins: n={n_e})",
    )
    plt.axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.xlabel("Temperature T")
    plt.ylabel("q(t_max)")
    plt.title("q(t_max) vs T (stderr from Sample-bin variance)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_FILE, dpi=170)
    plt.close()

    lines = []
    lines.append("q(t=t_max) vs T comparison")
    lines.append("stderr estimated from per-Sample variance (within one run)")
    lines.append(f"base_param={PARAM_BASE}")
    lines.append(f"n_normal={n_n}, n_exchange={n_e}")
    lines.append(
        "Columns: T, q_normal_mean, q_normal_stderr, q_exchange_mean, q_exchange_stderr"
    )
    for i, T in enumerate(temps):
        lines.append(
            f"{T:.6f}\t{mean_n[i]:.10f}\t{se_n[i]:.10f}\t{mean_e[i]:.10f}\t{se_e[i]:.10f}"
        )
    TXT_FILE.write_text("\n".join(lines) + "\n")

    print(f"Wrote {FIG_FILE}")
    print(f"Wrote {TXT_FILE}")
    print(f"Saved raw sample bins: {NORMAL_FILE}, {EXCHANGE_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
