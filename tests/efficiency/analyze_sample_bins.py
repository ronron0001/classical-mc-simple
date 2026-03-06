#!/usr/bin/env python3
"""
Analyze MC_simple_sample_bins.dat for Normal vs Exchange.
Error bars are estimated from variance across internal Sample bins.
"""

import math
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
NORMAL_FILE = SCRIPT_DIR / "normal_sample_bins.dat"
EXCHANGE_FILE = SCRIPT_DIR / "exchange_sample_bins.dat"
OUTPUT_FILE = SCRIPT_DIR / "efficiency_comparison.txt"
FIG_OUTPUT = SCRIPT_DIR / "efficiency_comparison.png"

try:
    import matplotlib.pyplot as plt  # type: ignore[import]

    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


def parse_sample_bins(path):
    by_t = {}
    sample_ids = set()
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        if len(p) < 6:
            continue
        s = int(p[0])
        t = float(p[1])
        e = float(p[2])
        c = float(p[3])
        m2 = float(p[4])
        sample_ids.add(s)
        if t not in by_t:
            by_t[t] = {"E": [], "C": [], "M2": []}
        by_t[t]["E"].append(e)
        by_t[t]["C"].append(c)
        by_t[t]["M2"].append(m2)
    return by_t, len(sample_ids)


def mean_and_stderr(vals):
    n = len(vals)
    if n == 0:
        return 0.0, 0.0
    m = sum(vals) / float(n)
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in vals) / float(n - 1)
    return m, math.sqrt(var / float(n))


def main():
    if not NORMAL_FILE.exists() or not EXCHANGE_FILE.exists():
        raise SystemExit("Error: sample bin files not found. Run run_efficiency_test.sh first.")

    normal, n_n = parse_sample_bins(NORMAL_FILE)
    exchange, n_e = parse_sample_bins(EXCHANGE_FILE)
    temps = sorted(set(normal.keys()) | set(exchange.keys()))
    if not temps:
        raise SystemExit("Error: no data parsed from sample bin files.")

    lines = []
    lines.append("=" * 80)
    lines.append("Efficiency Comparison: Normal MC vs Exchange MC")
    lines.append("=" * 80)
    lines.append("")
    lines.append("Metrics: mean +/- standard_error (stderr from Sample-bin variance)")
    lines.append("")

    for t in temps:
        lines.append(f"--- T = {t:.4f} ---")
        for obs, label in [("E", "E_per_site"), ("C", "C_per_site"), ("M2", "M2")]:
            v_n = normal.get(t, {}).get(obs, [])
            v_e = exchange.get(t, {}).get(obs, [])
            m_n, se_n = mean_and_stderr(v_n)
            m_e, se_e = mean_and_stderr(v_e)
            lines.append(f"  {label}:")
            lines.append(f"    Normal:   {m_n:.6f} +/- {se_n:.6f}  (n={n_n})")
            lines.append(f"    Exchange: {m_e:.6f} +/- {se_e:.6f}  (n={n_e})")
            if se_n > 0:
                lines.append(f"    Ratio (Exchange/Normal stderr): {se_e / se_n:.3f}")
            lines.append("")
        lines.append("")

    lines.append("=" * 80)
    OUTPUT_FILE.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUTPUT_FILE}")

    if HAVE_MPL:
        fig, axes = plt.subplots(3, 1, figsize=(6.2, 9), sharex=True)
        for ax, (obs, label) in zip(axes, [("E", "E_per_site"), ("C", "C_per_site"), ("M2", "M2")]):
            m_n, se_n, m_e, se_e = [], [], [], []
            for t in temps:
                a, b = mean_and_stderr(normal.get(t, {}).get(obs, []))
                c, d = mean_and_stderr(exchange.get(t, {}).get(obs, []))
                m_n.append(a)
                se_n.append(b)
                m_e.append(c)
                se_e.append(d)
            ax.errorbar(temps, m_n, yerr=se_n, fmt="o-", capsize=3, label=f"Normal (n={n_n})")
            ax.errorbar(temps, m_e, yerr=se_e, fmt="s-", capsize=3, label=f"Exchange (n={n_e})")
            ax.set_ylabel(label)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[-1].set_xlabel("Temperature T")
        fig.suptitle("Normal vs Exchange MC (Sample-bin stderr)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(FIG_OUTPUT, dpi=150)
        plt.close(fig)
        print(f"Wrote {FIG_OUTPUT}")


if __name__ == "__main__":
    main()
