#!/usr/bin/env python3
"""
Plot overlap with initial (step 0) at same temperature: with EXMC vs without EXMC.
Overlap = (1/N) sum_i s_i(step=0) . s_i(step=t) at each T.
Reads overlap_exmc.dat and overlap_no_exmc.dat from the given directory.
Usage:
  python3 plot_overlap.py [directory]
  (default directory: samples/square_L16_Ising)
"""

import sys
import os

def load_overlap(path):
    """Load overlap.dat (columns: int_T, T, overlap); return (T_list, overlap_list)."""
    T_list = []
    overlap = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) >= 3:
                T_list.append(float(cols[1]))
                overlap.append(float(cols[2]))
    return T_list, overlap

def main():
    if len(sys.argv) >= 2:
        datadir = sys.argv[1]
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        datadir = os.path.join(os.path.dirname(script_dir), 'samples', 'square_L16_Ising')

    path_exmc = os.path.join(datadir, 'overlap_exmc.dat')
    path_no_exmc = os.path.join(datadir, 'overlap_no_exmc.dat')

    if not os.path.isfile(path_exmc):
        print(f"Missing {path_exmc}", file=sys.stderr)
        print("Run run_overlap_compare.sh first.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(path_no_exmc):
        print(f"Missing {path_no_exmc}", file=sys.stderr)
        sys.exit(1)

    T_ex, ov_ex = load_overlap(path_exmc)
    T_no, ov_no = load_overlap(path_no_exmc)

    # Always write TSV for external plotting (first, so we have data even if plot fails)
    tsvpath = os.path.join(datadir, 'overlap_compare.tsv')
    with open(tsvpath, 'w') as f:
        f.write("T\toverlap_with_EXMC\toverlap_without_EXMC\n")
        for i in range(len(T_ex)):
            f.write("{:.6f}\t{:.10f}\t{:.10f}\n".format(T_ex[i], ov_ex[i], ov_no[i]))
    print("Saved " + tsvpath)

    # Optional: matplotlib plot (skip if import fails or crashes)
    do_plot = os.environ.get('PLOT_OVERLAP_USE_MATPLOTLIB', '0') == '1'
    if do_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot(T_ex, ov_ex, 'o-', label='with EXMC', color='C0', markersize=8)
            ax.plot(T_no, ov_no, 's--', label='without EXMC', color='C1', markersize=8)
            ax.set_xlabel('T')
            ax.set_ylabel('Overlap with step 0 (same T)')
            ax.set_title('Overlap (1/N) sum_i s_i(0).s_i(t) at same temperature')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color='gray', linestyle=':', alpha=0.7)
            fig.tight_layout()
            outpath = os.path.join(datadir, 'overlap_compare.png')
            fig.savefig(outpath, dpi=150)
            plt.close()
            print("Saved " + outpath)
        except Exception as e:
            print("Plot failed: " + str(e), file=sys.stderr)
    else:
        print("Summary: T | overlap(with EXMC) | overlap(without EXMC)")
        for i in range(len(T_ex)):
            print("  {:.3f}  | {:+.4f} | {:+.4f}".format(T_ex[i], ov_ex[i], ov_no[i]))
        print("To plot: use overlap_compare.tsv in Excel/gnuplot, or run with PLOT_OVERLAP_USE_MATPLOTLIB=1")

if __name__ == '__main__':
    main()
