# gnuplot script: plot overlap with step 0 at same T
# Run from samples/square_L16_Ising: gnuplot ../../scripts/plot_overlap.gp

set terminal pngcairo size 640,480 font "Sans,12"
set output "overlap_compare.png"
set xlabel "T"
set ylabel "Overlap with step 0 (same T)"
set title "Overlap (1/N) sum_i s_i(0).s_i(t) at same temperature"
set key top right
set grid
set datafile separator "\t"
plot "overlap_compare.tsv" using 1:2 with linespoints pt 7 title "with EXMC", \
     "overlap_compare.tsv" using 1:3 with linespoints pt 5 title "without EXMC"
