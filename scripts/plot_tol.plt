set term postscript eps enhanced color
set output "tol.eps"
set logscale y
plot 'tol.out' u 1:2
