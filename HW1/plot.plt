set terminal latex
set output "output.tex"
set xlabel "$Deformation$"
set ylabel "$Energy$"
set title "V(q)"
plot "V_q" using 1:2 title "2-order Lagrange Interpolation" 



set xlabel "$Deformation$"
set ylabel "$Energy^{-1}$"
set title "M(q)"
plot "M_q" using 1:2 title "2-order Lagrange Interpolation"




