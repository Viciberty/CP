set terminal jpeg
set isosamples 50,50 
show palette colornames


set output 'nature_100.jpg' 
set xlabel 'r/fm'
set ylabel 'z/fm'
set zlabel 'Potential(C/fm)'
set xyplane 0.1
splot "nature_100" with lines


set output 'MC_1000.jpg' 
set xlabel 'r/fm'
set ylabel 'z/fm'
set zlabel 'Potential(C/fm)'
set xyplane 0.1
splot "MC_1000" with lines

set output 'Quadratic_1000.jpg' 
set xlabel 'r/fm'
set ylabel 'z/fm'
set zlabel 'Potential(C/fm)'
set xyplane 0.1
splot "Quadratic_1000" with lines

set output 'nature5000.jpg' 
set xlabel 'r/fm'
set ylabel 'z/fm'
set zlabel 'Potential(C/fm)'
set xyplane 0.1
splot "nature5000" with lines










