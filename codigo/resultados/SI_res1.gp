
reset
#set dgrid3d 30,30
set hidden3d

unset key

set xlabel 'Mnx(kg-cm)' font "Times,14" offset -5,0
set ylabel 'Mny(kg-cm)' font "Times,14" offset 5,0
set zlabel 'Pn(kg)' font "Times,14" rotate by 90 offset -1,-1
set xtics font "Times,14"
set ytics font "Times,14"
set ztics font "Times,14"
set xtics 2000000
set ytics 2000000
set ztics 200000

set view 60,45

#set contour base
#set cntrparam levels 30

#unset surface

#set yrange[-10:25]
#set zrange[-10:1000]
#splot "SI_res1.txt" with points lt rgb("#009999") ps 0.25

#set dgrid3d 10,10
#set style data lines
splot "SI_res1.txt" with pm3d

set term pdf dashed
set output "SI_res1.pdf"
replot

set term pop

#set parametric
#replot 0.33333, u, v
#replot u,v,0
