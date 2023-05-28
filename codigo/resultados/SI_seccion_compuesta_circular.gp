
reset
#set dgrid3d 30,30
set hidden3d

unset key

set xlabel 'Mnx(kg-cm)' font "Times,14" offset -7,0
set ylabel 'Mny(kg-cm)' font "Times,14" offset 7,0
set zlabel 'Pn(kg)' font "Times,14" rotate by 90 offset -4,-1
set xtics font "Times,14" offset -3
set ytics font "Times,14" offset 3
set ztics font "Times,14"
set xtics 2000000
set ytics 2000000
set ztics 200000
unset ztics

set view 0,45

#set contour base
#set cntrparam levels 30

#unset surface

#set size square # Esto para 2D
set view equal xy # Esto para 3D

#set yrange[-10:25]
#set zrange[-10:1000]
splot "SI_seccion_compuesta_circular.txt" with points lt rgb("#009999") ps 0.25, "sol_seccion_compuesta_circular.txt" with lines linetype 1 linecolor rgb("#FF0000") linewidth 12.5

set term pdf dashed
set output "SI_seccion_compuesta_circular.pdf"
replot

set term pop

#set parametric
#replot 0.33333, u, v
#replot u,v,0
