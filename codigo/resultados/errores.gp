
reset
#set dgrid3d 30,30
set hidden3d

unset key

set xlabel 'theta' font "Times,14"
set ylabel 'desp' font "Times,14"
set zlabel 'error' font "Times,14"
set xtics font "Times,14"
set ytics font "Times,14"
set ztics font "Times,14"
set ztics 1

set view 35,45

set contour base
set cntrparam levels 30

#unset surface

set yrange[-10:25]
#set zrange[-10:1000]
splot "errores_e.txt" w l

set term pdf dashed
set output "errores_3D.pdf"
replot

set term pop

#set parametric
#replot 0.33333, u, v
#replot u,v,0
