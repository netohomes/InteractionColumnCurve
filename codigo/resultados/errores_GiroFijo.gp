reset
set grid
set yrange [-10:10]

set term wxt 1
plot "errores_ex_GiroFijo.txt" u 1:2 with lines

set term wxt 2
plot "errores_ey_GiroFijo.txt" u 1:2 with lines

set term wxt 3
plot "errores_e_GiroFijo.txt" u 1:2 with lines
