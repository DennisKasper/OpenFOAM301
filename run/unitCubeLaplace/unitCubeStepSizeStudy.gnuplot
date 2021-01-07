set output "unitCubeStepSizeStudy.svg"

set terminal svg enhanced fsize 8  size 262,226
#set terminal postscript solid color eps enhanced  "Palatino" size 7.4cm,6cm
#set terminal epslatex dashed size 16cm,10cm color colortext

reset 
set border lw 0.5 

set xlabel "Step size"
set xtics font "Times-Roman, 10" 
#set xrange [0:800]
set logscale x 10

set grid back lw 0.05 ls -1

set ytics font "Times-Roman, 10" 
set ylabel "Sensitivity"
#set logscale y 10
#set yrange [-3.6:-4.5]
#set lmargin 10

set format y "%3.1E"
set format x '$10\^\{%T\}$'
#set format y '$10\^\{%T\}$'
#set format x "%3.0f"
set key right top box  lc -1 lt 1 lw 0.3 Left
                     
plot    'sensitivityDataAdj.txt'        using 1:5       every ::3::14   with linespoints        ps 0.5  pt 9 lt 1 lw 1 lc rgb "orange"    title   '  AdjX', \
        'sensitivityDataFD.txt'         using 1:5       every ::3::14   with linespoints        ps 0.5  pt 5 lt 1 lw 1 lc rgb "red"     title   '  FDX', \
        'sensitivityDataAdj.txt'        using 1:6       every ::3::14   with linespoints        ps 0.5  pt 9 lt 1 lw 1 lc rgb "blue"    title   '  AdjY', \
        'sensitivityDataFD.txt'         using 1:6       every ::3::14   with linespoints        ps 0.5  pt 5 lt 1 lw 1 lc rgb "yellow"     title   '  FDY', \
        'sensitivityDataAdj.txt'        using 1:7       every ::3::14   with linespoints        ps 0.5  pt 9 lt 1 lw 1 lc rgb "green"    title   '  AdjZ', \
        'sensitivityDataFD.txt'         using 1:7       every ::3::14   with linespoints        ps 0.5  pt 5 lt 1 lw 1 lc rgb "red"     title   '  FDZ', \
                #'airFoil2DForces.txt' using 1:3         with linespoints                   ps 0.5  pt 7 lt 1 lw 1 lc rgb "black" title '  Drag', \
                #'airFoil2DModForces.txt' using 1:3       with linespoints                   ps 0.5  pt 13 lt 1 lw 1 lc rgb "green" title '  Dragmod'
set output




# every command (plotting specific rows)
#every A:B:C:D:E:F
#A: line increment
#B: data block increment
#C: The first line
#D: The first data block
#E: The last line
#F: The last data block