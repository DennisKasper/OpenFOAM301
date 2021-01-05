set output "unitCubeStepSizeStudy.svg"

set terminal svg enhanced fsize 8  size 262,226
#set terminal postscript solid color eps enhanced  "Palatino" size 7.4cm,6cm
#set terminal epslatex dashed size 16cm,10cm color colortext

reset 
set border lw 0.5 

set xlabel 'Step size $ \\Delta \\boldsymbol\{x\} $'
set xtics font "Times-Roman, 8" 
#set xrange [0:800]
set logscale x 10

set grid back lw 0.05 ls -1

set ytics font "Times-Roman, 8" 
set ylabel 'Sensitivity $ \\text\{d\} L / \\text\{d\} \\boldsymbol\{x\} $ '
#set logscale y 10
#set yrange [-0.006:0.006]
#set lmargin 10

set format y '%3.1E'
set format x '$10\^\{%T\}$'
#set format y '$10\^\{%T\}$'
#set format x "%3.0f"
set key right center box  lc -1 lt 1 lw 0.2 Left

# Point 116
set label '-0.40' at 1.2E-9,-2.17E-05
set label '-0.12' at 1.2E-8,-2.17E-05
set label '-0.17' at 1.2E-7,-2.17E-05
set label '-0.17' at 1.2E-6,-2.17E-05
set label '-0.17' at 1.2E-5,-2.17E-05
set label '-0.17' at 1.2E-4,-2.17E-05
set label '-0.17' at 1.2E-3,-2.17E-05
set label '-0.07' at 1.2E-2,-2.17E-05
set label '-9.18' at 1.2E-1,-2.17E-05

# Point 117
set label '0.03' at 1.2E-9,1.3E-05
set label '0.02' at 1.2E-8,1.3E-05
set label '0.12' at 1.2E-7,1.3E-05
set label '0.11' at 1.2E-6,1.3E-05
set label '0.11' at 1.2E-5,1.3E-05
set label '0.11' at 1.2E-4,1.3E-05
set label '0.11' at 1.2E-3,1.3E-05
set label '0.10' at 1.2E-2,1.3E-05
set label '0.68' at 1.2E-1,1.3E-05

set pointsize 1.2
     
plot 'sensitivityDataAdj.dat' using 1:5 every ::62::70   with points  title '   Adj116', \
     'sensitivityDataFD.dat'  using 1:5 every ::62::70      with points  title '   FD116', \
     'sensitivityDataAdj.dat' using 1:5 every ::72::80   with points  title '   Adj117', \
     'sensitivityDataFD.dat'  using 1:5 every ::72::80     with points  title '   FD117'
     
set output



# every command (plotting specific rows)
#every A:B:C:D:E:F
#A: line increment
#B: data block increment
#C: The first line
#D: The first data block
#E: The last line
#F: The last data block