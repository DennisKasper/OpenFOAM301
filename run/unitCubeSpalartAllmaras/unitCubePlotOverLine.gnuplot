set output "unitCubeSensFieldPlotOverLine.svg"

set terminal svg enhanced fsize 8  size 262,226
#set terminal postscript solid color eps enhanced  "Palatino" size 7.4cm,6cm
#set terminal epslatex dashed size 16cm,10cm color colortext

reset 
set datafile separator ","
set border lw 0.5 

set xlabel '$ x $ in \\si\{\\metre\}'
set xtics font "Times-Roman, 8" 
#set xrange [0:800]
#set logscale x 10

set grid back lw 0.05 ls -1

set ytics font "Times-Roman, 8" 
set ylabel 'Sensitivity $ \\text\{d\} L / \\text\{d\} \\boldsymbol\{x\} $ '
#set logscale y 10
#set yrange [-0.006:0.006]
#set lmargin 10

set format y '%3.1E'
set format x '%2.1f'
#set format y '$10\^\{%T\}$'
#set format x "%3.0f"
set key right top box  lc -1 lt 1 lw 0.2 Left

set pointsize 1.0
     
plot 'fieldData.csv' using 27:23 every ::1::101   with linespoints  title '   Adj', \
     'fieldData.csv' using 27:24 every ::1::101   with linespoints  title '   FD'
     
set output



# every command (plotting specific rows)
#every A:B:C:D:E:F
#A: line increment
#B: data block increment
#C: The first line
#D: The first data block
#E: The last line
#F: The last data block