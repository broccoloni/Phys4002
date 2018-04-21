# gnuplot script file for plotting heatmap of data

#sets up the plot
set pm3d map
set xlabel "r[km]" font ',20'
set ylabel "z[km]" font ',20'
set xtics font ',18'
set ytics font ',18'
set key right bottom

# set the area you want the graph to show
set yrange [0:40]
set xrange [18.5:150]

#define new colour system
set palette defined (0 "blue", 1 "light-green", 2 "yellow", 3 "orange", 4 "red", 5 "magenta")

#creates multiple plots at once (rows, columns) if uncommented
#set multiplot layout 2,1

#changes terminal size if uncommented because sometimes the axes titles don't display because gnuplot is old
#set terminal wxt size 1680,1050 
#set terminal wxt size 1000,750

#set the title to the plot
set title "Electron Antineutrino Opacity [unitless]" font ',20'

#Uncomment the lines to choose which variable you want to plot -> all of these plot with the neutrino surface 
#Note: cblabel is colour bar label -> doesn't always show up because gnuplot is old

#Temperature -> this plot also has vertical lines at r = 12.7 and 18.5 to show Risco and inner radius of accretion disk
#set cblabel "Temperature [MeV]"  font ',20'
#splot 'simul_all.dat' u 1:2:3 with points ps 2 pt 5 palette, 'Riscoline.dat' with points pt 7 ps 0.5 lt -1 title '', 'v_surface2.dat' with points pt 7 ps 0.5 lt -1 title '','vA_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Density
#set cblabel "Density [g/cm^3]"  font ',20'
#splot 'simul_all.dat' u 1:2:4 with points ps 2 pt 5 palette title 'density', 'v_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Electron fraction
#set cblabel "Electron Fraction"  font ',20'
#splot 'simul_all.dat' u 1:2:5 with points ps 2 pt 5 palette title 'electron fraction', 'v_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Neutrino mfp
#set cblabel "Neutrino MFP [km]"  font ',20'
#splot 'simul_all.dat' u 1:2:6 with points ps 2 pt 5 palette title 'Neutrino MFP [km]', 'v_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Antineutrino mfp
#set cblabel "Antineutrino MFP [km]"  font ',20'
#splot 'simul_all.dat' u 1:2:7 with points ps 2 pt 5 palette title 'Antineutrino MFP [km]', 'vA_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#plot of opacities
# datafile in order radius, height, opac_v, opac_Va

#Opacity for antineutrinos
#set cblabel "Opacity neutrinos" font ',20'
#splot 'opac_all.dat' u 1:2:3 with points ps 2 pt 5 palette title 'Electron Neutrino Opacity', 'v_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Opacity for neutrinos
#set cblabel "Opacity antineutrinos"  font ',20'
splot 'opac_all.dat' u 1:2:4 with points ps 2 pt 5 palette title 'antineutrino opacity', 'vA_surface2.dat' with points pt 7 ps 0.5 lt -1 title ''

#Replot the plots and save as a png because the graphs often look better ... again, gnuplot is old software
set terminal png enhanced truecolor font times 14 size 700, 500
set output 'M3a08t20Aopac.png'
replot

#uncomment this if multiplot is used
#unset multiplot




