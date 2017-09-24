#!/usr/bin/gnuplot

#set style line 1  lw 0.1;

plot 'velocity_particle_1.txt' u 1:($2) every 1 with l ls 1

pause 1
reread
