clear all
NROWS = 64

  
  for i = 1:NROWS
   X(i) = 1/(2*NROWS) + (i-1)*(1/NROWS);
  end
  for i = 1:NROWS
   Y(i) = 1/(2*NROWS) + (i-1)*(1/NROWS);
  end
  
  load Velocity_x.dat
  U = zeros(NROWS,NROWS);
  a=1;
  for i=1:NROWS
    for j =1:NROWS
     U(i,j) = Velocity_x(a);
     a = a+1;
   end
  end
  figure(1)
  contourf(X,Y,U,15)
  title('u velocity contour at steady state Re = 100')
  xlabel('Length along x')
  ylabel('Length along y')
  colorbar
 
  load Velocity_y.dat
  V = zeros(NROWS,NROWS);
  a=1;
  for i=1:NROWS
    for j =1:NROWS
     V(i,j) = Velocity_y(a);
     a = a+1;
   end
  end
  figure(2)
  contourf(X,Y,V,15)
  title('v velocity contour at steady state Re = 100')
  xlabel('Length along x')
  ylabel('Length along y')
  colorbar
  
  load Pressure.dat
  P = zeros(NROWS,NROWS);
  a=1;
  for i=1:NROWS
    for j =1:NROWS
     P(i,j) = Pressure(a);
     a = a+1;
   end
  end
  figure(3)
  contourf(X,Y,P)
  title('Pressure contour at steady state Re = 100')
  xlabel('Length along x')
  ylabel('Length along y')
  colorbar
  
  load AVG_L2_c.dat
  avg_L2_u = AVG_L2_c(:,1);
  %avg_L2_v = AVG_L2(:,2);
  time = AVG_L2_c(:,3);
  figure(4)
  semilogy(avg_L2_u,time)
  hold on
  load AVG_L2_b.dat
  avg_L2_u = AVG_L2_b(:,1);
  %avg_L2_v = AVG_L2(:,2);
  time = AVG_L2_b(:,3);
  semilogy(avg_L2_u,time)
  hold on
  load AVG_L2_a.dat
  avg_L2_u = AVG_L2_a(:,1);
  %avg_L2_v = AVG_L2(:,2);
  time = AVG_L2_a(:,3);
  semilogy(avg_L2_u,time)
  hold on
  %semilogy(avg_L2_v,time)
  xlabel('Time (s)')
  ylabel('L2 norm')
  title('Error History for 64x64')
  legend('dt=0.0005','dt=0.003','dt=0.001')
  grid on
  hold off

  load Mid_u32.dat
  figure(5)
  for i = 1:32
   Ya(i) = 1/(2*32) + (i-1)*(1/32);
  end
  plot(Mid_u32,Ya)  
  xlabel('u velocity in midplane')
  ylabel('vertical length of cavity (y)')
  hold on;
  grid on;
  for i = 1:64
   Yb(i) = 1/(2*64) + (i-1)*(1/64);
  end
  load Mid_u64.dat
  plot(Mid_u64,Yb)
  for i = 1:128
   Yc(i) = 1/(2*128) + (i-1)*(1/128);
  end
  load Mid_u128.dat
  plot (Mid_u128,Yc)
  u_midplane_ghia = [0, -.03717, -.04192, -.04775, -.06434, -.10150, -.15662,...
  -.21090 ,-.20581, -.13641, .00332, .23151, 0.68717, .73722, .78871, .84123, 1];
%u_midplane_1000 = [0 -.18109 -.20196 -.22220 -.29730 -.38289 -.27805...
%-.10648 -.06080 .05702 .18719 .33304 .46604 .51117 .57492 .65928 1];
  y_midplane_ghia = [0, .0547, 0.0625, .0703, .1016, .1719, .2813, .4531, .5,...
  .6172, .7344, .8516, .9531, .9609, .9688, .9766, 1];
  plot(u_midplane_ghia,y_midplane_ghia,'o')
  title('Mid plane x velocity')
  legend('32x32','64x64','128x128','Ghia')
  
  hold off
  load Mid_v32.dat
  figure(6)
  for i = 1:32
   Xa(i) = 1/(2*32) + (i-1)*(1/32);
  end
  plot(Xa, Mid_v32)
  hold on
  xlabel('horizontal length of cavity (x)')
  ylabel('v velocity in midplane')
  grid on;hold on;
  load Mid_v64.dat
  for i = 1:64
   Xb(i) = 1/(2*64) + (i-1)*(1/64);
  end
  plot(Xb,Mid_v64)
  load Mid_v128.dat
  for i = 1:128
   Xc(i) = 1/(2*128) + (i-1)*(1/128);
  end
  plot (Xc,Mid_v128)
  v_midplane_ghia = [0 0.09233 .10091 .10890 .12317 .16077 .17507 .17527...
  .05454 -.24533 -.22445 -.16914 -.10313 -.08864 -.07391 -.05906 0];
%v_midplane_1000 = [0 0.27485 .29012 .30353 .32627 .37095 .33075 .32235 ...
%.02526 -.31966 -.42665 -.51550 -.39188 -.33714 -.27669 -.21388 0];
  x_midplane_ghia = [0 .0625 0.0703 0.0781 0.0938 0.1563 .2266 .2344 .5...
  .8047 .8594 .9063 .9453 .9531 .9609 .9688 1];
  plot(x_midplane_ghia,v_midplane_ghia,'o')
  title('Mid plane y velocity')
  legend('32x32','64x64','128x128','Ghia')
    

  
