close all;clear;clc

reactionCoordinate = [0,0.1041,0.2032,0.3021 ,0.4012,0.5017,0.6012,0.7008,0.8007,0.9017,1]; %reaction coordinate
Up_mev_b = [0,3.1537,9.8913,19.0456,24.9915,26.6873,25.0082,18.9957,9.8166,3.054,0]; %meV/b


burgers = 2.4829;
latticeConstant = burgers*2/sqrt(3);
h0 = latticeConstant*sqrt(6)/3;

x = reactionCoordinate*h0;%A
Up = Up_mev_b/1000/burgers; %eV/A

%% fitting
Up_fun=@(U0,xdata)U0/2*(1-cos(2*pi*xdata/h0));
U0_init = 1;
U0 = lsqcurvefit(Up_fun,U0_init,x,Up)



%% plotting
hold on
x_plot = 0:0.01*h0:h0;
scatter(x,Up)
plot(x_plot,Up_fun(U0,x_plot))
hold off
title('fitting for iron potential')
legend(strcat('U0 = ',num2str(U0),'eV/A'))
xlabel('x [A]')
ylabel('Up [eV/A]')


