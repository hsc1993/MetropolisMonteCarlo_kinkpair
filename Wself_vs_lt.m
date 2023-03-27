clear;close all;clc


a = [0.2,0.5,0.8,1.2,1.6,2]
lt = [0.0074,0.0049,0.0036,0.0026,0.0021,0.0017]
wself = [0.242,0.039,0.015,0.0067,0.0037,0.0024]


hold on

plot(a,lt,'LineWidth',2)
plot(a,wself,'LineWidth',2)
scatter(a,lt)
scatter(a,wself)
hold off
xlabel('core width a [b]')
ylabel('self energy [eV/A]')
LG = legend('line tension','non-singular anisotropic calculation');
LG.FontSize = 14;