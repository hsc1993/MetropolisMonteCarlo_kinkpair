clear; clc;close all


U0 = 0.0112;  
N_land = 5000;

h0 = 2.3409;
blength = 2.4829;
currentTau = 0.00125;%eV/A
x_equi = 0.0777179;
h_predict = 1.69294;
x_predict = 1.77066;

for i = 1:N_land
    x_landscape(i) = h0/N_land*i;
    U_land(i) = U0/2*(1-cos(2*pi*x_landscape(i)/h0));
    U_land_gradiant(i) = pi*U0/h0*sin(2*pi*x_landscape(i)/h0);
    UW_land(i) = U_land(i)-currentTau*blength*x_landscape(i);
    
    
end


subplot(2,2,1)
hold on
plot(x_landscape,U_land)
plot([x_equi,x_equi],[0,max(U_land)])
plot([x_predict,x_predict],[0,max(U_land)])

hold off

subplot(2,2,2)
hold on
tau_gradiant = currentTau*blength;
plot(x_landscape,U_land_gradiant)
plot([x_equi,x_equi],[0,max(U_land_gradiant)])
plot([x_predict,x_predict],[0,max(U_land_gradiant)])
plot([0,1],[tau_gradiant,tau_gradiant])

hold off

subplot(2,2,3)
hold on
plot(x_landscape,UW_land)
plot([x_equi,x_equi],[0,max(UW_land)])
plot([x_predict,x_predict],[0,max(UW_land)])

hold off

subplot(2,2,4)
hold on
plot(x_landscape,U_land)
plot([x_equi,x_equi],[0,max(U_land)])
plot([x_predict,x_predict],[0,max(U_land)])

hold off






