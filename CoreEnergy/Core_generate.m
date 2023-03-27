close all;clear;clc

b = 2.4829;%Angstrom  for iron

%screw core energy: 0.22 eV/Angstrom
%edge core energy: 0.80 eV/Angstrom

data = [0.000        0.22
        90.000       0.8
];
for i = 1:2
    angle(i) = data(i,1)/180*pi;
    ecore(i) = data(i,2);
end


fun_p = @(c,angle_data)(c(1)+c(2)*cos(2*angle_data));
c0 = [0.5,0.5];
c = lsqcurvefit(fun_p,c0,angle,ecore)

%{
fun_p = @(c,angle_data)(c(1)+c(2)*sin(2*angle_data)+c(3)*cos(2*angle_data)+...
    c(4)*sin(4*angle_data)+c(5)*cos(4*angle_data)+(-c(2)/3-c(4)*2/3)*sin(6*angle_data)+c(6)*cos(6*angle_data));
c0 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5];
c = lsqcurvefit(fun_p,c0,angle,ecore);
%}

%% Ecore expression

fun_ecore = @(c,d,a,angle_data)(c(1)+c(2)*sin(2*angle_data)+c(3)*cos(2*angle_data)+c(4)*sin(4*angle_data)...
    +c(5)*cos(4*angle_data)+(-c(2)/3-c(4)*2/3)*sin(6*angle_data)+c(6)*cos(6*angle_data))+...
    +(d(1)+d(2)*sin(2*angle_data)+d(3)*cos(2*angle_data)+...
    d(4)*sin(4*angle_data)+d(5)*cos(4*angle_data)+(-d(2)/3-d(4)*2/3)*sin(6*angle_data)+d(6)*cos(6*angle_data))*log(a/b);
%%
angle_fourier = 0:angle(end)/1000:angle(end)*2;

ECORE = fun_p(c,angle_fourier);

%ECORE(:,j) = fun_ecore(c_final,d_final,a(j),angle_fourier);

hold on
plot(angle_fourier,fun_p(c,angle_fourier),'LineWidth',3)
scatter(data(:,1)/180*pi,data(:,2))
legend('fitting','data points')
hold off
xlabel('angle [radian]')
ylabel('energy [eV/A]')
title('iron Ecore fitting')
