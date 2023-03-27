clc;clear;close all;

A = [-1/sqrt(6) 2/sqrt(6) -1/sqrt(6); -1/sqrt(2) 0 1/sqrt(2); 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];

x_coor = [1 0 0];
y_coor = [0 1 0];
z_coor = [0 0 1];

x_equi = 0;
h_saddle = 10;
N = 100;

for i = 1:20%fraction of a1 in AD
    x(i) = x_equi;
end
for i = 21:40%fraction of a1 in AD
    x(i) = x_equi+h_saddle/(N*0.2)*i-h_saddle;
end
for i = 41:60%fraction of w in AD
    x(i) = x_equi+h_saddle;
end
for i = 61:80%fraction of a2 in AD
    x(i) = x_equi-h_saddle/(N*0.2)*i+h_saddle*(1)/0.25;
end
for i = 81:100
    x(i) = x_equi;
end


for i = 1:100
    y(i) = 0;
    z(i) = i*2.7;
end

for i = 1:99
    dL_old(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2+(z(i+1)-z(i))^2);
end
sum_dL_old = sum(dL_old)


for i = 1:100
    trapezoid_old(i,:) = [x(i),y(i),z(i)];
    x_old(i) = trapezoid_old(i,1);
    y_old(i) = trapezoid_old(i,2);
    z_old(i) = trapezoid_old(i,3);
end
x_axis_old = [1 0 0]*50;
y_axis_old = [0 1 0]*50;
z_axis_old = [0 0 1]*50;
for i = 1:100
   hold on  
   scatter3(x_old(i),y_old(i),z_old(i),'b','filled')
   %{
   plot3([0,x_axis_old(1)],[0,x_axis_old(2)],[0,x_axis_old(3)],'b')
   plot3([0,y_axis_old(1)],[0,y_axis_old(2)],[0,y_axis_old(3)],'b')
   plot3([0,z_axis_old(1)],[0,z_axis_old(2)],[0,z_axis_old(3)],'b')
   %}
   hold off

end





%%  new coordinate
x_coor_trans = [-1 2 -1]/sqrt(6);
y_coor_trans = [-1 0 1]/sqrt(2);
z_coor_trans = [1 1 1]/sqrt(3);



%% transformation matrix
% A' = T*A T(1,1) = cos(i'*i)
T(3,3) = zeros;

T(1,1) = x_coor_trans*x_coor';
T(1,2) = x_coor_trans*y_coor';
T(1,3) = x_coor_trans*z_coor';

T(2,1) = y_coor_trans*x_coor';
T(2,2) = y_coor_trans*y_coor';
T(2,3) = y_coor_trans*z_coor';

T(3,1) = z_coor_trans*x_coor';
T(3,2) = z_coor_trans*y_coor';
T(3,3) = z_coor_trans*z_coor';

Tinv = inv(T);
for i = 1:100
    trapezoid_new(i,:) = trapezoid_old(i,:)*T;
    x_new(i) = trapezoid_new(i,1);
    y_new(i) = trapezoid_new(i,2);
    z_new(i) = trapezoid_new(i,3);
end

for i = 1:99
    dL_new(i) = sqrt((x_new(i+1)-x_new(i))^2+(y_new(i+1)-y_new(i))^2+(z_new(i+1)-z_new(i))^2);
end
sum_dL_new = sum(dL_new)
x_axis_new = [-1 2 -1]*50;
y_axis_new = [-1 0 1]*50;
z_axis_new = [1 1 1]*50;

for i = 1:100
   hold on  
   scatter3(x_new(i),y_new(i),z_new(i),'k','filled')
   %{
   plot3([0,x_axis_new(1)],[0,x_axis_new(2)],[0,x_axis_new(3)],'k')
   plot3([0,y_axis_new(1)],[0,y_axis_new(2)],[0,y_axis_new(3)],'k')
   plot3([0,z_axis_new(1)],[0,z_axis_new(2)],[0,z_axis_new(3)],'k')
   %}
   %scatter3(x_axis_new(1),x_axis_new(2),x_axis_new(3),8)
   hold off

end
xlabel('[100]','FontSize',20)
ylabel('[010]','FontSize',20)
zlabel('[001]','FontSize',20)
%Alan W: x,y,z [-1 2 -1], [-1 0 1], [111]

grid on
%{
set(gca,'xtick',[0:20:100])
set(gca,'ytick',[0:20:100])
set(gca,'ztick',[0:20:100])
%}




