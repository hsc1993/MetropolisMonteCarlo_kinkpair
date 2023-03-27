close all; clear; clc;

%% Wself
%% importing data
Wself_an_file = fopen('Wself_an.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
Wself_an = fscanf(Wself_an_file,formatSpec,sizeA);
fclose(Wself_an_file);

%% plotting data
subplot(1,2,1)
hold on
%plot(Wself_an(1,:),Wself_an(2,:),'LineWidth',2)
scatter(Wself_an(1,1:end-1),Wself_an(2,1:end-1),'LineWidth',2)

hold off
legend('Wself an [eV/A]')
title('Wself an','FontSize',20)
xlabel('theta [degree]','FontSize',20)
ylabel('Wself an [eV/A]','FontSize',20)


%% Wself
%% importing data
Wself_ratio_file = fopen('Wself_ratio.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
Wself_ratio = fscanf(Wself_ratio_file,formatSpec,sizeA);
fclose(Wself_ratio_file);

%% plotting data
%{
subplot(1,2,1)
hold on
%plot(Wself_ratio(1,:),Wself_ratio(2,:),'LineWidth',2)
%scatter(Wself_ratio(1,:),Wself_ratio(2,:),'LineWidth',2)

hold off
legend('Wself an/Wself iso')
title('Wself_an/Wself_iso','FontSize',20)
xlabel('theta [degree]','FontSize',20)
ylabel('Wself an/Wself iso [ratio]','FontSize',20)
%}


%% r1 coordinate
%% importing data
r1coor = fopen('r1coor.txt','r');
formatSpec = '%f %f';
sizeA = [3 Inf];
r1coor_data = fscanf(r1coor,formatSpec,sizeA);
fclose(r1coor);



%% plotting data
subplot(1,2,2)
hold on
line_111 = 10000*2.7*[0 1;0 1;0 1];
plot3(line_111(1,:),line_111(2,:),line_111(3,:))
plot3(r1coor_data(1,1:end/2),r1coor_data(2,1:end/2),r1coor_data(3,1:end/2),'LineWidth',2)
%plot3(r1coor_data(1,end/2+1:end),r1coor_data(2,end/2+1:end),r1coor_data(3,end/2+1:end),'LineWidth',2)
scatter3(r1coor_data(1,end/2+1:end),r1coor_data(2,end/2+1:end),r1coor_data(3,end/2+1:end),'LineWidth',2)

hold off
title('dislocation loop','FontSize',20)

