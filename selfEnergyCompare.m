clear; clc;

%% Wself
%% importing data
Wself_file = fopen('Wself.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
Wself_data = fscanf(Wself_file,formatSpec,sizeA);
fclose(Wself_file);

%% Wself0
%% importing data
Wself0 = fopen('Wself0.txt','r');
formatSpec = '%f %f';
sizeA = [4 Inf];
Wself0_data = fscanf(Wself0,formatSpec,sizeA);
fclose(Wself0);

%% plotting data

hold on
%plot(Wself0_data(1,:),Wself0_data(4,:),'LineWidth',2)
hold off


%% line tension
%% importing data
lt_file = fopen('lt.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
lt_data = fscanf(lt_file,formatSpec,sizeA);
fclose(lt_file);

%% core 
%% importing data
core_file = fopen('core.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
core_data = fscanf(core_file,formatSpec,sizeA);
fclose(core_file);



%% plotting data
hold on
scatter(Wself_data(1,:),Wself_data(2,:),'d')
plot(lt_data(1,:),lt_data(2,:),'LineWidth',2)
plot(core_data(1,:),core_data(2,:),'LineWidth',4)

legend('Wself','lt','core')
hold off
title('Wself comparison','FontSize',20)
xlabel('theta [degree]')
ylabel('eV/A')



