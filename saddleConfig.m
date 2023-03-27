clear; clc;close all

%% initial config
initial_config_file = fopen('initial_config.txt','r');
formatSpec = '%f %f';
sizeA = [3 Inf];
initial_config_data = fscanf(initial_config_file,formatSpec,sizeA);
fclose(initial_config_file);

%% plotting data
subplot(2,2,1)
hold on
initial_config_data_shrink = initial_config_data/1;
%scatter3(initial_config_data_shrink(1,:),initial_config_data_shrink(2,:),initial_config_data_shrink(3,:),3,'b')
hold off
title('initial config')

%% initial config ref
%% importing data
initial_config_ref_file = fopen('initial_config_ref.txt','r');
formatSpec = '%f %f';
sizeA = [3 Inf];
initial_config_data_ref = fscanf(initial_config_ref_file,formatSpec,sizeA);
fclose(initial_config_ref_file);
%% plotting data
subplot(2,2,2)
hold on
%scatter(initial_config_data_ref(3,:),initial_config_data_ref(1,:),3,'b')
hold off
title('initial config')


%% saddle config ref
%% importing data
saddle_config_file = fopen('saddle_config.txt','r');
formatSpec = '%f %f';
sizeA = [3 Inf];
saddle_config_data = fscanf(saddle_config_file,formatSpec,sizeA);
fclose(saddle_config_file);
%% plotting data
subplot(2,2,3)
hold on
saddle_config_data_shrink = saddle_config_data/1;
pureScrew = [[0; 0; 0],[saddle_config_data_shrink(1,end); saddle_config_data_shrink(2,end); saddle_config_data_shrink(3,end)]];
%plot3(pureScrew(1,:),pureScrew(2,:),pureScrew(3,:),'LineWidth',2)
%scatter3(saddle_config_data_shrink(1,:),saddle_config_data_shrink(2,:),saddle_config_data_shrink(3,:),6,'k')
hold off
title('saddle config')


%% importing data
saddle_config_ref_file = fopen('saddle_config_ref.txt','r');
formatSpec = '%f %f';
sizeA = [3 Inf];
saddle_config_ref_data = fscanf(saddle_config_ref_file,formatSpec,sizeA);
fclose(saddle_config_ref_file);

%% importing data
stressDependents_file = fopen('stressDependents.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
stressDependents_data = fscanf(stressDependents_file,formatSpec,sizeA);
fclose(stressDependents_file);
x_equi = stressDependents_data(1)
x_predict = stressDependents_data(2)
%% plotting data
figure
hold on


plot([0,250],[2.3409,2.3409]) % h0
plot([0,250],[x_predict,x_predict]) % x_predict
plot([0,250],[x_equi,x_equi]) % x_equi


scatter(saddle_config_ref_data(3,:),saddle_config_ref_data(1,:),6,'k')
plot(initial_config_data_ref(3,:),initial_config_data_ref(1,:),'b')


hold off
title('saddle config')
LG = legend('h0 ','x predict ','x equi ');
%LG = legend(strcat('a = ',num2str(0.3),'b', '  enthalpy = ',num2str(0.541695),'eV'));
LG.FontSize = 14;
grid on;


