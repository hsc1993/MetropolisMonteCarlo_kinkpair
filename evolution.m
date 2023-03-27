close all;clear;clc

%%
tension_anisotropic_evolution_file = fopen('tension_anisotropic_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
tension_anisotropic_evolution_data = fscanf(tension_anisotropic_evolution_file,formatSpec,sizeA);
fclose(tension_anisotropic_evolution_file);

iterator = 1:numel(tension_anisotropic_evolution_data);

subplot(2,4,1)
plot(iterator,tension_anisotropic_evolution_data)
title('tension anisotropic')


%%
tension_evolution_file = fopen('tension_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
tension_evolution_data = fscanf(tension_evolution_file,formatSpec,sizeA);
fclose(tension_evolution_file);

iterator = 1:numel(tension_evolution_data);

subplot(2,4,2)
plot(iterator,tension_evolution_data)
title('tension')


%%
core_evolution_file = fopen('core_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
core_evolution_data = fscanf(core_evolution_file,formatSpec,sizeA);
fclose(core_evolution_file);


subplot(2,2,2)
plot(iterator,core_evolution_data)
title('core')

%%
U_evolution_file = fopen('U_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
U_evolution_data = fscanf(U_evolution_file,formatSpec,sizeA);
fclose(U_evolution_file);


subplot(2,4,5)
plot(iterator,U_evolution_data)
title('U')

%%
W_evolution_file = fopen('W_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
W_evolution_data = fscanf(W_evolution_file,formatSpec,sizeA);
fclose(W_evolution_file);

iterator = 1:numel(W_evolution_data);

subplot(2,4,6)
plot(iterator,W_evolution_data)
title('W')

%%
enthalpy_evolution_file = fopen('enthalpy_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
enthalpy_evolution_data = fscanf(enthalpy_evolution_file,formatSpec,sizeA);
fclose(enthalpy_evolution_file);

iterator = 1:numel(enthalpy_evolution_data);

subplot(2,4,7)
plot(iterator,enthalpy_evolution_data)
title('enthalpy')
legend(strcat('enthalpy = ',num2str(enthalpy_evolution_data(end))))

enthalpy_anisotropic_evolution_file = fopen('enthalpy_anisotropic_evolution.txt','r');
formatSpec = '%f %f';
sizeA = [1 Inf];
enthalpy_anisotropic_evolution_data = fscanf(enthalpy_anisotropic_evolution_file,formatSpec,sizeA);
fclose(enthalpy_anisotropic_evolution_file);

iterator = 1:numel(enthalpy_evolution_data);

subplot(2,4,8)
plot(iterator,enthalpy_anisotropic_evolution_data)
legend(strcat('enthalpy anisotropic = ',num2str(enthalpy_anisotropic_evolution_data(end))))
title('enthalpy anisotropic')






