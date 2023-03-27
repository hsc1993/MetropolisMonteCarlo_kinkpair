% elastic constants rotation

%C11=169.9; C12=122.6; C44=76.2; %Cu Mishin in GPa
%C11=243.4; C12=145.0; C44=116.0; % Fe Mendelev in GPa
%C11=516.0; C12=198.0; C44=161.0; % W MEAM in GPa
C11=5.160; C12=1.980; C44=1.610; % W MEAM in GPa


C = cubic_elast_stiff(C11,C12,C44);

disp(sprintf('C(1,1,1,1) = %g',C(1,1,1,1)));
disp(sprintf('C(2,2,2,2) = %g',C(2,2,2,2)));
disp(sprintf('C(3,3,3,3) = %g',C(3,3,3,3)));

disp(sprintf('C(1,1,2,2) = %g',C(1,1,2,2)));
disp(sprintf('C(1,1,3,3) = %g',C(1,1,3,3)));
disp(sprintf('C(2,2,3,3) = %g',C(2,2,3,3)));

disp(sprintf('C(2,3,2,3) = %g',C(2,3,2,3)));
disp(sprintf('C(3,1,3,1) = %g',C(3,1,3,1)));
disp(sprintf('C(1,2,1,2) = %g',C(1,3,1,3)));

e1=[1 0 0]; e2=[0 1 0]; e3=[0 0 1];

e2p=[-1 0 1]; e3p=[1 1 1]; e1p=[-1 2 -1];
e1p=e1p/norm(e1p);e2p=e2p/norm(e2p);e3p=e3p/norm(e3p);

Q = [dot(e1p,e1) dot(e1p,e2) dot(e1p,e3)
    dot(e2p,e1) dot(e2p,e2) dot(e2p,e3)
    dot(e3p,e1) dot(e3p,e2) dot(e3p,e3)];

Cp = rotate_tensor_4th (C, Q);

disp(sprintf('Cp(1,1,1,1) = %g',Cp(1,1,1,1)));
disp(sprintf('Cp(2,2,2,2) = %g',Cp(2,2,2,2)));
disp(sprintf('Cp(3,3,3,3) = %g',Cp(3,3,3,3)));

disp(sprintf('Cp(1,1,2,2) = %g',Cp(1,1,2,2)));
disp(sprintf('Cp(1,1,3,3) = %g',Cp(1,1,3,3)));
disp(sprintf('Cp(2,2,3,3) = %g',Cp(2,2,3,3)));

disp(sprintf('Cp(2,3,2,3) = %g',Cp(2,3,2,3)));
disp(sprintf('Cp(3,1,3,1) = %g',Cp(3,1,3,1)));
disp(sprintf('Cp(1,2,1,2) = %g',Cp(1,3,1,3)));



%% ALAN added 

fileID = fopen('Cijkl_ext.out','w');
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
               B = [i; j; k; l; Cp(i,j,k,l)];
               fprintf(fileID,'%d %2d %2d %2d %12.8e\n',B);
            end
        end
    end
end
fclose(fileID);



burgers = 2.7223;
fileID = fopen('Burgers.out','w');
fprintf(fileID,'%.8f',burgers);

fclose(fileID);








