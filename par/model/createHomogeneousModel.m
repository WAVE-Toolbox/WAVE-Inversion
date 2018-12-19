clearvars; %close all;

%% Input parameter
velocityP=3500; % P-wave velocity in m/s
velocityS=2000; % S-wave velocity in m/s
density=2000; % Density in Kg/m^3
pertub=0;

NX=100;
NY=100;

N=100*100; % Number of grid points
filename='model'; % Base filename

%% Calculation

% Set homogeneous vectors
VP=velocityP*ones(1,N);
VS=velocityS*ones(1,N);
RHO=density*ones(1,N);

if pertub==1;
index=1;
for jj=1:NY
for ii=1:NX
    
    if ii>30 && ii<70 && jj>40 && jj<60;
        VP(1,index)=4500;
        VS(1,index)=2500;
        RHO(1,index)=2500;
    end
    
    index=index+1;
end
end
end

% write to file
writeVector2mtx([filename '.vp.mtx'],VP);
writeVector2mtx([filename '.vs.mtx'],VS);
writeVector2mtx([filename '.density.mtx'],RHO);