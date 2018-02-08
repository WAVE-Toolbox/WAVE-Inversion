clearvars; close all;

%% Define input parameter
filename='grad_vp.It2.mtx'; % File name of the model
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model
gradient=readGradientsFromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
imagesc(X,Y,gradient(:,:,LAYER)/max(max(abs(gradient))))
colorbar
caxis([-1 1])
xlabel('X in meter')
ylabel('Y in meter')
title('gradient')