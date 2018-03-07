clearvars; close all;

%% Define input parameter
filename='grad.It_1.vp.mtx'; % File name of the gradient
shot = 1;
% filenameHessian=['Hessian.shot_' num2str(shot) '.mtx']; % File name of the gradient
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D gradient to display as 2D slice

%% Read gradient
gradient=readGradientFromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot gradient
figure
imagesc(X,Y,gradient(:,:,LAYER)/max(max(abs(gradient))))
colorbar
% caxis([-1 1])
xlabel('X in meter')
ylabel('Y in meter')
title('gradient')

%% Read Hessian
% hessian=readGradientFromMtx(filenameHessian,NX,NY,NZ);

%% Plot Hessian
% figure
% imagesc(X,Y,hessian(:,:,LAYER))
% colorbar
% xlabel('X in meter')
% ylabel('Y in meter')
% title('Hessian')

