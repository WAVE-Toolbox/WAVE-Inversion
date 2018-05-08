function plotGradient(parameter,stage, iteration,geometry,gradientName)

load 'seismic.map'


%% Define input parameter
filename=[gradientName '.stage_' num2str(stage) '.It_' num2str(iteration)  '.' parameter '.mtx']; % File name of the gradient
NX=geometry.NX;  % Number of grid points in X
NY=geometry.NY;  % Number of grid points in Y
NZ=geometry.NZ;  % Number of grid points in Z
DH=geometry.DH;   % Spatial grid sampling
LAYER=geometry.LAYER; % Define layer of 3D gradient to display as 2D slice

%% Read gradient
gradient=readGradientFromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot gradient
figure
colormap(seismic);
imagesc(X,Y,gradient(:,:,LAYER)/max(max(abs(gradient))))
colorbar
 caxis([-1 1])
xlabel('X in meter')
ylabel('Y in meter')
title('normed gradient')


