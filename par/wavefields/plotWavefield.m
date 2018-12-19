clearvars; close all;

%% Define input parameter
filename='test.mtx'; % File name of the model
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
NT=1000;
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model

fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
A=fscanf(fileID,'%e',[1 size(1)*size(2)]);
model=reshape(A(:),[size(2), size(1)]);
model=permute(reshape(A,[NX,NY,NZ,NT]),[2 1 3 4]);
%model=permute((reshape(model(:),[NX, NY, NZ])),[2 1 3]);
%model=readWavefieldfromMtx(filename,NX,NY,NZ,NT);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);


caxis_value=5.0e-8;
load 'seismic.map'
%% Plot

figure
colormap(seismic);

xlabel('X in meter')
ylabel('Y in meter')
for ii=1:4:NT;
imagesc(X,Y,model(:,:,LAYER,ii))
title( num2str(ii))
colorbar
caxis([-caxis_value caxis_value])
pause(0.1)

end