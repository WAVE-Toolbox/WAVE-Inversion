clearvars; close all;

Min=2500;
Max=4500;

parameter='vp';
iteration=2;

load ../acquisition/sources.mtx;
sources=spconvert(sources(2:end,:));
load '../acquisition/receiver.mtx';
receiver=spconvert(receiver(2:end,:));



%% Define input parameter

NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

filename=['model.It' num2str(iteration) '.' parameter '.mtx']; % File name of the model
%% Read model
model=readModelfromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
subplot(1,3,3)
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
title([parameter ' model at iteration ' num2str(iteration)])
% 
%% Define input parameter
filename=['model.' parameter '.mtx']; % File name of the model


%% Read model
model=readModelfromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);


%% Plot

subplot(1,3,2)
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,1)*DH,sources(:,2)*DH,'*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'.')
title([parameter ' starting model'])

%% Define input parameter
filename=['model_true.' parameter '.mtx']; % File name of the model

%% Read model
model=readModelfromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
subplot(1,3,1)
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,1)*DH,sources(:,2)*DH,'*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'.')
title([parameter ' true model'])