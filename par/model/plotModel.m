function plotModel(parameter,colorbarRange,iteration,geometry,acquisition,inversionModel,startingModel,trueModel)

nargin

Min=colorbarRange.min;
Max=colorbarRange.max;

sources=acquisition.sources;
receiver=acquisition.receiver;

%% Define input parameter

NX=geometry.NX;  % Number of grid points in X
NY=geometry.NY;  % Number of grid points in Y
NZ=geometry.NZ;  % Number of grid points in Z
DH=geometry.DH;   % Spatial grid sampling
LAYER=geometry.LAYER; % Define layer of 3D model to display as 2D slice


%% Read model
filename=[inversionModel '.It_' num2str(iteration) '.' parameter '.mtx']; % File name of the model
model=readModelfromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
switch nargin
case 8
subplot(1,3,3)
case 7
subplot(1,2,2)
case default
 subplot(1,1,1)     
end
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
title([parameter ' model at iteration ' num2str(iteration)])



if nargin > 6
%% Read model
filename=[startingModel '.' parameter '.mtx']; % File name of the model
model=readModelfromMtx(filename,NX,NY,NZ);

%% Plot
if nargin > 7
subplot(1,3,2)
else
subplot(1,2,1)
end
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,1)*DH,sources(:,2)*DH,'*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'.')
title([parameter ' starting model'])
end

if nargin > 7
%% Read model
filename=[trueModel '.' parameter '.mtx']; % File name of the model
model=readModelfromMtx(filename,NX,NY,NZ);

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
end
