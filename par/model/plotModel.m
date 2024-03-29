%% Function Name: plotModel
%
% Inputs:
%   parameter(string): model parameter eg. vp,vs,density
%
%   colorbarRange(struct): struct with entries min and max for a fixed colorbar
%
%   iteration(int): Nr. of the iterarion to be displayed
%
%   geometry(struct): struct with entries: NX,NY,NZ,DH,LAYER (layer of 3D model to
%   display as 2D slice)
%
%   acquisition(struct): struct with entries sources and receivers matrices
%
%   inversionModel(string): location of the model at itaration 'itaration'
%
%   startingModel(string): location  of the starting model 
%
%   true model(string): location of the "true" model in synthetic tests
%
% Modes:
%   inversion model, starting model + inversion model, true model +
%   starting model + inversion model
%
% Outputs:
%   None

%
% $Date: March 20, 2018
% ________________________________________
function plotModel(parameter,colorbarRange,stage,iteration,geometry,acquisition,inversionModel,startingModel,trueModel)


Min=colorbarRange.min;
Max=colorbarRange.max;

sources=acquisition.sources;
receiver=acquisition.receiver;

%% Define input parameter

NX=geometry.NX;  % Number of grid points in X
NY=geometry.NY;  % Number of grid points in Y
NZ=geometry.NZ;  % Number of grid points in Z
DH=geometry.DH;   % Spatial grid sampling
LAYER=geometry.LAYER+1; % Define layer of 3D model to display as 2D slice


%% Read model
if iteration > 0
model=readModelfromMtx(inversionModel,NX,NY,NZ);
end
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
switch nargin
case 9
subplot(3,1,3)
case 8
subplot(1,2,2)
case default
 subplot(1,1,1)     
end

if iteration > 0
imagesc(X,Y,model(:,:,LAYER))
end
caxis([Min Max])
colorbar('location','eastoutside');
colormap('jet')
xlabel('X in meter')
ylabel('Y in meter')
title([parameter ' model at iteration ' num2str(iteration)])
hold on
% plot(sources(:,2)*DH,sources(:,3)*DH,'r*')
% plot(receiver(:,1)*DH,receiver(:,2)*DH,'kv')
axis equal
axis([min(X) max(X) min(Y) max(Y)]);

if nargin > 7
%% Read model
model=readModelfromMtx(startingModel,NX,NY,NZ);

%% Plot
if nargin > 7
subplot(3,1,2)
else
subplot(1,2,1)
end
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','eastoutside');
xlabel('X in meter')
ylabel('Y in meter')
title([parameter ' starting model'])
hold on
plot(sources(:,2)*DH,sources(:,3)*DH,'r*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'kv')
axis equal
axis([min(X) max(X) min(Y) max(Y)]);
end

if nargin > 8
%% Read model
model=readModelfromMtx(trueModel,NX,NY,NZ);

%% Plot
subplot(3,1,1)
imagesc(X,Y,model(:,:,LAYER))
caxis([Min Max])
colorbar('location','eastoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,2)*DH,sources(:,3)*DH,'r*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'kv')
axis equal
axis([min(X) max(X) min(Y) max(Y)]);
title([parameter ' true model'])
numbOFcont=8;
subplot(3,1,3)
hold on
contour(X,Y,model(:,:,LAYER),numbOFcont,'k-','LineWidth',1);
hold off
end
