%% Function Name: plotModel_el
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
function plotModel_el(parameter,colorbarRange,stage,iteration,geometry,acquisition,configuration_inv,profile_layer,inversionModel,startingModel,trueModel)


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

X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Read model
if iteration > 0
    filename=[inversionModel '.stage_' num2str(stage) '.It_' num2str(iteration) '.' parameter '.mtx']; % File name of the model
    if str2double(configuration_inv.getString("useVariableGrid"))
        [X,Y,model]=readmeshmodel(filename,configuration_inv.getString("coordinateFilename"),DH,LAYER,true);
    else
        model=readModelfromMtx(filename,NX,NY,NZ);

    end
    model_inv=model;
end

%% Plot
figure
switch nargin
case 11
subplot(1,3,3)
case 10
subplot(1,2,2)
case default
 subplot(1,1,1)     
end

if iteration > 0
imagesc(X,Y,model(:,:,LAYER+1))
end
caxis([Min Max])
colorbar('location','southoutside');
colormap('jet')
xlabel('X in meter')
ylabel('Y in meter')
title([parameter ' model at iteration ' num2str(iteration)])


if nargin > 9
%% Read model
filename=[startingModel '.' parameter '.mtx']; % File name of the model
model=readModelfromMtx(filename,NX,NY,NZ);
model_start=model;

%% Plot
if nargin > 10
subplot(1,3,2)
else
subplot(1,2,1)
end
imagesc(X,Y,model(:,:,LAYER+1))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,2)*DH,sources(:,3)*DH,'*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'.')
title([parameter ' starting model'])
end

if nargin > 10
%% Read model
filename=[trueModel '.' parameter '.mtx']; % File name of the model
model=readModelfromMtx(filename,NX,NY,NZ);
model_true=model;

%% Plot
subplot(1,3,1)
imagesc(X,Y,model(:,:,LAYER+1))
caxis([Min Max])
colorbar('location','southoutside');
xlabel('X in meter')
ylabel('Y in meter')
hold on
plot(sources(:,2)*DH,sources(:,3)*DH,'*')
plot(receiver(:,1)*DH,receiver(:,2)*DH,'.')
title([parameter ' true model'])

if iteration > 0
    figure
    plot(model_start(:,profile_layer,LAYER+1))
    hold on
    plot(model_true(:,profile_layer,LAYER+1))
    plot(model_inv(:,profile_layer,LAYER+1))
    legend('starting model','true model','model after inversion')
    xlabel('Y in meter')
    ylabel('velocity in m/s')
    title(['Velocity profile at x = ',num2str(profile_layer*DH)])
end
end
