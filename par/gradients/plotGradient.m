function plotGradient(parameter,stage, iteration,geometry,gradientName, configuration_inv)

load 'seismic.map'


%% Define input parameter
filename=[gradientName '.stage_' num2str(stage) '.It_' num2str(iteration) '.' parameter '.mtx']; % File name of the gradient
NX=geometry.NX;  % Number of grid points in X
NY=geometry.NY;  % Number of grid points in Y
NZ=geometry.NZ;  % Number of grid points in Z
DH=geometry.DH;   % Spatial grid sampling
LAYER=geometry.LAYER+1; % Define layer of 3D gradient to display as 2D slice

%% Read gradient
if str2double(configuration_inv.getString("useVariableGrid"))
        [X,Y,gradient]=readmeshmodel(filename,configuration_inv.getString("coordinateFilename"),DH,LAYER,1);
else
    gradient=readGradientFromMtx(filename,NX,NY,NZ);
    X=0:DH:(NX*DH-DH);
    Y=0:DH:(NY*DH-DH);
end


%% Plot gradient
figure
colormap(seismic);
imagesc(X,Y,gradient(:,:,LAYER))
colorbar
caxis([-max(max(max(abs(gradient(:,:,LAYER))))) max(max(max(abs(gradient(:,:,LAYER)))))])
xlabel('X in meter')
ylabel('Y in meter')
title('gradient')
axis equal
axis([min(X) max(X) min(Y) max(Y)]);
