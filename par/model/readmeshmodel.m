%% Function Name: readmeshvariablemodel
%
% Inputs:
%   filename: filename
%   
%
%
% Outputs:
%   regular grid model

%
% $Date: October 1, 2019
% ________________________________________
function [x_reg,y_reg,model_reg]=readmeshmodel(filename,coordinates_base,DH,layerZ,interpolate_bool)


%% Read model
% 
x=readVectorfromMtx([coordinates_base,'X.mtx']);
y=readVectorfromMtx([coordinates_base,'Y.mtx']);
z=readVectorfromMtx([coordinates_base,'Z.mtx']);
model=readVectorfromMtx(filename);

if length(x) ~= length(model)
    error('Error size of coordinate vector differs from size of model')
end

indices=find(z==layerZ);
model_layer=model(indices);
x_layer=x(indices);
y_layer=y(indices);

if ~interpolate_bool
    
    model_reg=zeros(max(y),max(x));
    for ii=1:length(y_layer)
       xi=x_layer(ii);
       yi=y_layer(ii);
    
       model_reg(yi+1,xi+1)=model_layer(ii);
    end
    
else

    [xq,yq] = meshgrid(0:1:max(x), 0:1:max(y));
    model_reg = griddata(x_layer,y_layer,model_layer,xq,yq);
    
end
 
x_reg=(0:max(x))*DH;
y_reg=(0:max(y))*DH;
