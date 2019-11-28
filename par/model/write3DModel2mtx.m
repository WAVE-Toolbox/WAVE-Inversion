function [model1D]=write3DModel2mtx(filename,model)
% write3DModel2mtx  writes 3D array models to 1D mtx files and returns
% model vector

%   model1D = write3DModel2mtx(filename,model)
%   filename: filename of the model
%   model: 3D array
%   input must be model(Y,X,Z)
%   major order of output is X,Z,Y
%   

model1D=reshape(permute(model,[2 3 1]),[],1);
writeVector2mtx(filename,model1D);

end