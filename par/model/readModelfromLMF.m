function [model]=readModelfromLMF(filename,NX,NY,NZ)
% readModelfromLMF  read 1D models into 3D array

%   model = readModelfromLMF(filename,NX,NY,NZ)
%   filename: filename of the model
%   NX,NY,NZ number of gridpoints in each direction+
%
%   major order of input should by X,Z,Y
%   output is model(Y,X,Z)

modelVec=readVectorfromLMF(filename);
% major order of vector is X, Z, Y
model=permute((reshape(modelVec(:),[NX, NZ, NY])),[3 1 2]);



end
