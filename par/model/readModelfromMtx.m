function [model]=readModelfromMtx(filename,NX,NY,NZ)
% readModelfromMtx  read 1D models into 3D array

%   model = readModelfromMtx(filename,NX,NY,NZ)
%   filename: filename of the model
%   NX,NY,NZ number of gridpoints in each direction+
%
%   major order of input should by X,Z,Y
%   output is model(Y,X,Z)



modelVec=readVectorfromMtx(filename);
model=permute((reshape(modelVec(:),[NX, NZ, NY])),[3 1 2]);

end