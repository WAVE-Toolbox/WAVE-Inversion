%% Converts a raw binary float matrix model to a mtx vector
% This script can be used to convert model that are saved as
% raw floating number in binary format to the mtx format.

%% Clear memory
clearvars;
close all;

%% Input parameter
NX=100; % Number of grid points in X
NY=100; % Number of grid points in Y
NZ=100; % Number of grid points in Z (set to 1 for 2-D case)
filenameRawBin='model.bin'; % Filename input raw binary
filenameMtx='model.density.mtx'; % Filename output mtx

%% Convert
% Open raw binary file
[fid]=fopen(filenameRawBin,'r');
matrix=fread(fid,NX*NY*NZ,'float');
matrix=reshape(matrix,NY,NX,NZ);
matrix=permute(matrix,[2 1 3]);
% Write matrix to mtx format
writeVector2mtx(filenameMtx,matrix(:));
