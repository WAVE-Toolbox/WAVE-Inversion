% plotGradientTests
clc
close all
clear all

iteration=1;

parameter='vp';   % model parameter

%% Usually, there is no need to change anything below this line

addpath('../configuration')
config=conf('../ci/configuration_ci.2D.acoustic.txt');

geometry.NX=config.getValue('NX');  % Number of grid points in X
geometry.NY=config.getValue('NY');  % Number of grid points in Y
geometry.NZ=config.getValue('NZ');  % Number of grid points in Z
geometry.DH=config.getValue('DH');   % Spatial grid sampling
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice

gradientName=config.getString('gradientFilename');

plotGradient(parameter,iteration,geometry,gradientName);