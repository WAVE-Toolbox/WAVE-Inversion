clc
close all
clear all

iteration=2;

parameter='vp';   % model parameter

colorbarRange.min=2500; %lower clip of the colorbar 
colorbarRange.max=4500; %upper clip of the colorbar

%% Usually, there is no need to change anything below this line

currentDir=pwd;
cd ../ % change to par directory
addpath('configuration')
addpath(currentDir)
config=conf('./ci/configuration_ci.2D.acoustic.txt'); % filename of configuration

geometry.NX=config.getValue('NX');  % Number of grid points in X
geometry.NY=config.getValue('NY');  % Number of grid points in Y
geometry.NZ=config.getValue('NZ');  % Number of grid points in Z
geometry.DH=config.getValue('DH');   % Spatial grid sampling
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice


sources = load (config.getString('SourceFilename'));
acquisition.sources=spconvert(sources(2:end,:));
receiver = load (config.getString('ReceiverFilename'));
acquisition.receiver=spconvert(receiver(2:end,:));

inversionModel=config.getString('ModelFilename');
startingModel=[config.getString('ModelFilename') '.out'];

plotModel (parameter,colorbarRange,iteration,geometry,acquisition,inversionModel,startingModel)

cd (currentDir)
rmpath(currentDir)
rmpath('configuration')