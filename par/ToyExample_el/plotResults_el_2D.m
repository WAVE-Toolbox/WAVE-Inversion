clc
close all
clear all

addpath('../configuration')
config=conf('../ToyExample_el/Input/configuration.txt');
configTrue=conf('../ToyExample_el/Input/configuration_true.txt');

% general parameter ----------------------------------------------

stage=1;
iteration=3; % iteration number

% model parameter-------------------------------------------------
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice
parameter='vs';   % model parameter

colorbarRange.min=1800; %lower clip of the colorbar 
colorbarRange.max=2800; %upper clip of the colorbar

% seismogram parameter--------------------------------------------
skipTraces=3;  % every 'skipTraces' trace will be displayed
shot=0;        % shot number
component='vy'; % seismogram component


%% Usually, there is no need to change anything below this line

addpath('../model')
addpath('../seismograms')
addpath('../gradients')

geometry.NX=config.getValue('NX');  % Number of grid points in X
geometry.NY=config.getValue('NY');  % Number of grid points in Y
geometry.NZ=config.getValue('NZ');  % Number of grid points in Z
geometry.DH=config.getValue('DH');   % Spatial grid sampling


%% plot model

inversionModel=config.getString('ModelFilename');
startingModel=[config.getString('ModelFilename') '.out'];
trueModel=configTrue.getString('ModelFilename');

load ([config.getString('SourceFilename') '.mtx']);
acquisition.sources=spconvert(sources(2:end,:));
load ([config.getString('ReceiverFilename') '.mtx']);
acquisition.receiver=spconvert(receiver(2:end,:));

plotModel (parameter,colorbarRange,stage,iteration,geometry,acquisition,...
    inversionModel,startingModel,trueModel)

%% plot seismograms


DT=config.getValue('DT');
fieldData=config.getString('fieldSeisName');
syntheticData=config.getString('SeismogramFilename');

plotSeismogram(DT,stage,iteration,shot,component,skipTraces,syntheticData,fieldData);

%% plot gradient

gradientName=config.getString('gradientFilename');

if iteration > 0
plotGradient(parameter,stage,iteration,geometry,gradientName);
end

rmpath('../model')
rmpath('../seismograms')
rmpath('../gradients')
rmpath('../configuration')