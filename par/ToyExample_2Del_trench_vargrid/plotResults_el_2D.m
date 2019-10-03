clc
close all
clear all

addpath('../configuration')
config=conf('../ToyExample_2Del_trench_vargrid/Input/configuration.txt');
configTrue=conf('../ToyExample_2Del_trench_vargrid/Input/configuration_true.txt');

% general parameter ----------------------------------------------

stage=4;
iteration=4; % iteration number

% model parameter-------------------------------------------------
geometry.LAYER=0; % Define layer of 3D model to display as 2D slice
profile_layer=140; % Display velocity profile at x layer
parameter='vs';   % model parameter

colorbarRange.min=100; %lower clip of the colorbar
colorbarRange.max=320; %upper clip of the colorbar

% seismogram parameter--------------------------------------------
skipTraces=3;  % every 'skipTraces' trace will be displayed
shot=1;        % shot number
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
startingModel=config.getString('ModelFilename');
trueModel=configTrue.getString('ModelFilename');

delimiterIn = ' ';
headerlinesIn = 1;

filename_sources=[config.getString('SourceFilename') '.txt'];

A = importdata(filename_sources,delimiterIn,headerlinesIn);
acquisition.sources=A.data;

filename_receiver=[config.getString('ReceiverFilename') '.txt'];
B = importdata(filename_receiver,delimiterIn,headerlinesIn);
acquisition.receiver=B.data;

plotModel_el(parameter,colorbarRange,stage,iteration,geometry,acquisition,...
    config,profile_layer,inversionModel,startingModel,trueModel)

%% plot seismograms
DT=config.getValue('DT');
fieldData=[config.getString('fieldSeisName'),'.stage_',num2str(stage)];
syntheticData=config.getString('SeismogramFilename');

plotSeismogram(DT,stage,iteration,shot,component,skipTraces,syntheticData,fieldData);

%% plot gradient

gradientName=config.getString('gradientFilename');

if iteration > 0
plotGradient(parameter,stage,iteration,geometry,gradientName,config);
end

rmpath('../model')
rmpath('../seismograms')
rmpath('../gradients')
rmpath('../configuration')