clc
close all
clear all

addpath('../configuration','../model')
config=conf('Input/configuration.txt');
configTrue=conf('Input/configuration_true.txt');

% general parameter ----------------------------------------------

stage=3;
iteration=1; % iteration number

% model parameter-------------------------------------------------
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice
parameter='vp';   % model parameter

colorbarRange.min=1900; %lower clip of the colorbar 
colorbarRange.max=4000; %upper clip of the colorbar

% seismogram parameter--------------------------------------------
skipTraces=3;  % every 'skipTraces' trace will be displayed
shot=8;        % shot number
component='p'; % seismogram component

filename_source_true='sourceInversion/source_signal_ricker_shot_8.mtx';
filename_source_guess='sourceInversion/source_signal_sin3_shot_8.mtx';
filename_source_inv='sourceInversion/invSource.stage_4.shot_8.p.mtx';

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


delimiterIn = ' ';
headerlinesIn = 1;

filename_sources=[config.getString('SourceFilename') '.txt'];

A = importdata(filename_sources,delimiterIn,headerlinesIn);
acquisition.sources=A.data;

filename_receiver=[config.getString('ReceiverFilename') '.txt'];
B = importdata(filename_receiver,delimiterIn,headerlinesIn);
acquisition.receiver=B.data;

%%
plotModel (parameter,colorbarRange,stage,iteration,geometry,acquisition,...
    inversionModel,startingModel,trueModel)

%% plot seismograms


DT=config.getValue('DT');
T=config.getValue('T');
fieldData=config.getString('fieldSeisName');
syntheticData=config.getString('SeismogramFilename');

plotSeismogram(DT,stage,iteration,shot,component,skipTraces,syntheticData,fieldData);

%% plot gradient

gradientName=config.getString('gradientFilename');

if iteration > 0
plotGradient(parameter,stage,iteration,geometry,gradientName);
end

%% plot source
source_true=readVectorfromMtx2(filename_source_true);
source_guess=readVectorfromMtx2(filename_source_guess);
source_inv=readVectorfromMtx2(filename_source_inv);
t=DT:DT:T;
figure
title('Source signal')
plot(t,source_true)
hold on
plot(t,source_guess)
plot(t,source_inv)
hold off
legend('True source signal','Guessed source signal','Inverted source signal')
xlabel('time in s')
ylabel('amplitude')

rmpath('../model')
rmpath('../seismograms')
rmpath('../gradients')
rmpath('../configuration')
