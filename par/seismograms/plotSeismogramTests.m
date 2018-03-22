clc
close all
clear all

iteration=0;
skipTraces=3;  % every 'skipTraces' trace will be displayed
shot=5;        % shot number
component='p'; % seismogram component

%% Usually, there is no need to change anything below this line

currentDir=pwd;
addpath(currentDir)
cd ../ % change to par directory
addpath('./configuration')
config=conf('./ci/configuration_ci.2D.acoustic.txt');


DT=config.getValue('DT');
fieldData=config.getString('fieldSeisName');
syntheticData=config.getString('SeismogramFilename');


plotSeismogram(DT,iteration,shot,component,skipTraces,syntheticData,fieldData);

cd (currentDir)