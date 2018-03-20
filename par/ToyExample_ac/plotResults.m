clc
close all
clear all
addpath('../model')

geometry.NX=100;  % Number of grid points in X
geometry.NY=100;  % Number of grid points in Y
geometry.NZ=1;  % Number of grid points in Z
geometry.DH=50;   % Spatial grid sampling
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice

colorbarRange.min=2500;
colorbarRange.max=4500;

load 'Input/sources.mtx';
acquisition.sources=spconvert(sources(2:end,:));
load 'Input/receiver.mtx';
acquisition.receiver=spconvert(receiver(2:end,:));

plotModel ('vp',colorbarRange,15,geometry,acquisition,'ToyModel','ToyModel.out','ToyModel_true')