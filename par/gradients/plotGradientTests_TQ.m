close all; clear all;
addpath('../configuration');
addpath('../common');

modelName = 'Horstwalde';
observationType = 'Surface_Stream';
equationType = 'TMEM';
NoiseType = '';
modelType = 'Inv';
HPCType = 'HPC';
dimension=cellMerge({equationType,'2D'},0);
configFilename=cellMerge({'configuration',modelName,observationType,dimension...
    ,NoiseType,modelType,HPCType},1);
modelType = 'True';
configTrueFilename=cellMerge({'configuration',modelName,observationType,dimension...
    ,modelType},1);
configFilename=addfileSuffix(configFilename,5);
configTrueFilename=addfileSuffix(configTrueFilename,5);
config=conf(configFilename);
configTrue=conf(configTrueFilename);

imagesave = 0;
copy_inv = 0;
parameter='epsilonEMr';   % model parameter
stage=5;
iteration=1;
shotnr=0;
gradientType=1; % 1=gradient,2=crossGradient,3=crossGradientDerivative
gradientPerShot=0; % 1 = gradientPerShot, other = gradientSum;

DIR_PATH_NEW = 'data/';
invertParameterType = 'PorositySaturation50';
bandPass = 'BP520Hz15MHz';
timeGain = '';
sourceType = '';
depthGain = '';
NoisedB = '';
NoisedB = cellMerge({NoiseType,NoisedB},0);
Insert_name = cellMerge({invertParameterType,...
    bandPass,NoisedB,sourceType,depthGain,timeGain},1); 

%% Usually, there is no need to change anything below this line
geometry.LAYER=1; % Define layer of 3D model to display as 2D slice
fileFormat=config.getValue('FileFormat');
valueMax = 0e-1;
if valueMax == 0    
    clim=[ ];
else
    clim=[-valueMax valueMax];
end
plotGradient_TQ(clim,parameter,stage,shotnr,iteration,geometry...
    ,config,configTrue,imagesave,gradientPerShot,gradientType);
% cope the file to a defined directory
if copy_inv == 1
    % get the filename
    gradientName=['../' config.getString('gradientFilename')];
    if gradientPerShot == 1
        filename=[gradientName '.stage_' num2str(stage) '.It_' num2str(iteration)...
            '.Shot_' num2str(shotnr) '.' parameter]; % File name of the gradient
    else
        filename=[gradientName '.stage_' num2str(stage) '.It_' num2str(iteration)...
            '.' parameter]; % File name of the gradient
    end
    copyfile_image(filename,fileFormat,DIR_PATH_NEW,Insert_name,config);
end
rmpath('../common');
rmpath('../configuration');