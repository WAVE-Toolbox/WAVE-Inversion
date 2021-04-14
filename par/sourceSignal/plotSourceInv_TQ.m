clear all;close all;
addpath('../configuration');
addpath('../common');

dhType = '';
NoiseType = '';
NoisedB = '';
blockType = '';  

modelName = 'Jiangwan_Wall3';
observationType = 'Crosshole';
equationType = 'EMEM';
dimension = '2D';
modelType = '_Inv';
HPCType = '_HPC';
config=conf(['configuration_' modelName '_' ...
    observationType '_' equationType dimension NoiseType modelType HPCType '.txt']);
modelType = '_True';
configTrue=conf(['configuration_' modelName '_' ...
    observationType '_' equationType dimension modelType '.txt']);

% Output file
% switch for saving snapshots to picture file 1=yes (jpg) 2= yes (png) other=no
imagesave=0;
writefiles=1;

stage=5;
showLegend = 0;
showResidual = 0;
showmax=30;
legendType = 2; % 1 = parameter symblicName; 2 = inversionTypeName
% 3 = parameter symblicName + inversionTypeName
showTitle=1; showXlabel=1; showYlabel=1;

normalize=config.getValue('normalizeTraces');
inversionType = [-1 config.getValue('inversionType')];
parameterisation = [0 config.getValue('parameterisation')];
exchangeStrategy = [0 config.getValue('exchangeStrategy')];
orientation = 'vertical';
wiggleType='vararea';
writeSource=configTrue.getValue('writeSource');
writeSource=0;
sources = readSourcesfromConfig(configTrue);
DT=config.getValue('seismoDT');
fileFormat=config.getValue('SeismogramFormat');
if DT > 1e-8 
    labelTime='Time (ms)';
else
    labelTime='Time (ns)';
end
[Nshot n]=size(sources);
seismogram=[];seismogramTrue=[];
T0damp=10e-9; T1damp=30e-9; damp_factor=1e-2/DT;
for ishot=1:Nshot
    SOURCE_TYPE=sources(ishot,5);% Source Type (1=P,2=vX,3=vY,4=vZ)
    component = getSeismogramComponent(equationType,SOURCE_TYPE);
    
    %% Read seismogram
    sourceSeismogramFilename=config.getString('sourceSeismogramFilename');
    filenameInv=['../' sourceSeismogramFilename '.stage_' num2str(stage)...
        '.shot_' num2str(ishot-1) '.' component];
    sourceInv=readSeismogram(filenameInv,fileFormat);
    T=1*DT:DT:size(sourceInv,2)*DT;
    time_damping=ones(size(sourceInv));
    time_damping(T<T0damp)=time_damping(T<T0damp).*exp(damp_factor*(T(T<T0damp)-T0damp));
    time_damping(T>T1damp)=time_damping(T>T1damp).*exp(-damp_factor*(T(T>T1damp)-T1damp));
    sourceInv=sourceInv.*time_damping;
    writeSourceFilename=configTrue.getString('writeSourceFilename');
    if writeSource~=0        
        filenameTrue=['../' writeSourceFilename '_shot_' num2str(ishot-1)];
        sourceTrue=readSeismogram(filenameTrue,fileFormat);
        seismogramTrue=[seismogramTrue;sourceTrue];
    end
    seismogram=[seismogram;sourceInv];
end
%% write
filenameInv=['../' sourceSeismogramFilename '.stage_' num2str(stage) '.' component];
filenameTrue=['../' writeSourceFilename '.' component];
filenameStart=['../' writeSourceFilename '.' component];
if writefiles~=0
    answer=overwritedlg(filenameInv,fileFormat);
    if answer==0
        return;
    else
        filenameInv=filenameInv
    end
    writeSeismogram(filenameInv,seismogram,fileFormat);
    if writeSource~=0  
        answer=overwritedlg(filenameTrue,fileFormat);
        if answer==0
            return;
        else
            filenameTrue=filenameTrue
        end
        writeSeismogram(filenameTrue,seismogramTrue,fileFormat);
    end
end
%% plot
imageScale=1;
symblicName='';
legendNameAll = getLegendName(legendType,equationType,...
    inversionType,parameterisation,exchangeStrategy,symblicName);
lineSettingAll = getLinesettingInv(inversionType,parameterisation,exchangeStrategy); 
titleName = '';
[titleLabelSettingAll] = getTitleLabelSettingAll(showTitle,showXlabel,showYlabel...
    ,inversionType,parameterisation,exchangeStrategy,titleName,equationType,legendType);
offset.value=[1:Nshot];
offset.labelName='Source number';
syntheticDatanameAll = getFilenameModelAll(filenameTrue,filenameStart,{filenameInv}...
    ,inversionType);
residualScale=1;
[h answer] = plotSeismogram_multi(DT,showmax,...
    syntheticDatanameAll,imagesave,orientation,lineSettingAll,legendNameAll,offset,...
    titleLabelSettingAll,showLegend,showResidual,normalize,wiggleType...
    ,imageScale,fileFormat,residualScale);
% ylim([0 max(T)*2/3]);

addpath('../configuration');
addpath('../common');