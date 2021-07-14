clear all;close all;
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

% Output file
% switch for saving snapshots to picture file 1=yes (jpg) 2= yes (png) other=no
imagesave=0;
writefiles=1;
stage=1;
useDamping=0; T0damp=30e-9; T1damp=100e-9; 
showLegend = 0;
showResidual = 0;
showmax=30;
legendType = 2; % 1 = parameter symblicName; 2 = inversionTypeName
% 3 = parameter symblicName + inversionTypeName
showTitle=0; showXlabel=1; showYlabel=1;

normalize=1;
inversionType = [0 1];
parameterisation = [0 0];
exchangeStrategy = [0 0];
orientation = 'vertical';
wiggleType='vararea';
writeSource=0;
source = readSourcesfromConfig(configTrue);
DT=config.getValue('seismoDT');
fileFormat=config.getValue('SeismogramFormat');
if DT > 1e-8 
    labelTime='Time (ms)';
else
    labelTime='Time (ns)';
end
[Nshot n]=size(source);
seismogram=[];seismogramTrue=[];
dampFactor=1e-2/DT;
for ishot=1:Nshot
    shotnr=source(ishot,1);
    SOURCE_TYPE=source(ishot,5);% Source Type (1=P,2=vX,3=vY,4=vZ)
    component = getSeismogramComponent(equationType,SOURCE_TYPE);
    
    %% Read seismogram
    sourceSeismogramFilename=config.getString('sourceSeismogramFilename');
    filenameInv=['../' sourceSeismogramFilename '.stage_' num2str(stage)...
        '.shot_' num2str(shotnr) '.' component];
    sourceInv=readSeismogram(filenameInv,fileFormat);
    T=1*DT:DT:size(sourceInv,2)*DT;
    timeDamping=ones(size(sourceInv));
    if useDamping==1
        timeDamping(T<T0damp)=timeDamping(T<T0damp).*exp(dampFactor*(T(T<T0damp)-T0damp));
        timeDamping(T>T1damp)=timeDamping(T>T1damp).*exp(-dampFactor*(T(T>T1damp)-T1damp));
    end
    sourceInv=sourceInv.*timeDamping;
    writeSourceFilename=configTrue.getString('writeSourceFilename');
    if writeSource~=0        
        filenameTrue=['../' writeSourceFilename '_shot_' num2str(shotnr)];
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
    end
    writeSeismogram(filenameInv,seismogram,fileFormat);
    if writeSource~=0  
        answer=overwritedlg(filenameTrue,fileFormat);
        if answer==0
            return;
        end
        writeSeismogram(filenameTrue,seismogramTrue,fileFormat);
    end
end
%% plot
imageScale=1.2;
symblicName='';
legendNameAll = getLegendName(legendType,equationType,...
    inversionType,parameterisation,exchangeStrategy,symblicName);
lineSettingAll = getLinesettingInv(equationType,inversionType,parameterisation,exchangeStrategy); 
lineSettingAll{2}.Color='k';
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

addpath('../configuration');
addpath('../common');