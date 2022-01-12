clear all;close all;
addpath('../configuration');
addpath('../common');

modelName = 'EttlingerCB';
observationType = 'Surface_Stream';
equationType = 'TMEM';
NoiseType = '';
modelType = 'Inv';
HPCType = 'HPC';
dimension=cellMerge({equationType,'2D'},0);
configFilename=cellMerge({'configuration',modelName,observationType,dimension...
    ,modelType,NoiseType,HPCType},1);
modelType = 'True';
configTrueFilename=cellMerge({'configuration',modelName,observationType,dimension...
    ,modelType},1);
configFilename=addfileSuffix(configFilename,5);
configTrueFilename=addfileSuffix(configTrueFilename,5);
config=conf(configFilename);
configTrue=conf(configTrueFilename);

copy_inv = 1; copy_true_start = 0;
imagesave = 0;
DIR_PATH_NEW = 'data/';
invertParameterType = 'EpsilonEMSigmaEM_Field';
bandPass = 'BP550MHz';
NoisedB = '';
NoisedB = cellMerge({NoiseType,NoisedB},0);
sourceType = '';
timeGain = '';
depthGain = '';
Insert_name = cellMerge({invertParameterType,...
    bandPass,NoisedB,sourceType,depthGain,timeGain},1); 

writefiles=1;
stage=5;
useDamp=0; T0damp=15e-9; T1damp=25e-9; 
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
wiggleType='wiggle';
writeSourceTrue=1;
shotIncr=config.getAndCatch('ShotIncr',0);
source = readSourcesfromConfig(config,shotIncr);
DT=config.getValue('seismoDT');
fileFormat=config.getValue('SeismogramFormat');
if DT > 1e-8 
    labelTime='Time (ms)';
else
    labelTime='Time (ns)';
end
[Nshot n]=size(source);
seismogramInv=[];seismogramTrue=[];
dampFactor=2e-2/DT;
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
    if useDamp==1
        timeDamping(T<T0damp)=timeDamping(T<T0damp).*exp(dampFactor*(T(T<T0damp)-T0damp));
        timeDamping(T>T1damp)=timeDamping(T>T1damp).*exp(-dampFactor*(T(T>T1damp)-T1damp));
    end
    sourceInv=sourceInv.*timeDamping;
    writeSourceFilename=configTrue.getString('writeSourceFilename');
    if writeSourceTrue~=0        
        filenameTrue=['../' writeSourceFilename '.shot_' num2str(shotnr)];
        sourceTrue=readSeismogram(filenameTrue,fileFormat);
        seismogramTrue=[seismogramTrue;sourceTrue];
    end
    seismogramInv=[seismogramInv;sourceInv];
end
filenameInv=['../' sourceSeismogramFilename '.stage_' num2str(stage) '.' component];
filenameTrue=['../' writeSourceFilename '.' component];

%% plot
imageScale=1.2;
symblicName='';
showMarker=0;
showText=[1,1]*0; % 0=no, 1st:1=white, 2=black; 2nd:1=left bottom, 2=left top, 3=right top, 4=center
textName={'',''};
legendNameAll = getLegendName(legendType,equationType,...
    inversionType,parameterisation,exchangeStrategy,symblicName);
lineSettingAll = getLinesettingInv(equationType,inversionType,parameterisation,exchangeStrategy,showMarker); 
%lineSettingAll{2}.Color='k';
titleName = '';
[titleLabelSettingAll] = getTitleLabelSettingAll(showTitle,showXlabel,showYlabel...
    ,inversionType,parameterisation,exchangeStrategy,titleName,equationType,legendType,showText,textName);
offset.value=[1:Nshot];
offset.labelName='Source number';
residualScale=1;
syntheticDatanameAll={filenameTrue;filenameInv};
syntheticDataAll = {seismogramTrue;seismogramInv};
[imagename] = getImagenameFromfileList(syntheticDatanameAll);
[h answer] = plotSeismogramData_multi(DT,showmax,...
    syntheticDataAll,imagesave,orientation,lineSettingAll,legendNameAll,offset,...
    titleLabelSettingAll,showLegend,showResidual,normalize,wiggleType...
    ,imageScale,imagename,residualScale);
 
% plot spectra
FC=config.getValue('CenterFrequencyCPML');
NT=size(seismogramInv,2); NT2=2*NT-1;
Fs = 1/DT; f = (0:NT-1)*Fs/(NT*2); t=[0:(NT-1)]*DT;
for ir=1:Nshot
    SeismogramTrace = fft(seismogramTrue(ir,:),NT2);
    SeismogramTrue(ir,:) = SeismogramTrace(1:NT);
    SeismogramTrace = fft(seismogramInv(ir,:),NT2);
    SeismogramInv(ir,:) = SeismogramTrace(1:NT);
end
SeismogramTrue=mean(SeismogramTrue,1);
SeismogramInv=mean(SeismogramInv,1);
maxValue = max(abs(SeismogramTrue));
fshowIndex = [];
for it=1:NT
    if abs(SeismogramTrue(it)) > 0.01*maxValue
        fshowIndex = [fshowIndex it];
    end
end
flim=[0 f(fshowIndex(end))];
SeismogramInv=SeismogramInv/max(abs(SeismogramInv))*max(abs(SeismogramTrue));
figure
isSeismic = checkEquationType(equationType);
if isSeismic == 0
    X0=f*1e-6;
    flim=flim*1e-6;
    ylabelName='Frequency (MHz)';
else
    X0=f;
    ylabelName='Frequency (Hz)';
end
plot(X0(1:fshowIndex(end)),abs(SeismogramTrue(1:fshowIndex(end))),lineSettingAll{1}); hold on;
plot(X0(1:fshowIndex(end)),abs(SeismogramInv(1:fshowIndex(end))),lineSettingAll{2}); hold on;
xlabel(ylabelName);
ylabel('Amplitude');
legend('Source true','Source inverted','Location','best');
legend('boxoff')
axis([flim 0 maxValue*1.1]);
ImageBold(gca);
pause(0.1);

%% write
if writefiles~=0
    answer=overwritedlg(filenameInv,fileFormat);
    if answer==0
        return;
    end
    writeSeismogram(filenameInv,seismogramInv,fileFormat);
    if writeSourceTrue~=0  
        answer=overwritedlg(filenameTrue,fileFormat);
        if answer==0
            return;
        end
        writeSeismogram(filenameTrue,seismogramTrue,fileFormat);
    end
end

if copy_inv == 1
    answer = copyfile_image(filenameInv,fileFormat,DIR_PATH_NEW,Insert_name,config); 
    if answer==0
        return;
    end 
end
if copy_true_start == 1
    answer = copyfile_image_forward(filenameTrue,fileFormat,DIR_PATH_NEW,Insert_name,config);  
    if answer==0
        return;
    end 
end

addpath('../configuration');
addpath('../common');