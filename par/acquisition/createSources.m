clearvars; close all;

OUTPUT_FILENAME='sources.txt';
NSOURCES=18;
SOURCE_FILE=zeros(NSOURCES,10);
for ii=1:NSOURCES
%% Requiered Parameters
shotNo=[ii-1];
X=[ii*4+11]; % Coordinates in X (Grid points)
Y=[49]; % Coordinates in Y (Grid points) (depth)
Z=[0]; % Coordinates in Z (Grid points)
SOURCE_TYPE=[1]; % Source Type (1=P,2=vX,3=vY,4=vZ)
WAVELET_TYPE=[1]; % Wavelet Type (1=Synthetic)

%% Optional Parameters
WAVELET_SHAPE=[3]; % Wavelet Shape (1=Ricker,2=Sinw,3=sin^3,4=FGaussian,5=Spike/Delta,6=integral sin^3)
FC=[10]; % Center Frequency in Hz
AMP=[1]; % Amplitude
TShift=[1*10^-1]; % Time shift in s



%% Write to file

% Create Matrix
SOURCE_FILE(ii,:)=[shotNo' X' Y' Z' SOURCE_TYPE' WAVELET_TYPE' WAVELET_SHAPE' FC' AMP' TShift'];
end

% Create header
headerline=["#shot_number","source_coordinate_(x)","source_coordinate_(y)",...
    "source_coordinate_(z)","source_type","wavelet_type","wavelet_shape",...
    "center_frequency","amplitude","time_shift"];

% Write txt file
fileID = fopen(OUTPUT_FILENAME,'w');
fprintf(fileID,'%s %s %s %s %s %s %s %s %s %s\n',headerline);
fprintf(fileID,'%i %i %i %i %i %i %i %g %g %g\n',SOURCE_FILE');
fclose(fileID);
%dlmwrite(OUTPUT_FILENAME,SOURCE_FILE,'delimiter',' ')
% writematrix(SOURCE_FILE,OUTPUT_FILENAME,'Delimiter','space') % introduced
% in Matlab 2019