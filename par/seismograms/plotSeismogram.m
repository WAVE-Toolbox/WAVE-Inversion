%% Function Name: plotSeismogram
%
% Inputs:
%   DT(float): time sampling intervall in seconds
%
%   iteration(int): Nr. of the iterarion to be displayed
%
%   shot(int): Nr. of the shot gather to be displayed
%
%   component(string): receiver component eg. p,vx,vy,vz
%
%   skiptraces(int): display only every skiptraces trace (skiptraces=1 for
%   all traces)
%
%   syntheticData(string): location of the inverted data
%
%   syntheticData(string): location of the field data
%
%
% Outputs:
%   None

%
% $Date: March 21, 2018
% ________________________________________
function plotSeismogram(DT,stage,iteration,shot,component,skipTraces,syntheticData,fieldData)



%% Read seismogram
filename=[fieldData '.stage_' num2str(stage) '.shot_' num2str(shot) '.' component '.mtx'];
seismogramtrue=readSeismogram(filename);

T=1*DT:DT:size(seismogramtrue,2)*DT;
%% Plot seismogram
figure
for trace=1:skipTraces:size(seismogramtrue,1)
plot(T,seismogramtrue(trace,:)/max(abs(seismogramtrue(trace,:)))+trace,'black');
hold on
end
title(['field data (black) and modeled data (red): iteration ' num2str(iteration) ' shot ' num2str(shot)])
xlabel('Time in seconds')
ylabel('Traces')
axis([0.0 size(seismogramtrue,2)*DT 0 size(seismogramtrue,1)+1])


%% Read seismogram
filename=[syntheticData '.stage_' num2str(stage) '.It_' num2str(iteration) '.shot_' num2str(shot) '.' component '.mtx'];
seismogram=readSeismogram(filename);

%% Plot seismogram

for trace=1:skipTraces:size(seismogram,1)
plot(T,seismogram(trace,:)/max(abs(seismogram(trace,:)))+trace,'red');
hold on
end







