clearvars; close all;

DT=2e-3;

shot=5;
iteration=0;

fieldData='../ci/rectangle.true.shot_';
syntheticData='seismogram';


%% Read seismogram
filename=[fieldData num2str(shot) '.p.mtx'];
seismogramtrue=readSeismogram(filename);

T=1*DT:DT:size(seismogramtrue,2)*DT;
%% Plot seismogram
figure
for(trace=1:3:size(seismogramtrue,1))
plot(T,seismogramtrue(trace,:)/max(abs(seismogramtrue(trace,:)))+trace,'black');
hold on
end
title(['field data (black) and modeled data (red): iteration ' num2str(iteration) ' shot ' num2str(shot)])
xlabel('Time in seconds')
ylabel('Traces')
axis([0.7 size(seismogramtrue,2)*DT 0 size(seismogramtrue,1)+1])


%% Read seismogram
filename=[syntheticData '.It_' num2str(iteration) '.shot_' num2str(shot) '.p.mtx'];
seismogram=readSeismogram(filename);

%% Plot seismogram

for(trace=1:3:size(seismogram,1))
plot(T,seismogram(trace,:)/max(abs(seismogram(trace,:)))+trace,'red');
hold on
end







