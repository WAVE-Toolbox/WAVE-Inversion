clearvars; close all;

OUTPUT_FILENAME='receiver.mtx';

%% Requiered Parameters
X=[70]; % Coordinates in X (Grid points)
Y=[70]; % Coordinates in Y (Grid points) (depth)
Z=[70]; % Coordinates in Z (Grid points)
RECEIVER_TYPE=[1]; % RECEIVER Type (1=P, 2=vX, 3=vY, 4=vZ)

%% Write to file

% Create Matrix
RECEIVER_FILE=[X' Y' Z' RECEIVER_TYPE'];

% Write mtx file
writeMatrix2mtx(OUTPUT_FILENAME,RECEIVER_FILE);