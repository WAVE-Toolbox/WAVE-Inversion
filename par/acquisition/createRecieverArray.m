clearvars; close all;

OUTPUT_FILENAME='receiver.mtx';

%% Requiered Parameters

% Grid information
NX=100; % Total number of grid points in X
NY=100; % Total number of grid points in Y (depth)
NZ=0; % Total number of grid points in Z

% Define receiver array
X1=80; X2=80; % X1: First receiver pos in X, X2: Last receiber pos in X
Y1=10; Y2=89; % Y1: First receiver pos in Y, Y2: Last receiber pos in Y
Z1=0; Z2=0; % Z1: First receiver pos in Z, Z2: Last receiber pos in Z
NGEOPH=80; % Place a geophone on each GP (NGEOPH=1) or every other (NGEOPH=2)


% Define receiver type
RECEIVER_TYPE=[1]; % (1=P, 2=vX, 3=vY, 4=vZ), if you want all write [1 2 3 4]


%% Create receiver array
X=ones(1,NGEOPH);
Y=ones(1,NGEOPH);
Z=ones(1,NGEOPH);

for ng=1:NGEOPH
X(ng)=floor(X1+(X2-X1)/(NGEOPH-1)*(ng-1));
Y(ng)=floor(Y1+(Y2-Y1)/(NGEOPH-1)*(ng-1));
Z(ng)=floor(Z1+(Z2-Z1)/(NGEOPH-1)*(ng-1));
end
% Y=80*ones(size(X));
% Z=zeros(size(X));

numel_X=numel(X);
RECEIVER_TYPE_vec=RECEIVER_TYPE(1)*ones(1,numel_X);
for(i=2:size(RECEIVER_TYPE,2))
    if(size(RECEIVER_TYPE,2)>1)
        X=[X X(1:numel_X)]; Y=[Y Y(1:numel_X)]; Z=[Z Z(1:numel_X)];
        RECEIVER_TYPE_vec=[RECEIVER_TYPE_vec RECEIVER_TYPE(i)*ones(1,numel_X)];
    end
end


%% Write to file

% Create Matrix
RECEIVER_FILE=[X' Y' Z' RECEIVER_TYPE_vec'];

% Write mtx file
writeMatrix2mtx(OUTPUT_FILENAME,RECEIVER_FILE);