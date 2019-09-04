function [vector]=readVectorfromLMF(filename)


fileID = fopen(filename,'r');

% read the header [ID, indexType, valueType]
HEADER = fread(fileID, [1 3], 'int');
if HEADER(1) ~= hex2dec( '4711E01' )
   error('HEADER=%s not 4711E01 (dense vector)', dec2hex(HEADER(1)))
end
if HEADER(2) ~= 0
   error('data type for index values not int')
end
if HEADER(3) ~= 2
   error('data type for values not float')
end

NDIMS  = fread(fileID, [1 1], 'int');
if NDIMS~=1
    error('readModelFromLMF, NDIMS=%d must be 1', NDIMS);
end

SIZE = fread(fileID, [1 1], 'int');

vector = fread(fileID, [1 SIZE], 'float');
fclose(fileID);

end
