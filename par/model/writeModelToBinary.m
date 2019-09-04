function writeModelToBinary(model,filename)
% writeModelToBinary  writes models to a binary (float,le) files

%   writeModelToBinary(filename,model)
%   filename: filename of the model
%   model: 
%   input must be model(Y,X,Z) to use the binary with seismic unix 
%   with Y(Depth): fast dimension
%   Transform to su with
%   suaddhead ns=NY < filename.bin | sushw key=gy,gx,dt a=0,0,DH b=DH,0,0 c=0,DH,0 j=NX,NX,0  > filename.su


fileID = fopen(filename,'w');
fwrite(fileID,model,'float32','ieee-le');
fclose(fileID);
end


