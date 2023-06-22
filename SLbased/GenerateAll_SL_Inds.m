%% Returns a matrix containing the linear indices for all voxels within 
% the specified radius (radiusInVox) for each linear index (/voxel
% location) which is 1 in the provided mask. The dimensions of the space
% are derived from the mask. Invalid indices and  

% This will use up quite some memory. If this is
% an issue, use the GenerateSingle_SL_Ind function 

function [SLindsLinear, LinearInds, NrOfValidVoxelsPerSL]=GenerateAll_SL_Inds(radiusInVox,  mask)
% radiusInVox = 4;
% XYZdims=[100,100,100];
% mask=ones(XYZdims)==1;
XYZdims=size(mask);
[x,y,z] = ndgrid(1:XYZdims(1),1:XYZdims(2),1:XYZdims(3));
refCoord = round([XYZdims(1)/2,XYZdims(2)/2,XYZdims(3)/2]);
refCoordLin = sub2ind( XYZdims, refCoord(1), refCoord(2), refCoord(3) );
linArray = find(sqrt((x-refCoord(1)).^2+(y-refCoord(2)).^2+(z-refCoord(3)).^2) <= radiusInVox)-refCoordLin;
clear x y z;
LinearInds = find( mask == 1 );
SLindsLinear=repmat(single(linArray),1,numel(LinearInds))+...
    repmat(single(LinearInds'),numel(linArray),1);
SLindsLinear(SLindsLinear<1|SLindsLinear>numel(mask))=NaN;
SLindsLinear(ismember(SLindsLinear,LinearInds)==0)=NaN;
NrOfValidVoxelsPerSL=single(sum(isnan(SLindsLinear)==0));