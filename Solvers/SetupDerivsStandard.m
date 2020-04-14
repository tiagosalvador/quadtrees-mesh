function [Dxx,Dyy,Dxy] = SetupDerivsStandard(numBdy,numInt,numTot,ind,axx,axy,ayy)

ind(ind==0) = 1;
% Matrix for Dxx.
Dxx = sparse(repmat((1:numInt)',1,size(ind,2)),ind,axx,numInt,numTot) - ...
    sparse(1:numInt,numBdy+(1:numInt),sum(axx,2),numInt,numTot);

% Matrix for Dyy.
Dyy = sparse(repmat((1:numInt)',1,size(ind,2)),ind,ayy,numInt,numTot) - ...
    sparse(1:numInt,numBdy+(1:numInt),sum(ayy,2),numInt,numTot);

% Matrix for Dxy.
Dxy = sparse(repmat((1:numInt)',1,size(ind,2)),ind,axy,numInt,numTot) - ...
    sparse(1:numInt,numBdy+(1:numInt),sum(axy,2),numInt,numTot);

% % Matrix for Dxx.
% Dxx = sparse(repmat((1:numInt)',1,9),ind,axx,numInt,numTot) - ...
%     sparse(1:numInt,numBdy+(1:numInt),sum(axx,2),numInt,numTot);
% 
% % Matrix for Dyy.
% Dyy = sparse(repmat((1:numInt)',1,9),ind,ayy,numInt,numTot) - ...
%     sparse(1:numInt,numBdy+(1:numInt),sum(ayy,2),numInt,numTot);
% 
% % Matrix for Dxy.
% Dxy = sparse(repmat((1:numInt)',1,9),ind,axy,numInt,numTot) - ...
%     sparse(1:numInt,numBdy+(1:numInt),sum(axy,2),numInt,numTot);