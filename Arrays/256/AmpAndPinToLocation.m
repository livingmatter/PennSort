function C=AmpAndPinToLocation(i,j)
load('IndexMatrix1.mat');
load('SpatialXY.mat');
B=Amp{i};
for k=1:16
    for l=1:16
        if B(k,l)==j
            C=[X(k,l) Y(k,l) k l];
           
        end
    end
end
