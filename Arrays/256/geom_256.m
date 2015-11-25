function [G, GroupNum] = geom_256(n)
GroupNum = 4;
G = zeros(64,2);
for i=1:64
    C = AmpAndPinToLocation(n,i);
    G(i,:) = [C(1),C(2)];
end