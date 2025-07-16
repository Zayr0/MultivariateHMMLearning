d = 3;


P = rand(d, d, d);
P = P ./sum(P(:));

P_cumsum = cumsum(vec(P));

T_seq = 10;

for t = 1:T_seq
    [x, y, z] = index_2_3D(find(P_cumsum >= rand, 1, "first"), d)
    % index = index_2_3D(indices(1));
end


function [x, y, z] = index_2_3D(index, d)
 z = mod(index, d);
 y = mod(rem(index, d), d);
 x = rem(index, (d * d)); 
end