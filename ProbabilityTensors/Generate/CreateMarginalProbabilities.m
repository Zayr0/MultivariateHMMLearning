function [P1, P21, P31, P32, P321, P312] = CreateMarginalProbabilities(seq, d)
[Y, ~] = discretize(seq, d);

P1 = zeros(d, 1);
P21 = zeros(d, d);
P31 = zeros(d, d);
P32 = zeros(d, d);
P321 = zeros(d,d,d);
P312 = zeros(d,d,d);

for t = 1:(length(Y) - 3)
    I1 = Y(t);
    I2 = Y(t+1);
    I3 = Y(t+2);
    
    P1(I1) = P1(I1) + 1;
    P21(I2, I1) = P21(I2, I1) + 1;
    P31(I3, I1) = P31(I3, I1) + 1;
    P32(I3, I2) = P32(I3, I2) + 1;
    P321(I3, I2, I1) = P321(I3, I2, I1) + 1;
    P312(I3, I1, I2) = P312(I3, I1, I2) + 1;
end

P1 = P1 ./ sum(P1);
P21 = P21 ./ sum(P21, 2);
P31 = P31 ./ sum(P31, 2);
P32 = P32 ./ sum(P32, 2);

P321 = P321 ./ sum(P321, 2);
P312 = P312 ./ sum(P312, 2);
end