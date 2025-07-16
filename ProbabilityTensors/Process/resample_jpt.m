function [P3_n] = resample_jpt(P3_inf, T_seq)

P_cumsum= vec(P3_inf);

for t = 1:T_seq
    val = rand;
    index = find(P_cumsum > rand)
end

end