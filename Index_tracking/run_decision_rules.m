function [  ] = run_decision_rules( ts_iter, ins_num )

Tset = [ 1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
% number of instances generated for each time stage

k = 3;
n = 4;
m = n+1;
f = ones(m,1);
s0 = 1.0;

strmat = dec2bin(0:(2^k - 1));
  
G = ones(2^k, k); 
for i = 1:size(strmat, 1)
    for j = 1:k
        if strmat(i, j) == '0'
            G(i, j) = -1;
        end
    end
end

T = Tset(ts_iter);

if T <= 3
    run_time = zeros(3,1);
    return_code = zeros(3,1);
    obj_val = zeros(3,1);
else
    run_time = zeros(2,1);
    return_code = zeros(2,1);
    obj_val = zeros(2,1);
end

data.T = T;
data.G = G;
data.f = f;
data.s0 = s0;


cnt = 1;
seed_cnt = 1;
while cnt <= ins_num
    seed_cnt = seed_cnt + 1 ;
    seed = ts_iter*1000 + seed_cnt;
    rng(seed);
    F = (2*rand(m,k)-1)/k;
    data.F = F;
    
    [val_sl, tm_sl, rt_sl] = approxSL(data);
    run_time(1) = tm_sl;
    return_code(1) = rt_sl;
    obj_val(1) = val_sl;
    if rt_sl ~= 0
        continue;
    end
    
    [val_ldr, tm_ldr, rt_ldr] = lqdr(data);
    run_time(2) = tm_ldr;
    return_code(2) = rt_ldr;
    obj_val(2) = val_ldr;
    if rt_ldr ~= 0
        continue;
    end
    
    if T <= 3
        [val_pdr3, tm_pdr3, rt_pdr3] = pdr(data,3);
        run_time(3) = tm_pdr3;
        return_code(3) = rt_pdr3;
        obj_val(3) = val_pdr3;
        if rt_pdr3 ~= 0
            continue;
        end
    end
    
    foldername = sprintf('results/%d_%d', T, cnt);
    mkdir(foldername);
    filename = sprintf('results/%d_%d/Stage_%d_Seed_%d.mat', T, cnt, T, cnt);
    save(filename,'run_time','return_code','obj_val');
    cnt = cnt + 1;
end

end

