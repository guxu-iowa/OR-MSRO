function [  ] = run_decision_rules( ts_iter, seed_iter )

Tset = [ 1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
% number of instances generated for each time stage

seed = ts_iter*100 + seed_iter;
rng(seed);

P = 4;
k = 4;
cb = 0.2;
ch = 0.2;
rD = 1;
rP = 1;
I0 = 3*ones(P,1);
Imx = 24;


T = Tset(ts_iter);
if T <= 3
    run_time = zeros(5,1);
    return_code = zeros(5,1);
    obj_val = zeros(5,1);
else
    run_time = zeros(3,1);
    return_code = zeros(3,1);
    obj_val = zeros(3,1);
end

alp = 2*rand(P,k) - 1;
bet = 2*rand(P,k) - 1;

data.T = T;
data.alp = alp;
data.bet = bet;
data.cb = cb;
data.ch = ch;
data.rD = rD;
data.rP = rP;
data.I0 = I0;
data.Imx = Imx;

foldername = sprintf('results/%d_%d', T, seed_iter);
mkdir(foldername);
filename = sprintf('results/%d_%d/Stage_%d_Seed_%d.mat', T, seed_iter, T, seed_iter);

fprintf('Now we are solving instance %d at time stage %d by using approximate S lemma method \n', seed_iter, Tset(ts_iter));
yalmip('clear');
[val_sl, tm_sl, rt_sl] = approxSL(data);
run_time(1) = tm_sl;
return_code(1) = rt_sl;
obj_val(1) = val_sl;

save(filename,'run_time','return_code','obj_val');
fprintf('Now we are done with solving instance %d at time stage %d by using approximate S lemma method \n', seed_iter, Tset(ts_iter));


fprintf('Now we are solving instance %d at time stage %d by using our LDR \n', seed_iter, Tset(ts_iter));
yalmip('clear');
[val_ldr, tm_ldr, rt_ldr] = ldr(data);
run_time(2) = tm_ldr;
return_code(2) = rt_ldr;
obj_val(2) = val_ldr;
save(filename,'run_time','return_code','obj_val');
fprintf('Now we are done with solving instance %d at time stage %d by using our LDR \n', seed_iter, Tset(ts_iter));

fprintf('Now we are solving instance %d at time stage %d by using our piecewise LDR \n', seed_iter, Tset(ts_iter));
yalmip('clear');
[val_piece, tm_piece, rt_piece] = piecewiseLDR(data);
run_time(3) = tm_piece;
return_code(3) = rt_piece;
obj_val(3) = val_piece;
save(filename,'run_time','return_code','obj_val');
fprintf('Now we are done with solving instance %d at time stage %d by using our piecewise LDR \n', seed_iter, Tset(ts_iter));

if T <= 3
    fprintf('Now we are solving instance %d at time stage %d by using PDR 2 \n', seed_iter, Tset(ts_iter));
    yalmip('clear');
    [val_pdr2, tm_pdr2, rt_pdr2] = pdr(data,2);
    run_time(4) = tm_pdr2;
    return_code(4) = rt_pdr2;
    obj_val(4) = val_pdr2;
    save(filename,'run_time','return_code','obj_val');
    fprintf('Now we are done with solving instance %d at time stage %d by using PDR 2 \n', seed_iter, Tset(ts_iter));

    fprintf('Now we are solving instance %d at time stage %d by using PDR 3 \n', seed_iter, Tset(ts_iter));
    yalmip('clear');
    [val_pdr3, tm_pdr3, rt_pdr3] = pdr(data,3);
    run_time(5) = tm_pdr3;
    return_code(5) = rt_pdr3;
    obj_val(5) = val_pdr3;
    save(filename,'run_time','return_code','obj_val');
    fprintf('Now we are done with solving instance %d at time stage %d by using PDR 3 \n', seed_iter, Tset(ts_iter));
end


end

