
clear
clc

% time stage set
Tset = [ 1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
% number of instances generated for each time stage
NN  = 25;

TL = length(Tset);

for ts_iter = 1 : TL
    for seed_iter = 1 : NN  
        run_decision_rules( ts_iter, seed_iter );
    end
end