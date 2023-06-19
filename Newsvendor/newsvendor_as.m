% This is Example 1 from the paper entitled "Linearized Robust
% Counterparts of Two-Stage Robust Optimization Problems with Applications
% in Operations Management by Amir Ardestani-Jaafari and Erick Delage.
%
% In particular, this is a multi-item newsvendor problem:
% max_{x \in X} min_{zeta \in U} \sum_j r_j min(x_j, zeta_j) - c_jx_j + s_j max(x_j - zeta_j, 0) - p_j max(zeta_j - x_j, 0),
% where r_j, c_j, s_j <= r_j, and p_j denote sale price, ordering cost,
% salvage price, and shortage cost of a unit of the j-th item, j \in J,
% respectively, and zeta_j denotes the demand for item j for each j. 
%
% This problem can be reformulated to a two-stage adjustable robust linear
% optimization problem:
% max_{x \in X, y} min_{zeta \in \U} sum_j y_j(zeta)
% s.t.             y_j(zeta) <= (r_j - c_j)x_j - (r_j - s_j)(x_j - zeta_j)  \forall j \in J, \forall zeta \in \U
%                  y_j(zeta) <= (r_j - c_j)x_j - p_j(zeta_j - x_j) \forall j \in J, \forall zeta \in \U
% 
% The paper consider the following instance:
% 3 items, r = [80; 80; 80]; c = [70; 50; 20]; s = [20; 15; 10], and p = [60; 60; 50]
% Demand vector zeta is defined in the following uncertainty set U:
% 
% U = { zeta : exit deltap >= 0, 
%              deltan >= 0, deltap + deltan <= e,
%              e^T*deltap + e^T*deltan = Gamma, 
%              zeta_1 = zbar_1 + zhat_1*(deltap_1 + deltap_2 - deltan_1 - deltan_2)/2,
%              zeta_2 = zbar_2 + zhat_2*(deltap_2 + deltap_3 - deltan_2 - deltan_3)/2,
%              zeta_3 = zbar_3 + zhat_3*(deltap_3 + deltap_1 - deltan_3 - deltan_1)/2
%     }
% where Gamma = 2, zbar = [80; 80; 60], and zhat = [60; 60; 40].
% 

function [sdpval, sdpx, run_time] = newsvendor_as(A, B, d, F, P)


myset = sdpsettings;
myset.verbose                      =  1;
myset.solver                       = 'mosek';
myset.MaxFunEvals                  = 10000;
myset.debug                        = 1;


% Take the indices
n1 = size(A,2);
n2 = size(B,2);
m  = size(F,1);
k  = size(F,2);
s  = size(P,1);

e1 = zeros(k, 1);
e1(1) = 1;


%% The following part is for using the quadratic decision rules:
%  We essentially solve a semidefinite-programming based problem,
%  which is computationally tractable.
pi = sdpvar(m,1);

Q  = sdpvar(k,k,n2,'symm');
% N  = sdpvar(s, s, 'symm');
Tht = sdpvar(s, 1);
S  = sdpvar(k, k, 'symm');
tau = sdpvar(1,1);

% NN = sdpvar(s, s, m, 'symm');
SS = sdpvar(k, k, m, 'symm');
Thtm = sdpvar(s, m, 'full');
taum = sdpvar(s,1);

x = sdpvar(n1, 1);
ll = sdpvar(1);

L = -eye(k);
L(1,1) = k*k+1;

con = [];
con = [ con; x >= 0 ];

T1 = zeros(k);

for j = 1 : n2
    T1 = T1 + d(j)*Q(:,:,j);
end
con = [ con; ll*(e1)*(e1)' - T1 == S + 0.5*(P'*Tht*e1' + e1*Tht'*P) + tau*L; S >= 0; Tht >= 0; tau == 0 ];

con = [ con; pi >= 0 ];

for i = 1 : m
    Ti = zeros(k);
    for j = 1 : n2
        Ti = Ti + B(i,j)*Q(:,:,j);
        
    end
    con = [ con; (A(i,:)*x-pi(i))*(e1)*(e1)' + Ti - 0.5*F(i,:)'*(e1)' - 0.5*(e1)*F(i,:) == SS(:,:,i) + 0.5*(P'*Thtm(:,i)*e1' + e1*Thtm(:,i)'*P) + taum(i)*L; SS(:,:,i) >= 0; Thtm(:,i) >= 0; taum(i) == 0 ];
end


obj = ll;

diagnostics = solvesdp(con,obj,myset);
run_time = diagnostics.solvertime; 

sdpx = double(x);
sdpval = double(obj);

% The true value is negating sdpval
sdpval = -sdpval;

end
  