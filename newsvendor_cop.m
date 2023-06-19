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
function [sdpval,sdpx,run_time] = newsvendor_cop(A, B, d, F, P)

myset = sdpsettings;
myset.verbose                      =  1;
myset.solver                       = 'mosek';
myset.MaxFunEvals                  = 10000;
myset.debug                        = 1;


% Take the indices
n1 = size(A, 2);
n2 = size(B, 2);
m = size(F, 1);
k = size(F, 2);
s = size(P, 1);

e1 = zeros(k, 1);
e1(1) = 1;

g1 = zeros(k+m, 1);
g1(1) = 1;
E = [ -d*e1', B' ];

x = sdpvar(n1, 1);
G = [ zeros(k), (F - A*x*e1')'; F - A*x*e1', zeros(m) ];

ll = sdpvar(1);

alp = sdpvar(k, 1);



LL = sdpvar(k+m, n2, 'full');

S11 = sdpvar(k, k, 'symm');
S21 = sdpvar(m, k, 'full');
S22 = sdpvar(m, m, 'symm');
S = [S11, S21'; S21, S22];

R11 = sdpvar(k, k, 'symm');
R21 = sdpvar(m, k, 'full');
R22 = sdpvar(m, m, 'symm');
R = [R11, R21'; R21, R22];
  
M = sdpvar(k + m, k + m, 'symm');

con = [];
con = [ con; x >= 0 ];

% Add the constraints, which are independent of unc_norm

con = [con; S22(:) >= 0];

con = [con; R21(:) == 0];
con = [con; R22(:) == 0];

con = [con; M >= 0];

con = [con; ll*(g1)*(g1)' - 0.5*G + 0.5*(E'*LL' + LL*E) == S + M + R];


N = sdpvar(s, s, 'symm');

% Add constraints for R11

con = [con; N(:) >= 0; R11 == P'*N*P];
    
tt = sdpvar(s, 1);
con = [con; S11 == (e1)*alp' + alp*(e1)'; alp == P'*tt; tt >= 0]; 

% Add polyhedral constraints

W = sdpvar(s, m);
con = [con; W(:) >= 0];
for i = 1:m
   con = [con; S21(i, :)' == P'*W(:,i)];
end
  
obj = ll;
diagnostics = solvesdp(con,obj,myset);
run_time = diagnostics.solvertime;

sdpx = double(x);
sdpx
sdpval = double(obj);

% The true value is negating sdpval
sdpval = -sdpval;

end