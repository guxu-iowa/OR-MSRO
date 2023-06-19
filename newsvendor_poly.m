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

function [poly_obj, poly_x, run_time] = newsvendor_poly(A, B, d, F, f, P, p, DEG)

myset = sdpsettings;
myset.verbose                      =  1;
myset.solver                       = 'mosek';
myset.MaxFunEvals                  = 10000;
myset.debug                        = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set base dimension

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization model preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up the degree of the polynomial decision rule
% DEG = 2;

n1 = size(A,2);
n2 = size(B,2);
m = size(F,1);
k = size(F,2);
s = size(P,1);

dim_poly = nchoosek(k+DEG,DEG);
% degree of the SDP problem
dd = ceil(DEG/2);

dim_x = nchoosek(k+dd,dd);
dim_palpha = nchoosek(k+2*dd,2*dd);

dz = dd-1; % for polyhedral uncertainty set
dim_z = nchoosek(k+dz,dz);

ll = sdpvar(1);
x  = sdpvar(n1,1);

% polynomial policy with degree d
polyfunc = sdpvar(n2,dim_poly,'full');
palpha = sdpvar(m+1,dim_palpha,'full');

X = sdpvar(dim_x,dim_x,m+1,'symm');

% Z matrices correspond to constraints Pu + p >= 0
% degree of Pu + p is 1. 
% quadratic 
Z = sdpvar(dim_z,dim_z,s,m+1,'symm');

[Balpha, Calpha] = genBC(dd, p, P);

cons = [];

cons = [ cons; x >= 0 ];

for i = 1 : m+1
    cons = [ cons; X(:,:,i) >= 0 ];
    
    for j = 1 : s
        cons = [ cons; Z(:,:,j,i) >= 0 ];
    end
end

cons = [ cons; ll - d'*polyfunc(:,1) == palpha(1,1) ];

for j = 2 : dim_poly
    cons = [ cons; palpha(1,j) == -d'*polyfunc(:,j) ];
end

for j = dim_poly+1 : dim_palpha
    cons = [ cons; palpha(1,j) == 0 ];
end

for i = 1 : m
    cons = [ cons; palpha(i+1,1) == A(i,:)*x - f(i) + B(i,:)*polyfunc(:,1) ];
    for j = 2 : k+1
        cons = [ cons; palpha(i+1,j) == B(i,:)*polyfunc(:,j) - F(i,j-1) ]; 
    end
    for j = k+2 : dim_poly
        cons = [ cons; palpha(i+1,j) == B(i,:)*polyfunc(:,j) ];
    end
    for j = dim_poly+1 : dim_palpha
        cons = [ cons; palpha(i+1,j) == 0 ];
    end
end


for i = 1 : m+1
    tmp_mtx = sparse(Balpha(:,:,1));
    tmp = palpha(i,1) - trace(tmp_mtx*X(:,:,i));
%     tmp = palpha(i,1) - trace(Balpha(:,:,1)*X(:,:,i));
    for j = 1 : s
        tmp_mtx = sparse(Calpha(:,:,j,1));
        tmp = tmp - trace(Calpha(:,:,j,1)*Z(:,:,j,i));
%         tmp = tmp - trace(Calpha(:,:,j,1)*Z(:,:,j,i));
    end
    cons = [ cons; tmp >= 0 ];
    
    for j = 2 : dim_palpha
        tmp_mtx = sparse(Balpha(:,:,j));
        tmp = trace(tmp_mtx*X(:,:,i));
%         tmp = trace(Balpha(:,:,j)*X(:,:,i));
        for l = 1 : s
            tmp_mtx = sparse(Calpha(:,:,l,j));
            tmp = tmp + trace(tmp_mtx*Z(:,:,l,i));
%             tmp = tmp + trace(Calpha(:,:,l,j)*Z(:,:,l,i));
        end
        cons = [ cons; tmp == palpha(i,j) ];
    end
end

obj = ll;
diagnostics = solvesdp(cons,obj,myset);
run_time = diagnostics.solvertime;  

poly_x = double(x);

poly_obj = -double(obj);


end




   