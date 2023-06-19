function [obj_val, run_time, ret_code] =  lqdr(data)  
yalmip('clear')

myset = sdpsettings;
myset.verbose                      =  1;
myset.solver                       = 'mosek';
myset.debug                        = 1;
myset.savesolveroutput             = 1;
myset.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 7200;

% myset.mosek.MSK_IPAR_INTPNT_OFF_COL_TRH = 0;

T = data.T;
F = data.F;
G = data.G;
f = data.f;
s0 = data.s0;

m = size(F,1);
k = size(F,2);
n = m - 1;
K = k*T;
 
Umin =  -ones(k,1);
Umax =   ones(k,1);

PP = [ eye(k) ; 
      -eye(k) ;
           G 
           ];    
       
p = [ -Umin  ; 
       Umax  ;
       (k-1)*ones(2^k,1) 
       ]; 
   
P = [ kron(eye(T),PP),kron(ones(T,1),p) ];
L = size(P,1);

eKp1 = zeros(K+1,1);
eKp1(K+1) = 1;

ll = sdpvar(1,1);
u0 = sdpvar(n,1);

N  = sdpvar(L,L,'symm');
V  = sdpvar(K+1,K+1,'symm');

piWT  = sdpvar(T,1,'full');
NWT   = sdpvar(L,L,T,'symm');
VWT   = sdpvar(K+1,K+1,T,'symm');

piWTm = sdpvar(T,1,'full');
NWTm  = sdpvar(L,L,T,'symm');
VWTm  = sdpvar(K+1,K+1,T,'symm');

piUkT = sdpvar(n,T,'full');
NUkT  = sdpvar(L,n,T,'full');

piSUT = sdpvar(T,1,'full');
NSUT  = sdpvar(L,L,T,'symm');
VSUT  = sdpvar(K+1,K+1,T,'symm');

idx_u = n*ones(1,T);

U = sdpvar(idx_u, k*(1:T)+1);
W = sdpvar(k*(1:T)+1, k*(1:T)+1);
S = sdpvar(k*(1:T)+1, k*(1:T)+1);

if T == 1
    UU = cell(1,1);
    WW = cell(1,1);
    UU{1,1} = U;
    WW{1,1} = W;
    SS{1,1} = S;
    U = UU;
    W = WW;
    S = SS;
end

obj  = ll;
cons = [ ];
Mtx = 0;

cons = [ cons ; u0 >= 0 ];
cons = [ cons ; sum(u0) <= s0 ];

for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ]; 
    Mtx = Mtx + PJt'*W{1,t}*PJt;
end

cons = [ cons ; ll*(eKp1)*eKp1' - Mtx == P'*N*P + V ];
cons = [ cons ; V >= 0 ; N(:) >= 0 ];

%% for w_t(xi^t) >= xi_t - s_t(xi^t) and  w_t(xi^t) >= -xi_t + s_t(xi^t)
for t = 1 : T
    PJt = [ eye(k*t), zeros(k*t,K+1-k*t); zeros(1,K),1 ];
    FF = [ zeros((t-1)*k,1); F(m,:)'; zeros(K-k*t,1); f(m) ];
    
    cons = [ cons ; PJt'*W{1,t}*PJt + PJt'*S{1,t}*PJt - 0.5*(eKp1*FF' + FF*eKp1') - piWT(t)*(eKp1)*eKp1' == P'*NWT(:,:,t)*P + VWT(:,:,t); VWT(:,:,t) >= 0; piWT(t) >= 0 ];
    for l = 1 : L
        cons = [ cons ; NWT(:,l,t) >= 0 ];
    end
    
    cons = [ cons ; PJt'*W{1,t}*PJt - PJt'*S{1,t}*PJt + 0.5*(eKp1*FF' + FF*eKp1') - piWTm(t)*(eKp1)*eKp1' == P'*NWTm(:,:,t)*P + VWTm(:,:,t) ; VWTm(:,:,t) >= 0; piWTm(t) >= 0 ];
    for l = 1 : L
        cons = [ cons ; NWTm(:,l,t) >= 0 ];
    end
end

%% for u_t(xi^t) >= 0
for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ];
    for i = 1 : n
        cons = [ cons ; PJt'*U{1,t}(i,:)' == piUkT(i,t)*eKp1 + P'*NUkT(:,i,t) ; piUkT(i,t) >= 0 ; NUkT(:,i,t) >= 0 ];
    end
end

%% for s_t(xi^t) >= e'*u_t(xi^t)
for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ];
    
    Mtx = 0;
    for i = 1 : n
        Mtx = Mtx + 0.5*(PJt'*U{1,t}(i,:)'*eKp1' + eKp1*U{1,t}(i,:)*PJt);
    end
    
    cons = [ cons ; PJt'*S{1,t}*PJt - Mtx - piSUT(t)*(eKp1)*eKp1' == P'*NSUT(:,:,t)*P + VSUT(:,:,t) ; VSUT(:,:,t) >= 0; piSUT(t) >= 0 ];
    
    for l = 1 : L
        cons = [ cons ; NSUT(:,l,t) >= 0 ];
    end
end

for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ];
    PJtm1 = [ eye(k*(t-1)),zeros(k*(t-1),K+1-k*(t-1));zeros(1,K),1 ];
    
    Mtx = 0;
    for i = 1 : n
        FFit = [ zeros((t-1)*k,1); F(i,:)'; zeros(K-k*t,1); f(i) ];
        if t <= 1
            Mtx = Mtx + 0.5*u0(i)*(FFit*eKp1' + eKp1*FFit');
        else
            Mtx = Mtx + 0.5*(FFit*U{1,t-1}(i,:)*PJtm1 + PJtm1'*U{1,t-1}(i,:)'*FFit');
        end
    end
    
    cons = [ cons ; PJt'*S{1,t}*PJt - Mtx == 0 ];

end

diagnostics = solvesdp(cons,obj,myset);

% ret_code = diagnostics.problem;
run_time = diagnostics.solvertime;

primalobj = diagnostics.solveroutput.res.sol.itr.pobjval;
dualobj   = diagnostics.solveroutput.res.sol.itr.dobjval;
      
abs_obj_gap = abs(dualobj - primalobj);
rel_obj_gap = abs_obj_gap/(1.0 + min( abs(primalobj), abs(dualobj)));
  
accepted = 1;

if ( rel_obj_gap > 1e-6 )
    fprintf('Warning: The relative objective gap is LARGE.\n');
    accepted = 0;
end

if accepted == 1
    ret_code = 0;
else
    ret_code = 9;
end
      
obj_val = double(primalobj);

end




