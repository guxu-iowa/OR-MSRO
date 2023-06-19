function [ obj_val, run_time, ret_code ] = piecewiseLDR( data )

myset = sdpsettings;
myset.verbose                      =  0;
myset.solver                       = 'mosek';
myset.debug                        = 1;
myset.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600;


T = data.T;
alp = data.alp;
bet = data.bet;
cb = data.cb;
ch = data.ch;
rD = data.rD;
rP = data.rP;
I0 = data.I0;
Imx = data.Imx;

P = size(alp,1);
m = size(alp,2);

% This parameter should change in response to the specific piecewise scheme
% n = m;
n = m ;

k = m+n;
K = k*T;

Q = n*T;

R = [ zeros(P,n),rP*alp,4*ones(P,1) ];
D = [ zeros(P,n),0.5*rD*bet,2*ones(P,1) ];

idx = 1 : T;

d = [ sin(2*pi*(idx-1)/12);
      sin(2*pi*(idx-1)/12);
      cos(2*pi*(idx-1)/12);
      cos(2*pi*(idx-1)/12) ];

  
Umin =  -ones(m,1);
Umax =   ones(m,1);
Wmax =   ones(m,1);
% Wmax =  [ ones(m,1); 2*ones(m*(m-1),1) ];


E = [ 1,  0,  0,  0  ;
      0,  1,  0,  0  ;
      0,  0,  1,  0  ;
      0,  0,  0,  1 ];
  
UU = [ zeros(m,n),  eye(m) ;
       zeros(m,n), -eye(m) ;
       eye(n),    zeros(n,m) ;
      -eye(n),    zeros(n,m) ;
       eye(n),    -E  ];
        
u = [ -Umin  ; 
       Umax  ;
       zeros(n,1) ;
       Wmax  ;
       zeros(n,1)];
   
U = [ kron(eye(T), UU), kron(ones(T,1), u) ];

L = size(U,1);

C = cell(1,Q);

for t = 1 : T
    for i = 1 : n 
        
        ff = zeros(n,1);
        ff(i) = 1;

        ee = -E(i,:)';
        
        zz1 = [ zeros((t-1)*k,1); ff; zeros(m,1); zeros(K+1-t*k,1) ];
        zz2 = [ zeros((t-1)*k,1) ; ff ; ee ; zeros(K+1-t*k,1) ];

        C{1,(t-1)*n+i} = 0.5*(zz1*zz2' + zz2*zz1');
    end
end

eKp1 = zeros(K+1,1);
eKp1(K+1) = 1;

ll  = sdpvar(1,1);

idx = P*ones(1,T);

X = sdpvar(idx,k*(1:T)+1);
Y = sdpvar(idx,k*(1:T)+1);
if T == 1
    XX = cell(1,1);
    YY = cell(1,1);
    XX{1,1} = X;
    YY{1,1} = Y;
    X = XX;
    Y = YY;
end


N  = sdpvar(L,L,'symm');
V  = sdpvar(K+1,K+1,'symm');
aa = sdpvar(Q,1,'full');

piYPT = sdpvar(P,T,'full');
NYPT  = sdpvar(L,P,T,'full');

piXPT = sdpvar(P,T,'full');
NXPT  = sdpvar(L,P,T,'full');

piBPT = sdpvar(P,T,'full');
NBPT  = sdpvar(L,P,T,'full');


piIPT = sdpvar(P,T,'full');
NIPT  = sdpvar(L,P,T,'full');


piIPTmx = sdpvar(P,T,'full');
NIPTmx  = sdpvar(L,P,T,'full');


obj  = ll;

cons = [ ];

Mtx = 0;

for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ]; 
    
    for p = 1 : P
        ee = [ zeros((t-1)*k,1); R(p,1:k)'; zeros(K-k*t,1); R(p,k+1) ];
        Mtx = Mtx - 0.5*(ee*X{1,t}(p,:)*PJt + PJt'*X{1,t}(p,:)'*ee');
        
        for s = 1 : t
            dd = [ zeros((s-1)*k,1); D(p,1:k)'; zeros(K-k*s,1); D(p,k+1) + d(p,s) ];
            PJs = [ eye(k*s),zeros(k*s,K+1-k*s);zeros(1,K),1 ]; 
            Mtx = Mtx + 0.5*cb*(dd*eKp1' + eKp1*dd' - eKp1*X{1,s}(p,:)*PJs - PJs'*X{1,s}(p,:)'*eKp1');      
        end
        
        Mtx = Mtx + ch*I0(p)*(eKp1)*eKp1';
        
        for s = 1 : t
            PJs = [ eye(k*s),zeros(k*s,K+1-k*s);zeros(1,K),1 ];
            Mtx = Mtx - 0.5*ch*(eKp1*X{1,s}(p,:)*PJs + PJs'*X{1,s}(p,:)'*eKp1');
            
            if s >= 2
                PJsm1 = [ eye(k*(s-1)),zeros(k*(s-1),K+1-k*(s-1));zeros(1,K),1 ];
                Mtx = Mtx + 0.5*ch*(eKp1*Y{1,s-1}(p,:)*PJsm1 + PJsm1'*Y{1,s-1}(p,:)'*eKp1');
            end
        end
    end
end

for q = 1 : Q
    Mtx = Mtx + aa(q)*C{1,q};
end
cons = [ cons ; ll*(eKp1)*eKp1' - Mtx == U'*N*U + V ];
cons = [ cons ; V >= 0 ; N(:) >= 0 ];

for t = 1 : T
    PJt = [ eye(k*t),zeros(k*t,K+1-k*t);zeros(1,K),1 ];
    for p = 1 : P
        cons = [ cons ; PJt'*X{1,t}(p,:)' == piXPT(p,t)*eKp1 + U'*NXPT(:,p,t) ; piXPT(p,t) >= 0 ; NXPT(:,p,t) >= 0 ];
        cons = [ cons ; PJt'*Y{1,t}(p,:)' == piYPT(p,t)*eKp1 +  U'*NYPT(:,p,t); piYPT(p,t) >= 0; NYPT(:,p,t) >= 0 ];
    end
end


for t = 1 : T
    for p = 1 : P
        vec = 0;
        for s = 1 : t
            PJs = [ eye(k*s),zeros(k*s,K+1-k*s);zeros(1,K),1 ];
            dd = [ zeros((s-1)*k,1); D(p,1:k)'; zeros(K-k*s,1); D(p,k+1) + d(p,s) ];
            vec = vec + dd - PJs'*X{1,s}(p,:)';
        end
        cons = [ cons ; vec == piBPT(p,t)*eKp1 + U'*NBPT(:,p,t) ; piBPT(p,t) >= 0 ; NBPT(:,p,t) >= 0 ];
    end
end

for t = 1 : T
    for p = 1 : P
        vec = 0;
        for s = 1 : t
            PJs = [ eye(k*s),zeros(k*s,K+1-k*s);zeros(1,K),1 ];
            
            vec = vec - PJs'* X{1,s}(p,:)';
            if s >= 2
                PJsm1 = [ eye(k*(s-1)),zeros(k*(s-1),K+1-k*(s-1));zeros(1,K),1 ];
                vec = vec + PJsm1'*Y{1,s-1}(p,:)';
            end
        end
                
        cons = [ cons ; vec == piIPT(p,t)*eKp1 + U'*NIPT(:,p,t); NIPT(:,p,t) >= 0; piIPT(p,t) + I0(p) >= 0  ];
        cons = [ cons ; -vec == piIPTmx(p,t)*eKp1 + U'*NIPTmx(:,p,t);  NIPTmx(:,p,t) >= 0; Imx - I0(p) + piIPTmx(p,t) >= 0 ];

    end
end

diagnostics = solvesdp(cons,obj,myset);

obj_val = double(obj);
run_time = diagnostics.solvertime;
ret_code = diagnostics.problem;


end

