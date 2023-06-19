function [ obj_val, run_time, ret_code ] = pdr(data, deg)

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
k = size(alp,2);
K = k*T;


R = [ rP*alp,     4*ones(P,1) ];
D = [ 0.5*rD*bet, 2*ones(P,1) ];

idx = 1 : T;
d = [ sin(2*pi*(idx-1)/12);
      sin(2*pi*(idx-1)/12);
      cos(2*pi*(idx-1)/12);
      cos(2*pi*(idx-1)/12) ];
  
ll  = sdpvar(1,1);
dvar = ll;

u = sdpvar(K,1);

Svar = cell(P,T);
Ovar = cell(P,T);

for t = 1 : T
    dim_s = nchoosek(k*t+deg-1,deg-1);
    
    dim_o = nchoosek(k*t+deg,deg);
    for p = 1 : P
        Svar{p,t} = sdpvar(dim_s,1);
        Ovar{p,t} = sdpvar(dim_o,1);
        dvar = [ dvar; Svar{p,t}; Ovar{p,t}];
    end
end


r_func = cell(P,T);
d_func = cell(P,T);

s_policy = cell(P,T);
o_policy = cell(P,T);

for t = 1 : T
    for p = 1 : P
        r_func{p,t} = R(p,1:k)*u((t-1)*k+1:t*k) + R(p,k+1);
        d_func{p,t} = D(p,1:k)*u((t-1)*k+1:t*k) + D(p,k+1) + d(p,t);
        
        s_policy{p,t} = Svar{p,t}'*monolist(u(1:k*t)',deg-1);
        o_policy{p,t} = Ovar{p,t}'*monolist(u(1:k*t)',deg);
    end
end

cons_func = [ u + ones(K,1) ;
             -u + ones(K,1)];
      
L = size(cons_func,1);

obj = ll;

F = [ ];

palp = ll;
for t = 1 : T
    for p = 1 : P
        palp = palp + r_func{p,t}*s_policy{p,t} - ch*I0(p);
        
        for s = 1 : t
            palp = palp - cb*(d_func{p,s} - s_policy{p,s}) + ch*s_policy{p,s};
            if s >= 2
                palp = palp - ch*o_policy{p,s-1};
            end
        end
    end
end

sig_func = cell(1,L);
sig_var  = cell(1,L);

for l = 1 : L
    [sig_func{1,l}, sig_var{1,l}] = polynomial(u',deg-1);
    
    dvar = [ dvar; sig_var{1,l}];
    
    palp = palp - sig_func{1,l}*cons_func(l,1);
    F = [ F ; sos(sig_func{1,l}) ];
end

F = [ F ; sos(palp) ];


%% For S and O policies
sig_func_s = cell(L,P,T);
sig_var_s  = cell(L,P,T);
sig_func_o = cell(L,P,T);
sig_var_o  = cell(L,P,T);
sig_func_i = cell(L,P,T);
sig_var_i  = cell(L,P,T);
sig_func_ix = cell(L,P,T);
sig_var_ix  = cell(L,P,T);
sig_func_b = cell(L,P,T);
sig_var_b  = cell(L,P,T);

for t = 1 : T
    for p = 1 : P
        palp_s = s_policy{p,t};
        palp_o = o_policy{p,t};
        palp_i = I0(p);
        palp_ix = Imx-I0(p);
        palp_b = 0;
        
        for s = 1 : t
            palp_i = palp_i - s_policy{p,s};
            palp_ix = palp_ix + s_policy{p,s};
            palp_b = palp_b + d_func{p,s} - s_policy{p,s};
            
            if s >= 2
                palp_i = palp_i + o_policy{p,s-1};
                palp_ix = palp_ix - o_policy{p,s-1};
            end
        end
        
        for l = 1 : L
            [sig_func_s{l,p,t}, sig_var_s{l,p,t}] = polynomial(u',deg-1);
            [sig_func_o{l,p,t}, sig_var_o{l,p,t}] = polynomial(u',deg-1);
            [sig_func_i{l,p,t}, sig_var_i{l,p,t}] = polynomial(u',deg-1);
            [sig_func_ix{l,p,t}, sig_var_ix{l,p,t}] = polynomial(u',deg-1);
            [sig_func_b{l,p,t}, sig_var_b{l,p,t}] = polynomial(u',deg-1);
            
            dvar = [ dvar; sig_var_s{l,p,t}; sig_var_o{l,p,t}; sig_var_i{l,p,t}; sig_var_ix{l,p,t}; sig_var_b{l,p,t} ];
            
            palp_s = palp_s - sig_func_s{l,p,t}*cons_func(l,1);
            palp_o = palp_o - sig_func_o{l,p,t}*cons_func(l,1);
            palp_i = palp_i - sig_func_i{l,p,t}*cons_func(l,1);
            palp_ix = palp_ix - sig_func_ix{l,p,t}*cons_func(l,1);
            palp_b = palp_b - sig_func_b{l,p,t}*cons_func(l,1);
            
            F = [ F ; sos(sig_func_s{l,p,t}) ; sos(sig_func_o{l,p,t}) ; sos(sig_func_i{l,p,t}) ; sos(sig_func_ix{l,p,t}) ; sos(sig_func_b{l,p,t}) ];
        end
        F = [ F ; sos(palp_s) ; sos(palp_o) ; sos(palp_i) ; sos(palp_ix) ; sos(palp_b) ];
    end
end

diagonosis = solvesos(F,obj,myset,dvar);

ret_code = diagonosis.problem;
run_time = diagonosis.solvertime;
obj_val =  double(obj);

end

