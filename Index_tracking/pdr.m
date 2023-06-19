function [obj_val, run_time, ret_code] = pdr(data, deg)

yalmip('clear')

myset = sdpsettings;
myset.verbose                      =  1;
myset.solver                       = 'mosek';
myset.debug                        = 1;
myset.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600;



T = data.T;
G = data.G;
FF = data.F;
ff = data.f;
s0 = data.s0;

m = size(FF,1);
k = size(FF,2);
n = m - 1;

K = k*T;

ll  = sdpvar(1,1);
dvar = ll;

u0 = sdpvar(n,1);
dvar = [ dvar; u0 ];

u = sdpvar(K,1);

Wvar = cell(1,T);
Svar = cell(1,T);
Uvar = cell(n,T);

for t = 1 : T
    dim_u = nchoosek(k*t+deg-1,deg-1);
    dim_w = nchoosek(k*t+deg,deg);
    dim_s = nchoosek(k*t+deg,deg);
    
    Wvar{1,t} = sdpvar(dim_w,1);
    Svar{1,t} = sdpvar(dim_s,1);
    dvar = [dvar; Wvar{1,t}; Svar{1,t} ];
    
    for i = 1 : n
        Uvar{i,t} = sdpvar(dim_u,1);
        dvar = [ dvar; Uvar{i,t} ];
    end
end

f_func = cell(m,T);
u_policy = cell(n,T);

w_policy = cell(1,T);
s_policy = cell(1,T);

for t = 1 : T
    w_policy{1,t} = Wvar{1,t}'*monolist(u(1:k*t)',deg);
    s_policy{1,t} = Svar{1,t}'*monolist(u(1:k*t)',deg);
    
    for j = 1 : m
        f_func{j,t} = FF(j,1:k)*u((t-1)*k+1:t*k) + ff(j);
    end
    
    for i = 1 : n
        u_policy{i,t} = Uvar{i,t}'*monolist(u(1:k*t)',deg-1);
    end        
end

P = kron(eye(T),G);
p = (k-1)*ones(T*2^k,1);
          
cons_func = [ u + ones(K,1) ;
             -u + ones(K,1) ;
              P*u + p 
             ];   
L = size(cons_func,1);

obj = ll;

cons = [ ];

cons = [ cons; ones(1,n)*u0 <= s0 ; u0 >= 0 ];

palp = ll;
for t = 1 : T
    palp = palp - w_policy{1,t};
end

sig_func = cell(1,L);
sig_var  = cell(1,L);
for l = 1 : L
    [sig_func{1,l}, sig_var{1,l}] = polynomial(u',deg-1);
    
    dvar = [ dvar; sig_var{1,l}];
    
    palp = palp - sig_func{1,l}*cons_func(l,1);
    cons = [ cons ; sos(sig_func{1,l}) ];
end
cons = [ cons ; sos(palp) ];

%% For S and O policies
sig_func_wp = cell(L,T);
sig_var_wp  = cell(L,T);
sig_func_wm = cell(L,T);
sig_var_wm  = cell(L,T);
sig_func_um = cell(L,T);
sig_var_um  = cell(L,T);
sig_func_sp = cell(L,T);
sig_var_sp  = cell(L,T);
sig_func_sm = cell(L,T);
sig_var_sm  = cell(L,T);
sig_func_u = cell(L,n,T);
sig_var_u = cell(L,n,T);

%% for w_t(xi^t) >= xi_t - s_t(xi^t) and  w_t(xi^t) >= -xi_t + s_t(xi^t)
for t = 1 : T
    palp_wp = w_policy{1,t} + s_policy{1,t} - f_func{m,t};
    palp_wm = w_policy{1,t} - s_policy{1,t} + f_func{m,t};
    
    for l = 1 : L
        [sig_func_wp{l,t}, sig_var_wp{l,t}] = polynomial(u',deg-1);
        [sig_func_wm{l,t}, sig_var_wm{l,t}] = polynomial(u',deg-1);
        
        dvar = [ dvar; sig_var_wp{l,t} ; sig_var_wm{l,t} ];
        
        palp_wp = palp_wp - sig_func_wp{l,t}*cons_func(l,1);
        palp_wm = palp_wm - sig_func_wm{l,t}*cons_func(l,1);    
        
        cons = [ cons ; sos(sig_func_wp{l,t}); sos(sig_func_wm{l,t}) ];
    end
    cons = [ cons ; sos(palp_wp); sos(palp_wm) ];
end

%% for s_t(xi^t) >= xi_t'*u_{t-1}(xi^{t-1}) and s_t(xi^t) <= xi_t'*u_{t-1}(xi^{t-1})
for t = 1 : T
    palp_sp = s_policy{1,t};
    palp_sm = -s_policy{1,t};
    
    for i = 1 : n
        if t <= 1
            palp_sp = palp_sp - f_func{i,t}*u0(i);
            palp_sm = palp_sm + f_func{i,t}*u0(i); 
        else
            palp_sp = palp_sp - f_func{i,t}*u_policy{i,t-1};
            palp_sm = palp_sm + f_func{i,t}*u_policy{i,t-1};
        end
    end
    
    for l = 1 : L     
        [sig_func_sp{l,t}, sig_var_sp{l,t}] = polynomial(u',deg-1);
        [sig_func_sm{l,t}, sig_var_sm{l,t}] = polynomial(u',deg-1);
        
        dvar = [ dvar; sig_var_sp{l,t} ; sig_var_sm{l,t}  ];

        palp_sp = palp_sp - sig_func_sp{l,t}*cons_func(l,1);
        palp_sm = palp_sm - sig_func_sm{l,t}*cons_func(l,1);
        
        
        cons = [ cons ; sos(sig_func_sp{l,t}); sos(sig_func_sm{l,t}) ];
    end
    cons = [ cons ; sos(palp_sp); sos(palp_sm)];
end

%% for u_t(xi^t) >= 0
for t = 1 : T
    for i = 1 : n
        [t,i]
        palp_u = u_policy{i,t};
        for l = 1 : L
            [sig_func_u{l,i,t}, sig_var_u{l,i,t}] = polynomial(u',deg-1);
            dvar = [ dvar; sig_var_u{l,i,t} ];
            palp_u = palp_u - sig_func_u{l,i,t}*cons_func(l,1);            
            cons = [ cons ; sos(sig_func_u{l,i,t}) ];
        end
        cons = [ cons; sos(palp_u) ];
    end
end

for t = 1 : T
    palp_um = s_policy{1,t};
    
    for i = 1 : n
        palp_um = palp_um - u_policy{i,t};
    end
 
    for l = 1 : L       
        [sig_func_um{l,t}, sig_var_um{l,t}] = polynomial(u',deg-1);
        
        dvar = [ dvar; sig_var_um{l,t} ];

        palp_um = palp_um - sig_func_um{l,t}*cons_func(l,1);
        cons = [ cons ; sos(sig_func_um{l,t}) ];
    end
    cons = [ cons ; sos(palp_um) ];
end

diagnostics = solvesos(cons,obj,myset,dvar);

ret_code = diagnostics.problem;
run_time = diagnostics.solvertime;
obj_val = double(obj);

end

