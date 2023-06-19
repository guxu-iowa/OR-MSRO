clear
clc

rslt = zeros(4,100);
tm   = zeros(4,100);

for seed = 1 : 100
    seed
    rng(seed);
    
    % Set base dimension
    n = 5;
    
    p = 60*ones(n,1);
    r = 80*ones(n,1);
    
    s = zeros(n,1);
    c = 0.5*r + 10*rand(n,1);
    Gamma = 4;
    
    dbar = 60*ones(n,1);
    dhat = 55*ones(n,1) + 5*(2*rand(n,1)-1);

    
    I = eye(n);
    
    M = (2*rand(n,n)-1)/n;
    
    d = -ones(n,1);
    A = zeros(2*n,n);
    B = zeros(2*n,n);
    F = zeros(2*n,2*n);
    f = zeros(2*n,1);
    
    PP = zeros(n,2*n);
    for nidx = 1 : n
        A(2*nidx-1,nidx) = s(nidx) - c(nidx);
        A(2*nidx,nidx) = r(nidx) - c(nidx) + p(nidx);
        B(2*nidx-1,:) = -I(nidx,:);
        B(2*nidx,:) = -I(nidx,:);
        
        F(2*nidx-1,:) = [(s(nidx)-r(nidx))*dhat(nidx)*M(nidx,:), - (s(nidx)-r(nidx))*dhat(nidx)*M(nidx,:)];
        F(2*nidx,:) = [p(nidx)*dhat(nidx)*M(nidx,:), -p(nidx)*dhat(nidx)*M(nidx,:)];
        
        PP(nidx,:) = [ -I(nidx,:), -I(nidx,:) ];
        
        f(2*nidx-1,1) = (s(nidx) - r(nidx))*dbar(nidx);
        f(2*nidx,1) = p(nidx)*dbar(nidx);
    end

  
  F = [f, F];
  
  p = [ Gamma; ones(n,1); zeros(2*n,1) ];
  
  PPP = [ -ones(1, 2*n) ;PP; eye(2*n) ];
  
  P = [p, PPP];
  
  [obj_qdr, x_qdr, tm_qdr] = newsvendor_qdr(A, B, d, F, P);
  tm(1,seed) = tm_qdr;
  
  
  [obj_cop,x_cop, tm_cop] = newsvendor_cop(A, B, d, F, P);
  tm(2,seed) = tm_cop;
  
  
  [obj_as, x_as, tm_as] = newsvendor_as(A, B, d, F, P);
  tm(3,seed) = tm_as;
  
  
  p = [ Gamma; ones(n,1); zeros(2*n,1) ];
  
  P = [ -ones(1, 2*n) ; PP; eye(2*n) ];
  
  DEG = 3;
  F = F(:,2:end);
  [obj_poly_3, x_poly_3, tm_poly_3] = newsvendor_poly(A, B, d, F, f, P, p, DEG);
  tm(4,seed) = tm_poly_3;
  

  %
  rslt(1,seed) = obj_qdr;
  rslt(2,seed) = obj_cop;
  rslt(3,seed) = obj_as;
  rslt(4,seed) = obj_poly_3;
  
  save('result.mat', 'rslt', 'tm');
  
end


 
 
 
 
 
