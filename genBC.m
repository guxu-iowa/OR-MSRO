function [ Balpha, Calpha ] = genBC( dd, p, P )

s = size(P,1);
k = size(P,2);

dim_palpha = nchoosek(k+2*dd,2*dd);
dim_x = nchoosek(k+dd,dd);

dz = dd-1; % for polyhedral uncertainty set Pu + p >= 0
dim_z = nchoosek(k+dz,dz);

[bin,pow,bas] = genind(k,2*dd);


Balpha = zeros(dim_x,dim_x,dim_palpha);
Calpha = zeros(dim_z,dim_z,s,dim_palpha);

for i = 1 : dim_palpha
    if mod(i,1000) == 0
        i
    end
    vv = pow(i,:);
    
    % process Balp
    T = zeros(dim_x,dim_x,k);
    for j = 1 : k
        T(:,:,j) = (bas(:,:,j)==vv(j));
    end

    TT = ones(dim_x,dim_x);
    for j = 1 : dim_x
        for jj = j : dim_x
            for jjj = 1 : k
                TT(j,jj) = TT(j,jj)*T(j,jj,jjj);
            end
            TT(jj,j) = TT(j,jj);
        end
    end
    Balpha(:,:,i) = TT;
    
    % process matrix Calp
    for j = 1 : s
        for l = 1 : k
            T = zeros(dim_z,dim_z,k);
            bas_new = bas(1:dim_z,1:dim_z,:);
            bas_new(:,:,l) = bas_new(:,:,l) + 1;
            for kk = 1 : k
                T(:,:,kk) = (bas_new(:,:,kk)==vv(kk)); 
            end
            TT = ones(dim_z,dim_z);
            for f_idx = 1 : dim_z
                for s_idx = f_idx : dim_z
                    for t_idx = 1 : k
                        TT(f_idx,s_idx) = TT(f_idx,s_idx)*T(f_idx,s_idx,t_idx);
                    end
                TT(s_idx,f_idx) = TT(f_idx,s_idx);
                end
            end
            Calpha(:,:,j,i) = Calpha(:,:,j,i) + P(j,l)*TT;
        end
        Calpha(:,:,j,i) = Calpha(:,:,j,i) + p(j)*Balpha(1:dim_z,1:dim_z,i); 
    end
end




end

