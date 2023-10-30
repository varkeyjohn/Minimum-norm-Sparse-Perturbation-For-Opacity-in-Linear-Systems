function rk_struct_val = struct_pert_fun(s_val,theta,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind)
    n = size(A,1);
    Mx = [A-s_val*exp(1i*theta)*eye(n),B;C,D];
    Rp = Mx(Rp_ind,:);
    Rpc = Mx(Rpc_ind,:);
    [Q,R] = qr(vpa(Rpc'));
    Q2 = Q(:,size(Rpc,1)+1:size(Q,2));
    Q2_beta = Q2(Cp_ind,:);
    M = vpa(Rp*Q2);
    N = vpa(Q2_beta);
    i = min(size(M,1),size(M,2));
    rk_max = -inf;
    for gamma_k = 1e-3:1e-1:1
        P1 = [real(M),-gamma_k*imag(M);inv(gamma_k)*imag(M),real(M)];
        P2 = [real(N),-gamma_k*imag(N);inv(gamma_k)*imag(N),real(N)];
        rk_temp = flip(vpa(gsvd(double(P1),double(P2))));
        rk_temp = rk_temp(2*i-1);
        if(rk_temp>rk_max)
            rk_max = rk_temp;
            gamma_val = gamma_k;
        end
    end
    %gamma_val
    rk_struct_val = rk_max;
end