function [rk_struct,gammak_temp] = min_struct_pert_fun(s_val,theta,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind)%change argument to [rk_struct,gammak_temp], arg as s_val, theta
    n = size(A,1);
    Mx = [A-s_val*exp(1i*theta)*eye(n),B;C,D];
    Rp = Mx(Rp_ind,:);
    Rpc = Mx(Rpc_ind,:);
    [Q,R] = qr(Rpc');
    Q2 = Q(:,size(Rpc,1)+1:size(Q,2));
    Q2_beta = Q2(Cp_ind,:);
    M = Rp*Q2;
    N = Q2_beta;
    %G_temp = [A-Rthetanext_temp*exp(1i*thetanext_temp)*eye(n),B]*Z;
    rk_max = -inf;
    i = min(size(M,1),size(M,2));
    for gamma_i = 1e-5:1e-2:1
        P1 = [real(M),-gamma_i*imag(M);inv(gamma_i)*imag(M),real(M)];
        P2 = [real(N),-gamma_i*imag(N);inv(gamma_i)*imag(N),real(N)];
        rk_temp = flip(gsvd(P1,P2));
        rk_temp = rk_temp(2*i-1);
        if(rk_temp>rk_max)
            rk_max = rk_temp;
            gammak_temp = gamma_i;
        end
    end
    rk_struct = rk_max;
end