function svals = svals_theorem_approx(gammak,A,B,C,D,E_pert,G_pert,epsilon,rk,theta)
    n = size(A,1);
    m = size(C,1);
    p = size(B,2);
    svd_E_pert = svd(E_pert);
    svd_E_pert(find(~svd_E_pert))=vpa(epsilon*ones(length(find(~svd_E_pert)),1));
    [UE_pert,SE_pert,VE_pert] = svd(E_pert);
    E_pert_approx = vpa(vpa(UE_pert)*diag(svd_E_pert)*vpa(VE_pert'));
    svd_G_pert = svd(G_pert);
    svd_G_pert(find(~svd_G_pert))=vpa(epsilon*ones(length(find(~svd_G_pert)),1));
    [UG_pert,SG_pert,VG_pert] = svd(G_pert);
    G_pert_approx = vpa(vpa(UG_pert)*diag(svd_G_pert)*vpa(VG_pert'));
    
    Pl = [eye(2*n),zeros(2*n,2*m),zeros(2*n,2*n),zeros(2*n,2*p);...
          zeros(2*n,2*n),zeros(2*n,2*m),eye(2*n),zeros(2*n,2*p);...
          zeros(2*m,2*n),eye(2*m),zeros(2*m,2*n),zeros(2*m,2*p);...
          zeros(2*p,2*n),zeros(2*p,2*m),zeros(2*p,2*n),eye(2*p)];%verified
    Pr = [zeros(2*n),eye(2*n),zeros(2*n,2*p),zeros(2*n,2*m);...
          zeros(2*m,2*n),zeros(2*m,2*n),zeros(2*m,2*p),eye(2*m);...
          eye(2*n),zeros(2*n),zeros(2*n,2*p),zeros(2*n,2*m);...
          zeros(2*p,2*n),zeros(2*p,2*n),eye(2*p),zeros(2*p,2*m)]; %verified % these are 2n, 2m etc. because they follow dimensions of A,B,C,D in H matrix, which are now blockdiagonal
    

    Tgnr = (1/(2*gammak))*([sqrt(1+gammak^2)*eye(n+m),1i*gammak*sqrt(1+gammak^2)*eye(n+m);...
                            sqrt(1+gammak^2)*eye(n+m),-1i*gammak*sqrt(1+gammak^2)*eye(n+m)]);%verified
    Tgnm = (1/(2*gammak))*([sqrt(1+gammak^2)*eye(n+p),1i*gammak*sqrt(1+gammak^2)*eye(n+p);...
                            sqrt(1+gammak^2)*eye(n+p),-1i*gammak*sqrt(1+gammak^2)*eye(n+p)]);%verified
    Pnr = [eye(n),zeros(n,n),zeros(n,m),zeros(n,m);...
           zeros(m,n),zeros(m,n),eye(m),zeros(m,m);...
           zeros(n,n),eye(n),zeros(n,m),zeros(n,m);...
           zeros(m,n),zeros(m,n),zeros(m,m),eye(m)];%verified
    Pnm = [eye(n),zeros(n,p),zeros(n,n),zeros(n,p);...
           zeros(n,n),zeros(n,p),eye(n),zeros(n,p);...
           zeros(p,n),eye(p),zeros(p,n),zeros(p,p);...
           zeros(p,n),zeros(p,p),zeros(p,n),eye(p)];%verified
    E_pert_approxinvcap = blkdiag(vpa(inv(E_pert_approx)),vpa(inv(E_pert_approx)));
    G_pert_approxinvcap = blkdiag(vpa(inv(G_pert_approx)),vpa(inv(G_pert_approx)));
    E = vpa(inv(Tgnr))*E_pert_approxinvcap*Pnr;%verified %added E_pert
    F = Pnm*G_pert_approxinvcap*Tgnm;%verified %added G_pert
    T = Pl*[-rk*inv(E'*E),zeros(2*(n+m),2*(n+p));...
            zeros(2*(n+p),2*(n+m)),-rk*inv(F*F')]*Pr;%verified - since E is 2(n+m)*2(n+m), F is 2(n+p)*2(n+p).
    Acap = blkdiag(A,A);%verified
    Bcap = blkdiag(B,B);%verified
    Ccap = blkdiag(C,C);%verified
    Dcap = blkdiag(D,D);%verified

    Atilde = [Acap,zeros(2*n);...
              zeros(2*n),Acap.']+T(1:2*2*n,1:2*2*n);%verified
    Btilde = [Bcap,zeros(2*n,2*m);...
              zeros(2*n,2*p),Ccap.']+T(1:2*2*n,2*2*n+1:size(T,2));%verified
    Ctilde = [Ccap,zeros(2*m,2*n);...
              zeros(2*p,2*n),Bcap.']+T(2*2*n+1:size(T,1),1:2*2*n);%verified
    Dtilde = [Dcap,zeros(2*m);...
              zeros(2*p),Dcap.']+T(2*2*n+1:size(T,1),2*2*n+1:size(T,2));%verified
    Deltacap = blkdiag(exp(1i*theta)*eye(n),exp(-1i*theta)*eye(n));%verified
    Deltatilde = [Deltacap,zeros(2*n);...
                  zeros(2*n),Deltacap'];%verified
    if(rank(Dtilde)<min(size(Dtilde,1),size(Dtilde,2)))
        %digits(1000);
        svals = eig(double([Atilde,Btilde;Ctilde,Dtilde]),double([Deltatilde,zeros(size(Btilde));zeros(size(Ctilde)),zeros(size(Dtilde))]));
        %disp("rank less")
        svals = svals(isfinite(svals)); % removing inf values
        %digits(32)
    else
        H = vpa(inv(Deltatilde))*(Atilde-Btilde*vpa(inv(Dtilde))*Ctilde);%verified
        svals = eig(H);%verified
        %disp("rank prop")
    end
end