clc;
clear all;
DeltaAo = [4.87e-3,4.626e-4;...
          -4.869e-2,-4.67e-3];
DeltaBo = [-2.289e-4,-2.464e-3;...
           2.33e-3,4.162e-4];
DeltaCo = [-6.427e-4,-1.53e-4];
DeltaDo = [1.155e-4,-4.916e-2];

Ao = [0.74,-0.69,-2.08;-0.12,1.62,0.63;-0.38,-0.21,0.14]; %nxn
Bo = [-1.23,-0.26;1.02,2.51;-0.66,1.13]; %nxm
Co = [1.06,0.71,0.61]; %rxn
Do = [1.33,-2.89]; %rxm
A = Ao';
B = Co';
C = Bo';
D = Do';

A = 1e4;
B = [0,0];
C = [1e-4;0];
D = [1e4,1e4;0,1e-4];
n = size(A,1);
m = size(C,1);
p = size(D,2);
%[CDt,CDtid] = licols([C,D]');
%CDt';
Rp_ind = [1,2]; %Rows to perturb
Rpc_ind = setdiff([1:n+m],Rp_ind);%Rows to not perturb
Cp_ind = [1,2];%Columns to perturb
Cpc_ind = setdiff([1:n+p],Cp_ind);%Columns to not perturb
Delta = zeros(n+m,n+p);
s = 18i;
Mx = [A-s*eye(n),B;C,D];
if(rank(Mx)<min(n+m,n+p))
    disp("Already low rank");
end

Rp = Mx(Rp_ind,:);
Rpc = Mx(Rpc_ind,:);
[Q,R] = qr(Rpc'); % Remember to remove the linear dependent rows here!
Q2 = Q(:,size(Rpc,1)+1:size(Q,2));
Q2_beta = Q2(Cp_ind,:);

M = Rp*Q2;
N = Q2_beta;

svd_max = -inf;
for gamma = 1e-5:1e-4:1
    P1 = [real(M),-gamma*imag(M);inv(gamma)*imag(M),real(M)];
    P2 = [real(N),-gamma*imag(N);inv(gamma)*imag(N),real(N)];
    svd_val = flip(gsvd(P1,P2));%for decreasing ordering
    svd_val = svd_val(2*(min(size(M,1),size(M,2)))-1);
    if(svd_val>svd_max)
        svd_max = svd_val;
    end
end
t = svd_max;

H = t^2*(N'*N)-M'*M;
S = t^2*(N.'*N)-M.'*M;
Cp = inv(S)*conj(H);
[V_c,Mu_c]=eig(Cp*conj(Cp));
mu_c = eig(Cp*conj(Cp)); 
ind_u = find(abs(imag(mu_c))<(10^(-20)));
ind_v = find(imag(mu_c)>10^(-20)); 
Qu = V_c(:,ind_u);
Qv1 = V_c(:,ind_v);
Qv2 = Cp*conj(V_c(:,ind_v));
Qv = [zeros(size(Qv2)),zeros(size(Qv1))];
Qv(:,1:2:end-1) = Qv1;
Qv(:,2:2:end) = Qv2;
Qp = [Qu,Qv];
U = mat2cell(Qu,size(Qu,1),ones(size(Qu,2),1)');
U = blkdiag(U{:});
V = mat2cell(Qv1,size(Qv1,1),ones(size(Qv1,2),1)');
V = blkdiag(V{:});

alpha_inv = diag(vpa(inv(sqrt(U.'*kron(eye(size(Qu,2)),S)*U))));% S in blkdiag form. No. of blocks=no. of vectors in Qu
beta_inv =  diag(inv(sqrt(V.'*kron(eye(size(Qv1,2)),S)*V)));
gamma = sqrt(mu_c(ind_v));
gb_inv = diag(conj((inv(diag(gamma)))*diag(beta_inv)));
P2 = [beta_inv.';gb_inv.'];
P2 = P2(:);
P = Qp*diag([alpha_inv;P2]);
X = [];%a1...al,b1...bs,c1...cw This is a vector and not a matrix if reducing by one rank.
Delta = -[real(M*X),imag(M*X)]*vpa(pinv([real(N*X),imag(N*X)]),500)

