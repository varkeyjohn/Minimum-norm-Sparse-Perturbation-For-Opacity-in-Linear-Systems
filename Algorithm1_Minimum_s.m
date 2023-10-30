clc;
clear all;
c0 = clock;
Ao = [0.74,-0.69,-2.08;-0.12,1.62,0.63;-0.38,-0.21,0.14]; %nxn
Bo = [-1.23,-0.26;1.02,2.51;-0.66,1.13]; %nxm
Co = [1.06,0.71,0.61]; %rxn
Do = [1.33,-2.89]; %rxm

A = Ao';
B = Co';
C = Bo';
D = Do';

n = size(A,1);
m = size(C,1);
p = size(D,2);

E_pert = [1,0,0,0,0;...
          0,0,0,0,0;...
          0,0,1,0,0;...
          0,0,0,0,0;...
          0,0,0,0,0];

G_pert = [1,0,0,0;...
          0,0,0,0;...
          0,0,1,0;...
          0,0,0,0];
TOL_Sk = 5;
TOL_Rk = 1;
TOL_rk = 0.001;%0.001

%% Step 1.1 %%
S0 = 0:0.01:2*pi;%0.008 %Step 1a
S = S0; % Step 2a

%% Step 1.2 %%
s0 = 0; %Step 1b(i)
s = s0; %s and s0 are only to find r0, gammak. Later s will be changed to vector of eigenvalues of H % Step 2a

Rp_ind = [1,3]; %Rows to perturb
Rpc_ind = setdiff([1:n+m],Rp_ind);%Rows to not perturb
Cp_ind = [1,3];%Columns to perturb
Cpc_ind = setdiff([1:n+p],Cp_ind);%Columns to not perturb
Delta = zeros(n+m,n+p);
Mx = [A-s*eye(n),B;C,D];
Rp = Mx(Rp_ind,:);
Rpc = Mx(Rpc_ind,:);
[Q,R] = qr(Rpc');
Q2 = Q(:,size(Rpc,1)+1:size(Q,2));
%Q2 = Q(:,n+1:size(Q,2));
Q2_beta = Q2(Cp_ind,:);
epsilon = 1e-20;
[r0,gamma0] = min_struct_pert_fun(abs(s0),angle(s0),A,B,C,D,Rp_ind,Rpc_ind,Cp_ind);
rk = r0; % Step 2a
gammak = gamma0; % Step 2a
R = []; 
Rtheta = [];
%disp("Line 57")
%% Step 2 %%
while(1)
    %% Step 2-a (Given) %%
    % r_{k-1}=rk, gamma_{k-1}=gammak, s_{k-1}=sk, S_{k-1}=Sk
    %% Step 2-b %%
    Theta = []; % Step 2b
    R = [];
    theta_arr = [];
    %flag = 0;
    
    %% Step 2-c %%
    for theta_num = 1:length(S) % Step 2c
        % Step 2c(i)_1 -> Computing H and s\in eig(H)
        flag = 0;
        theta = S(theta_num);
        %disp(theta)
        %disp("Line 70")
        %disp(theta)
        s = svals_theorem_approx(gammak,A,B,C,D,E_pert,G_pert,epsilon,rk,theta);
        Rtheta = [];
        for j = 1:length(s)
            [sigma,gammak_temp] = min_struct_pert_fun(s(j),theta,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind);
            if (abs(sigma - rk)<1e-6 && real(s(j))>=0 && abs(imag(s(j)))<1e-5) % is sigma=rk and s(j) is real, non-negative
                Rtheta = union(Rtheta,real(s(j)));
            end
        end
        % Step 2c(ii)
        if(~isempty(Rtheta))
            %Theta = union(Theta,theta);
        else
            % Step 2c - if there are no real non-negative eigenvalues, then
            % either Rtheta=phi or Rtheta is R+. This is checked by taking
            % a testpoint, say=1. If R+, then we get something, otherwise
            % Rtheta=phi and we move on.
            [sigma,gammak_temp] = min_struct_pert_fun(1,theta,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind);
            if(sigma < rk) %Rtheta is R+, not phi (Based on definition of Rktheta)
                %disp("Line 91")
                if(~isempty(intersect(S,theta)))
                    flag = 1;
                    Rthetanext=1;%arbitrary w=1
                    thetanext=theta;
                    rk_max = sigma;
                    rk = rk_max;
                    gammak = gammak_temp;
                end
            end
        end
        
        if(isempty(Rtheta) && flag==0)
            S(theta_num)=[];
        end
        R = union(R,Rtheta);
        if(flag == 0 && ~isempty(Rtheta))
            for l = 1:length(Rtheta)+1
                if l == 1
                    min_val = 0;
                else
                    min_val = Rtheta(l-1);
                end
                if l == length(Rtheta)+1
                    if(~isempty(Rtheta) && ~isempty(intersect(S,theta)))
                        Rthetanext_temp = Rtheta(l-1)+1;% Rtheta(l-1)+1 should be minimum since others werent
                    end
                else
                    max_val = Rtheta(l);
                    Rthetanext_temp = (max_val+min_val)/2;%(max_val-min_val)*rand(1) + min_val;
                end
                thetanext_temp = theta;
                [rk_max,gammak_temp] = min_struct_pert_fun(Rthetanext_temp,thetanext_temp,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind);
                if(rk_max<rk && ~isempty(intersect(S,theta)))
                    %disp("Line 122")
                    flag = 1;
                    Rthetanext=Rthetanext_temp;
                    thetanext=thetanext_temp;
                    rk = rk_max;
                    gammak = gammak_temp;
                    break;
                end
            end
        end
        if(flag == 1)
            break;
        end
        if(theta_num >= length(S))
            break;
        end
    end
    %% Step 2-d %%
    %S = intersect(S,Theta);
    %% Step 2-e %%
    if(size(S,2) < TOL_Sk && size(R,1) < TOL_Rk)
        break;
    end
    
    %% Step 2-f %%
    %[rk,gammak] = min_struct_pert_fun(Rthetanext,thetanext,A,B,C,D,Rp_ind,Rpc_ind,Cp_ind);
    disp(vpa(rk));
    c1 = clock;
    %% Step 3 %%
    if(rk < TOL_rk)
        break;
    end
    if(length(S) <= TOL_Sk)
        c1 = clock;
        break;
    end
    %% Step 4 %%
    % k is auto-updated
end
