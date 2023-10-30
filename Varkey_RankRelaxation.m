clc;clear;
%% Matrix initialization %%
Ao = [0.74,-0.69,-2.08;-0.12,1.62,0.63;-0.38,-0.21,0.14]; %nxn
Bo = [-1.23,-0.26;1.02,2.51;-0.66,1.13]; %nxm
Co = [1.06,0.71,0.61]; %rxn
Do = [1.33,-2.89]; %rxm

A = Ao';
B = Co';
C = Bo';
D = Do';
n = size(A,1);
m = size(B,2);
r = size(C,1);
epsilon = 1e-4;
Zeta = 1e-8;

%% Initialization k=0- %% %rand(20,1);
theta_tilde0 = [-0.0270,-0.2079];%[0.034109964669723415843918383980647,-0.068234818039511811478736810602993,0.20475353947158144942159090561862,0.030727269863612368401629180808088];%rand(20,1);%[3.72e-3,5.2e-3,4.49e-2,2.18e-2,-5.49e-3,-7.67e-3,-6.63e-2,-3.21e-2,3.28e-4,4.58e-4,3.95e-3,1.98e-2,-6.38e-2,-4.2e-2,1.98e-2,-2.2e-2,-2.86e-2,-1.87e-2,1.1e-2,-8.89e-3];
%theta_tilde0 = %[0.0090, 0.0081, 0.0334, 0.0192;...
     %0.0001,-0.0044,-0.0747,-0.0331;...
     %0.0222, 0.0153, 0.0047, 0.0128;...
    %-0.0641,-0.0420, 0.0214,-0.0216;...
    %-0.0208,-0.0137, 0.0056,-0.0076];
Adelta_tilde0 = [theta_tilde0(1),0,theta_tilde0(2);...
                 0,0,0;...
                 0,0,0];
Bdelta_tilde0 = zeros(3,2);
Cdelta_tilde0 = zeros(1,3);
Ddelta_tilde0 = zeros(1,2);
%Adelta_tilde0 = [theta_tilde0(1),theta_tilde0(5),theta_tilde0(9);...
%          theta_tilde0(2),theta_tilde0(6),theta_tilde0(10);...
%          theta_tilde0(3),theta_tilde0(7),theta_tilde0(11)];
%Bdelta_tilde0 = [theta_tilde0(13),theta_tilde0(17);...
%           theta_tilde0(14),theta_tilde0(18);...
%           theta_tilde0(15),theta_tilde0(19)];
%Cdelta_tilde0 = [theta_tilde0(4),theta_tilde0(8),theta_tilde0(12)];
%Ddelta_tilde0 = [theta_tilde0(16),theta_tilde0(20)];
Atheta_tilde0 = A + Adelta_tilde0';
Btheta_tilde0 = B + Cdelta_tilde0';
Ctheta_tilde0 = C + Bdelta_tilde0';
Dtheta_tilde0 = D + Ddelta_tilde0';

mu_tilde0 = 0.8108;%-2 + (2 + 2) * rand(1);%0.8209;
lambda_tilde0 = 0.5367;%-2 + (2 + 2) * rand(1);%0.2329;
Z_tilde0 = [Atheta_tilde0-mu_tilde0*eye(n),Btheta_tilde0,-lambda_tilde0*eye(n),zeros(n,m);...
            Ctheta_tilde0,Dtheta_tilde0,zeros(r,n),zeros(r,m);...
            lambda_tilde0*eye(n),zeros(n,m),Atheta_tilde0-mu_tilde0*eye(n),Btheta_tilde0;...
            zeros(r,n),zeros(r,m),Ctheta_tilde0,Dtheta_tilde0];
[U_tilde,S_tilde,V_tilde] = svd(Z_tilde0);
U_tilde1 = U_tilde(:,1:2*(n+min(r,m))-1);
V_tildet = V_tilde';
V_tilde1 = V_tildet(1:2*(n+min(r,m))-1,:)';
%norm(theta_tilde0)
norm([Adelta_tilde0,Bdelta_tilde0;Cdelta_tilde0,Ddelta_tilde0]')
%% Initialization k=0 %%
% Solve convex problem to obtain theta0, lambda0, mu0
cvx_begin quiet
    variable theta0(3)
    variable lambda0(1)
    variable mu0(1)
    variable Z0(2*(n+r),2*(n+m))
    Adelta0 = [theta0(1),0,theta0(3);...
                     0,0,0;...
                     theta0(2),0,0];
    Bdelta0 = zeros(3,2);
    Cdelta0 = zeros(1,3);
    Ddelta0 = zeros(1,2);
    %{
    Adelta0 = [theta0(1),theta0(2),theta0(3);...
               theta0(4),theta0(5),theta0(6);...
               theta0(7),theta0(8),theta0(9)];
    Bdelta0 = [theta0(10),theta0(11);...
               theta0(12),theta0(13);...
               theta0(14),theta0(15)];
    Cdelta0 = [theta0(16),theta0(17),theta0(18)];
    Ddelta0 = [theta0(19),theta0(20)];
    %}
    Atheta0 = A+Adelta0';
    Btheta0 = B+Cdelta0';
    Ctheta0 = C+Bdelta0';
    Dtheta0 = D+Ddelta0';
    minimize(norm_nuc(Z0)-trace(U_tilde1'*Z0*V_tilde1))
    subject to
        Z0 == [Atheta0-mu0*eye(n),Btheta0,-lambda0*eye(n),zeros(n,m);...
                   Ctheta0,Dtheta0,zeros(r,n),zeros(r,m);...
                   lambda0*eye(n),zeros(n,m),Atheta0-mu0*eye(n),Btheta0;...
                   zeros(r,n),zeros(r,m),Ctheta0,Dtheta0]; %#ok<EQEFF>
cvx_end
Zk = Z0;
gtheta0=norm([Adelta_tilde0,Bdelta_tilde0;Cdelta_tilde0,Ddelta_tilde0]');
%gtheta0 = norm([Adelta_tilde0,Bdelta_tilde0;Cdelta_tilde0,Ddelta_tilde0]);
gamma = gtheta0/epsilon;% min(5,gtheta0/epsilon);
gthetak = gtheta0;
Fk = gthetak+gamma*(nuclear_norm(Zk,2*(n+min(r,m)))-r_norm(Zk,2*(n+min(r,m))-1)); % since F(Z0) = g(theta0) according to paper.
%Fk
%% Iterative Algorithm %%
while(1)
    % Step 3: U1(k), V1(k) from SVD of Z(k)
    [Uk,Sk,Vk] = svds(Zk,2*(n+min(r,m))-1);
    Uk1 = Uk(:,1:2*(n+min(r,m))-1);
    Vkt = Vk';
    Vk1 = Vkt(1:2*(n+min(r,m))-1,:)';
    Fkminus = Fk;
    % Step 4: a. Set gamma (Set already to prevent reinitialization in loop)
    %         b. Solve convex problem to obtain theta(k+1), lambda(k+1), mu(k+1)
    cvx_begin quiet
        variable thetak(3)
        variable lambdak(1)
        variable muk(1)
        variable Zk(2*(n+r),2*(n+m))
        Adeltak = [thetak(1),0,thetak(3);...
                         0,0,0;...
                         thetak(2),0,0];
        Bdeltak = zeros(3,2);
        Cdeltak = zeros(1,3);
        Ddeltak = zeros(1,2);
        %{
        Adeltak = [thetak(1),thetak(5),thetak(9);...
                   thetak(2),thetak(6),thetak(10);...
                   thetak(3),thetak(7),thetak(11)];
        Bdeltak = [thetak(13),thetak(17);...
                   thetak(14),thetak(18);...
                   thetak(15),thetak(19)];
        Cdeltak = [thetak(4),thetak(8),thetak(12)];
        Ddeltak = [thetak(16),thetak(20)];
        %}
        Athetak = A+Adeltak';
        Bthetak = B+Cdeltak';
        Cthetak = C+Bdeltak';
        Dthetak = D+Ddeltak';
        minimize(norm([Adeltak,Bdeltak;Cdeltak,Ddeltak]')+gamma*(norm_nuc(Zk)-trace(Uk1'*Zk*Vk1)))
        subject to
            Zk == [Athetak-muk*eye(n),Bthetak,-lambdak*eye(n),zeros(n,m);...
                   Cthetak,Dthetak,zeros(r,n),zeros(r,m);...
                   lambdak*eye(n),zeros(n,m),Athetak-muk*eye(n),Bthetak;...
                   zeros(r,n),zeros(r,m),Cthetak,Dthetak]; %#ok<EQEFF>
    cvx_end
    Adeltak = [thetak(1),0,thetak(3);...
                 0,0,0;...
                 thetak(2),0,0];
    Bdeltak = zeros(3,2);
    Cdeltak = zeros(1,3);
    Ddeltak = zeros(1,2);
    %{
    Adeltak = [thetak(1),thetak(5),thetak(9);...
           thetak(2),thetak(6),thetak(10);...
           thetak(3),thetak(7),thetak(11)];
    Bdeltak = [thetak(13),thetak(17);...
               thetak(14),thetak(18);...
               thetak(15),thetak(19)];
    Cdeltak = [thetak(4),thetak(8),thetak(12)];
    Ddeltak = [thetak(16),thetak(20)];
    %}
    gthetak = norm([Adeltak,Bdeltak;Cdeltak,Ddeltak]');
    %gthetak = norm([Adelta_tilde0,Bdelta_tilde0;Cdelta_tilde0,Ddelta_tilde0]);
    Fk = gthetak+gamma*(nuclear_norm(Zk,2*(n+min(r,m)))-r_norm(Zk,2*(n+min(r,m))-1));
    vpa(Fk)
    %epsilon = 1e7*(nuclear_norm(Zk)-r_norm(Zk,2*n-1));
    %if(epsilon < epsilon_prev)
    %    gamma = gtheta0/epsilon;
    %    epsilon_prev = epsilon;
    %end
    %min(svd([Athetak-(muk+lambdak*1i)*eye(n),Bthetak]))
    %epsilon_prev
    %gamma
    % Step 5: Update k (Do nothing)
    % End on convergence
    if(abs(Fk - Fkminus) <= Zeta)
        break;
    end
    if(Fk > Fkminus)
        Fk - Fkminus
        disp("Fk > Fkmius!");
    end
    %min(svd([Athetak-(muk+lambdak*1i)*eye(n),Bthetak;Cthetak,Dthetak]))
end
disp("theta")
disp(thetak);
disp("norm(theta)")
Adeltak = [thetak(1),0,thetak(3);...
             0,0,0;...
             thetak(2),0,0];
Bdeltak = zeros(3,2);
Cdeltak = zeros(1,3);
Ddeltak = zeros(1,2);
%{
Adeltak = [thetak(1),thetak(5),thetak(9);...
       thetak(2),thetak(6),thetak(10);...
       thetak(3),thetak(7),thetak(11)];
Bdeltak = [thetak(13),thetak(17);...
           thetak(14),thetak(18);...
           thetak(15),thetak(19)];
Cdeltak = [thetak(4),thetak(8),thetak(12)];
Ddeltak = [thetak(16),thetak(20)];
%}
%disp(norm(thetak));
disp(vpa(norm([Adeltak,Bdeltak;Cdeltak,Ddeltak]')));
%rank(Zk)
min(svd([Athetak-(muk+lambdak*1i)*eye(n),Bthetak;Cthetak,Dthetak]))