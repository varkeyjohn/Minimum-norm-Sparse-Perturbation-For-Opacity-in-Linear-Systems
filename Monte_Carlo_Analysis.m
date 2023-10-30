clc;
clear all;
Nrun = 1e6;
count = 0;
for run = 1:Nrun
    theta = rand()*2*pi;
    rkminus = 1e230;%(b-a)*rand() + a;
    epsilon = 1e-230;
    digits(200);
    E_pert = diag([1,0,0,0,0]);
    G_pert = diag([0,1,0,0,0]);
    A = 1e3*rand(3,3);
    B = 1e3*rand(3,2);
    C = 1e3*rand(2,3);
    D = 1e3*rand(2,2);
    s = svals_theorem_approx(1e-3,A,B,C,D,E_pert,G_pert,epsilon,rkminus,theta);
    s = sort(s,'descend');
    flag = 0;
    for l = 1:length(s)
        if(abs(double(((struct_pert_fun(s(l),theta,A,B,C,D,[1],[2,3,4,5],2)-rkminus)/rkminus)*100))==0)
            abs(double(((struct_pert_fun(s(l),theta,A,B,C,D,[1],[2,3,4,5],2)-rkminus)/rkminus)*100));
            count = count+1
            flag = 1;
            break;
        end
    end
    if flag==0
        disp("phi")
    end
end
Percent = count/Nrun*100;
disp(Percent);
