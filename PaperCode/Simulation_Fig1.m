close all; clear; clc

theta_true = [-10, 10, 37]./180*pi; 
k = length(theta_true);
P = diag([7 3 1]);    
sigma2 = 1;
rho = diag(P)'/sigma2;  
CNT = 100;
coeff = 1:2:17;
err1 = zeros(CNT,length(coeff));
err2 = zeros(CNT,length(coeff));
err3 = zeros(CNT,length(coeff));

for it = 1:1:length(coeff)
    N = 128*coeff(it);
    T  = 4*N;
    c = N/T;  % c=0.25, sqrt(c) = 0.5

    clear i
    a = @(theta) exp(pi*1i*sin(theta)*(0:N-1)')/sqrt(N);
    A = [];
    for tmp_index=1:length(theta_true)
        A = [A a(theta_true(tmp_index))];
    end
    
    n = round(N/2);
    J = eye(N);
    delta = 1;
    J1 = J(1:n,:);
    J2 = J(1+delta:n+delta,:);   % delta = 1
    
    temp1 = n/N*(1-c*rho.^(-2))./(1+c*rho.^(-1)) + n/T*(1+rho.^(-1))./(c+rho);
    temp2 = n/N*(1-c*rho.^(-2))./(1+c*rho.^(-1)).*exp(1i*pi*delta*sin(theta_true));
    Phi1_bar = diag(temp1);
    Phi2_bar = diag(temp2);
    Phi_bar = inv(Phi1_bar)*Phi2_bar;
    
    for cnt = 1:1:CNT
        S = sqrtm(P)*randn(k,T);    % K*T 
        Z = complex(randn(N,T), randn(N,T));
        X = A*S + sqrt(sigma2/2)*Z;     % N*T
        SCM = X*(X')/T;
        [U,eigs_SCM] = eig(SCM,'vector');
        [eigs_SCM, index] = sort(eigs_SCM,'descend');
        U = U(:, index);
        U_S = U(:,1:k);

        Phi1 = U_S'*J1'*J1*U_S;
        Phi2 = U_S'*J1'*J2*U_S; 
        Phi = inv(Phi1)*Phi2;   % Phi1\Phi2

        err1(cnt,it) = norm(Phi1-Phi1_bar,2);
        err2(cnt,it) = norm(Phi2-Phi2_bar,2);
        err3(cnt,it) = norm(Phi-Phi_bar,2);
    end

end

fig2 = [mean(err1,1);mean(err2,1);mean(err3,1);std(err1,0,1);std(err2,0,1);std(err3,0,1)];

figure()
plot(coeff,log(mean(err1,1)))
hold on;
plot(coeff,log(mean(err2,1)))
hold off;
legend('Phi1','Phi2')

figure()
plot(coeff,std(err1,0,1))
hold on;
plot(coeff,std(err2,0,1))
hold off;
legend('Phi1','Phi2')