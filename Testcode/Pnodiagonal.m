clear; clc

%% 探究二阶矩阵 P非对角 时候得          特征值

coeff =40;
nb_Loop =1;
P = [1,0.4;0.4,1];
SNR =2;
sigma2 = 10.^(-SNR/10);

N = 40* coeff ;
T  = 80* coeff;
c = N/T;
theta_true = [0,pi/3];
k = length(theta_true);
clear i
a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
A = [];
for tmp_index=1:length(theta_true)
    A = [A a(theta_true(tmp_index))];
end

for jj = 1 : nb_Loop
    S = sqrtm(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
    Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
    X = A*S + Z;
    SCM = X*(X')/T;
    [U,eigs_SCM] = eig(SCM,'vector');
    [eigs_SCM, index] = sort(eigs_SCM,'descend');
    U = U(:, index);
    U_S = U(:,1:k);
    J_tmp = eye(N);
    n = N-1;
    J1 = J_tmp(1:n,:);
    J2 = J_tmp(2:n+1,:);
    Phi =  inv(U_S' * J1' * J1 * U_S) * (U_S' * J1' * J2 * U_S);
%     Phi1 = U_S' * J1' * J1 * U_S;
    Phi2 = U_S' * J1' * J2* U_S;
    [U,eigs_Phi2] = eig(Phi2,'vector');
    [eigs_Phi2, index] = sort(eigs_Phi2,'descend');
    eigs_Phi2
    eigs_Phi2_abs = abs(eigs_Phi2)
    eigs_Phi2_phase = angle(eigs_Phi2)

%     alpha1 = (g1+g2)/2 * n/N*(exp(1i*theta_true(1))+exp(1i*theta_true(2)));
%     alpha1 = n/N*(g1*exp(1i*theta_true(1)) + g2*exp(1i*theta_true(2)))
%     alpha2 = g1*g2*(n/N)*exp(1i*theta_true(1)) * (n/N) * exp(1i*theta_true(2))
% 
%     sqrtDelta  = sqrt(alpha1.^2 /4 - alpha2)
%     lameda_1 =  alpha1/2 + sqrtDelta
%     lameda_2  =  alpha1/2 - sqrtDelta

    [U_APA,eigs_APA] = eig((A*sqrtm(P))*(A*sqrtm(P))','vector');
    [eigs_APA, index] = sort(eigs_APA,'descend');
    U_APA = U_APA(:, index);
    eigs_APA(1:k);
    [U_P,eigs_P] = eig(P,'vector');
    [eigs_P, index] = sort(eigs_P,'descend');
    B = U_P(:, index);
    eigs_P = eigs_P(1:k)

    g1 = (1- c * (eigs_P(1)/sigma2)^(-2))/(1 + c * (eigs_P(1)/sigma2)^(-1));
    g2 = (1- c * (eigs_P(2)/sigma2)^(-2))/(1 + c * (eigs_P(2)/sigma2)^(-1));


    u1 = U_APA(:,1);
    u2 = U_APA(:,2);

    u1_hat = U_S(:,1);
    u2_hat = U_S(:,2);


    u1_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u1 * g1
    1/2 * (n/N *exp(1i*theta_true(1)) +n/N *exp(1i*theta_true(2))) * g1

    u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat
    u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2

    u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2

    u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat - u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2 -u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2
    g1*g2 *det(U_APA(:,1:2)'*J1'*J2*U_APA(:,1:2))
    g1*g2 *det(A'*J1'*J2*A)
    g1*g2 * (n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2))

    Alpha1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))+(g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(i*theta_true(2));
    Alpha2 = g1*g2*(n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2));

    Delta = Alpha1^2 - 4 * Alpha2;

    lameda1 = (Alpha1 + sqrt(Delta))/2
    lameda2 = (Alpha1 - sqrt(Delta))/2

    (g1*B(1,1)^2 + g2*B(2,1)^2)
    (g1*B(1,2)^2 + g2*B(2,2)^2)
    lameda1Abs = abs(lameda1)
    lameda2Abs = abs(lameda2)
    lameda1Angle = angle(lameda1)
    lameda2Angle = angle(lameda2)

    theta_true

    Guesslameda1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))
    Guesslameda2 =  (g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(1i*theta_true(2))
    abs(Guesslameda1)
    abs(Guesslameda2)

    g1*g2
    (g1*B(1,1)^2 + g2*B(2,1)^2)*(g1*B(1,2)^2+g2*B(2,2)^2)
    g1^2 * B(1,1)^2*B(1,2)^2 + g1*g2*B(2,1)^2*B(1,2)^2 +g1*g2*B(1,1)^2*B(2,2)^2 + g2^2*B(1,2)^2*B(1,1)^2
    g1*g2+B(1,1)^2 * B(1,2)^2 *(g1 - g2)^2
   
    g1*g2
    B(1,1)^2 * B(1,2)^2 *(g1 - g2)^2    
    figure;
    hold on ;
%     quiver(0,0,real(Alpha1)/2,imag(Alpha1)/2,0,'LineWidth',1,'Color','r','DisplayName','Alpha1');
%     quiver(0,0,real(sqrt(Delta))/2,imag(sqrt(Delta))/2,0,'LineWidth',1,'Color','r','DisplayName','Sqrt(Delta)');

    quiver(0,0,real(Alpha1+sqrt(Delta))/2,imag(Alpha1+sqrt(Delta))/2,0,'LineWidth',1,'Color','b','LineStyle','-');
    quiver(0,0,real(Alpha1-sqrt(Delta))/2,imag(Alpha1-sqrt(Delta))/2,0,'LineWidth',1,'Color','b','LineStyle','-');

    quiver(0,0,real(Guesslameda1),imag(Guesslameda1),0,'LineWidth',1,'Color','g','LineStyle','--');
    quiver(0,0,real(Guesslameda2),imag(Guesslameda2),0,'LineWidth',1,'Color','g','LineStyle','--');

    axis equal
    legend('$\overline{\lambda}_{1}$(true)','$\overline{\lambda}_{2}(true)$','$\overline{\lambda}_{1}$(Guess)','$\overline{\lambda}_{2}$(Guess)','interpreter','latex');

end


