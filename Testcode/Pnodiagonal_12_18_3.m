clear; clc

%% 探究重根情况下的元素极限和特征值极限

coeff =10;
nb_Loop =100;
P = [1,0;0,1];
% P = [1,0.4;0.4,1];
SNR =5;
sigma2 = 10.^(-SNR/10);
N = 40* coeff ;
T  = 80* coeff;
c = N/T;
theta_true = [0,2*pi/N * 0.4];
k = length(theta_true);
clear i
a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
A = [];
for tmp_index=1:length(theta_true)
    A = [A a(theta_true(tmp_index))];
end

EigenValue_Result = complex(zeros(k,nb_Loop),zeros(k,nb_Loop));
EigenValue_Result2 = complex(zeros(k,nb_Loop),zeros(k,nb_Loop));

% 求解g1 g2
[U_APA,eigs_APA] = eig((A*sqrtm(P))*(A*sqrtm(P))','vector');
[eigs_APA, index] = sort(eigs_APA,'descend');
% U_APA = U_APA(:, index);
eigs_APA = eigs_APA(1:k);
g1 = (1- c * (eigs_APA(1)/sigma2)^(-2))/(1 + c * (eigs_APA(1)/sigma2)^(-1));
g2 = (1- c * (eigs_APA(2)/sigma2)^(-2))/(1 + c * (eigs_APA(2)/sigma2)^(-1));

Phi1_Res = zeros(nb_Loop,k,k);
Phi2_Res = zeros(nb_Loop,k,k);

Det_Res = zeros(nb_Loop,1);
Tr_Res = zeros(nb_Loop,1);
Det_Res_R = zeros(nb_Loop,1);
Tr_Res_R = zeros(nb_Loop,1);
Test_1 = zeros(nb_Loop,1);
Test_2 = zeros(nb_Loop,1);

for jj = 1 : nb_Loop
    jj
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
    Phi1 = (U_S' * J1' * J1 * U_S) ;
    Phi1_Res(jj,:,:) = Phi1;
    Phi2 = (U_S' * J1' * J2 * U_S);
    Phi2_Res(jj,:,:) = Phi2;

    Phi =  inv(Phi1)*Phi2;                  % 传统ESPRIT方法
    [U,eigs_Phi] = eig(Phi,'vector');
    [~, index] = sort(angle(eigs_Phi),'ascend');
    eigs_Phi = eigs_Phi(index);
    EigenValue_Result(:,jj) = eigs_Phi;


    Det_Res(jj) = det(Phi);
    Tr_Res(jj)  = trace(Phi);

    Det_Res_R(jj) = det(Phi)/(g1*g2);
    Tr_Res_R(jj)  = Phi2(1,1)/g1 + Phi2(2,2)/g2;
% 
    u1 = U_APA(:,1);
    u2 = U_APA(:,2);
    
    pinv(A) * U_APA(:,1:2);
    u1_hat = U_S(:,1);
    u2_hat = U_S(:,2);

    Test_1(jj) = u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat;
    Test_2(jj) = u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat;


end


Phi1_Res_E = squeeze(mean(Phi1_Res,1)) 
Phi2_Res_E = squeeze(mean(Phi2_Res,1))


trace(Phi2_Res_E)
det(Phi2_Res_E)

Tr_Res_E   = mean(Tr_Res)
Det_Res_E  = mean(Det_Res)
n/N*(exp(1i*theta_true(1)) + exp(1i*theta_true(2)))
n/N*exp(1i*theta_true(1)) *n/N*exp(1i*theta_true(2))

Det_Res_R_E = mean(Det_Res_R)
Tr_Res_R_E  = mean(Tr_Res_R)

Test_1_E    =mean(Test_1)
Test_2_E    =mean(Test_2)

%     Test_1(jj) = u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat;
    u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2
%     
%     Test_2(jj) = u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2
% %     
% angle(eig(Phi2_Res_E))
% ！！特征值的极限的期望和期望的特征值 不等----
%     期望的行列式 和行列式的期望     不等
% trace(Phi2_Res_E)
% det(Phi2_Res_E)

% n/N*(exp(1i*theta_true(1)) + exp(1i*theta_true(2)))
% n/N*exp(1i*theta_true(1)) *n/N*exp(1i*theta_true(2))
% 
% n/N*(g1*exp(1i*theta_true(1))+g2*exp(1i*theta_true(2)))
% g1*n/N*exp(1i*theta_true(1)) *g2*n/N*exp(1i*theta_true(2))
% % 
% 
% % q = squeeze(mean(Phi2_Res,1))
% % 
% % angle(eig(Phi2_Res_E))
% % 
EigenValue_Result_E = mean(EigenValue_Result,2)
angle((EigenValue_Result_E))
theta_true
% % 
figure;
hold on ;
% 
quiver(0,0,real(g1*(n/N)*exp(1i*theta_true(1))),imag(g1*(n/N)*exp(1i*theta_true(1))),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
quiver(0,0,real(g1*(n/N)*exp(1i*theta_true(2))),imag(g1*(n/N)*exp(1i*theta_true(2))),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');

quiver(0,0,real(EigenValue_Result_E(1)),imag(EigenValue_Result_E(1)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
quiver(0,0,real(EigenValue_Result_E(2)),imag(EigenValue_Result_E(2)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
% 
% % 
axis equal;

legend('true_1','true_2','emperical_1','emperical_2')

