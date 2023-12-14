clear; clc

%% 探究二阶矩阵 P非对角
% ESPRIT算法GESPRIT算法的MSE 比较 随着SNR

coeff =5;
nb_Loop =100;
P = [1,0.4;0.4,1];

% SNR =2;
SNRList = 0:1:14;
% SNRList = 2;
% sigma2 = 10.^(-SNR/10);
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

EigenValue_Result = complex(zeros(k,nb_Loop),zeros(k,nb_Loop));
EigenValue_Result2 = complex(zeros(k,nb_Loop),zeros(k,nb_Loop));

% 求解g1 g2
[U_P,eigs_P] = eig(P,'vector');
[eigs_P, index] = sort(eigs_P,'descend');
B = U_P(:, index);
eigs_P = eigs_P(1:k);


ESPRIT_MSE = zeros(nb_Loop,length(SNRList));
GESPRIT_MSE = zeros(nb_Loop,length(SNRList));

for SNR_i = 1: length(SNRList)
    sigma2 = 10.^(-SNRList(SNR_i)/10);
    for jj = 1 : nb_Loop
        disp(['Snr' num2str(SNR_i) 'jj = ' num2str(jj)]);

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
        Phi2 = (U_S' * J1' * J2 * U_S) ;
        Phi =  inv(Phi1)*Phi2;       
    
        % 传统ESPRIT方法
        [~,eigs_Phi] = eig(Phi,'vector');
        [~, index] = sort(angle(eigs_Phi),'ascend');
        eigs_Phi = eigs_Phi(index);
        EigenValue_Result(:,jj) = eigs_Phi;
        ESPRIT_DoA = angle(eigs_Phi);
        ESPRIT_DoA = sort(ESPRIT_DoA).';
        ESPRIT_MSE(jj,SNR_i) = norm(ESPRIT_DoA - theta_true).^2 / k;


        % 修正参数 （ESPRIT-Modified）
        g1 = (1- c * (eigs_P(1)/sigma2)^(-2))/(1 + c * (eigs_P(1)/sigma2)^(-1));
        g2 = (1- c * (eigs_P(2)/sigma2)^(-2))/(1 + c * (eigs_P(2)/sigma2)^(-1));
        Pra1 = [1/g1 1/sqrt(g1*g2);1/sqrt(g1*g2) 1/g2 ];
        Phi2_test = Phi2 .* Pra1; 
        Phi =  inv(Phi1)*Phi2_test;
        [U,eigs_Phi] = eig(Phi,'vector');
        [~, index] = sort(angle(eigs_Phi),'ascend');
        eigs_Phi = eigs_Phi(index);
        EigenValue_Result2(:,jj) = eigs_Phi;
        GESPRIT_DoA = angle(eigs_Phi);
        GESPRIT_DoA = sort(GESPRIT_DoA).';
        GESPRIT_MSE(jj,SNR_i) = norm(GESPRIT_DoA - theta_true).^2 / k;
    end
end
ESPRIT_MSE_E = mean(ESPRIT_MSE,1);
GESPRIT_MSE_E = mean(GESPRIT_MSE,1);


figure;
hold on ;

xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(ESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
% xline(2)
legend('threshold','ESPRIT','GESPRIT')
% 特征值的渐进值(传统ESPRIT)

% [U_APA,eigs_APA] = eig((A*sqrtm(P))*(A*sqrtm(P))','vector');
% [eigs_APA, index] = sort(eigs_APA,'descend');
% U_APA = U_APA(:, index);
% [U_P,eigs_P] = eig(P,'vector');
% [eigs_P, index] = sort(eigs_P,'descend');
% B = U_P(:, index);
% eigs_P = eigs_P(1:k)
% 
% u1 = U_APA(:,1);
% u2 = U_APA(:,2);
% Alpha1 = g1  *  u1'*J1'*J2*u1 + g2 * u2'*J1'*J2*u2;
% Alpha2 = g1  *  g2 *(n/N).^2 * exp(1i * theta_true(1)) * exp(1i * theta_true(2));
% Delta = Alpha1^2 - 4 * Alpha2;
% lameda1_limit = (Alpha1 + sqrt(Delta))/2;
% lameda2_limit = (Alpha1 - sqrt(Delta))/2;
% 
% 
% % 特征值的渐进值(ESPRIT-modified)
% lameda1_modify_limit = (n/N) * exp(1i * theta_true(1));
% lameda2_modify_limit = (n/N) * exp(1i * theta_true(2));
% 
% % ESPRIT 的emprical eigenvalue E
% EigenValue_E = mean(EigenValue_Result,2).';
% 
% % ESPRIT-modified 的emprical eigenvalue E
% EigenValue_E2 = mean(EigenValue_Result2,2).';
% 
% 
% 
% 
% 
% figure;
% hold on ;
% 
% quiver(0,0,real(lameda1_limit),imag(lameda1_limit),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
% quiver(0,0,real(lameda2_limit),imag(lameda2_limit),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
% 
% quiver(0,0,real(EigenValue_E(1)),imag(EigenValue_E(1)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
% quiver(0,0,real(EigenValue_E(2)),imag(EigenValue_E(2)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
% 
% quiver(0,0,real(EigenValue_E2(1)),imag(EigenValue_E2(1)),0,'LineWidth',1,'Color','#EDB120','LineStyle','-');
% quiver(0,0,real(EigenValue_E2(2)),imag(EigenValue_E2(2)),0,'LineWidth',1,'Color','#EDB120','LineStyle','-');
% 
% quiver(0,0,real(lameda1_modify_limit),imag(lameda1_modify_limit),0,'LineWidth',1,'Color',	'#77AC30','LineStyle','--');
% quiver(0,0,real(lameda2_modify_limit),imag(lameda2_modify_limit),0,'LineWidth',1,'Color',	'#77AC30','LineStyle','--');
% axis equal
% legend('$\overline{\lambda}_{1}$(theory-ESPRIT)','$\overline{\lambda}_{2}$(theory-ESPRIT)',...
%     '$\overline{\lambda}_{1}$(Empirical-ESPRIT)','$\overline{\lambda}_{2}$(Empirical-ESPRIT)',...
%     '$\overline{\lambda}_{1}$(theory-GESPRIT)','$\overline{\lambda}_{2}$(theory-GESPRIT)',...
%     '$\overline{\lambda}_{1}$(Empirical-GESPRIT)','$\overline{\lambda}_{2}$(Empirical-GESPRIT)',...
%     'interpreter','latex');
% 
% 
% 













% 
% [U_APA,eigs_APA] = eig((A*sqrtm(P))*(A*sqrtm(P))','vector');
% [eigs_APA, index] = sort(eigs_APA,'descend');
% U_APA = U_APA(:, index);
% eigs_APA(1:k);
% [U_P,eigs_P] = eig(P,'vector');
% [eigs_P, index] = sort(eigs_P,'descend');
% B = U_P(:, index);
% eigs_P = eigs_P(1:k)
% 
% g1 = (1- c * (eigs_P(1)/sigma2)^(-2))/(1 + c * (eigs_P(1)/sigma2)^(-1));
% g2 = (1- c * (eigs_P(2)/sigma2)^(-2))/(1 + c * (eigs_P(2)/sigma2)^(-1));
% 
% u1 = U_APA(:,1);
% u2 = U_APA(:,2);
% 
% u1_hat = U_S(:,1);
% u2_hat = U_S(:,2);
% 
% u1_hat' * J1' *J2 * u1_hat;
% u1' * J1' *J2 * u1 * g1;
% 1/2 * (n/N *exp(1i*theta_true(1)) +n/N *exp(1i*theta_true(2))) * g1;
% 
% u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat;
% u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2;
% 
% u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat;
% u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2;;
% 
% u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat - u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat;
% u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2 -u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2;
% g1*g2 *det(U_APA(:,1:2)'*J1'*J2*U_APA(:,1:2));
% g1*g2 *det(A'*J1'*J2*A);
% g1*g2 * (n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2));
% Alpha1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))+(g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(i*theta_true(2));
% 
% (g1*g2*u1'*J1'*J2*u1*u2'*J1'*J2*u2 - g1*g2* u1'*J1'*J2*u2*u2'*J1'*J2*u1);;
% Alpha2 = g1*g2*(n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2));
% 
% Delta = Alpha1^2 - 4 * Alpha2;
% 
% EigenValue_E = mean(EigenValue_Result,2)
% lameda1 = (Alpha1 + sqrt(Delta))/2
% lameda2 = (Alpha1 - sqrt(Delta))/2
% 
% 
% Guesslameda1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))
% Guesslameda2 =  (g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(1i*theta_true(2))
% 
%     figure;
%     hold on ;
% %     quiver(0,0,real(Alpha1)/2,imag(Alpha1)/2,0,'LineWidth',1,'Color','r','DisplayName','Alpha1');
% %     quiver(0,0,real(sqrt(Delta))/2,imag(sqrt(Delta))/2,0,'LineWidth',1,'Color','r','DisplayName','Sqrt(Delta)');
% 
%     quiver(0,0,real(lameda1),imag(lameda1),0,'LineWidth',1,'Color','b','LineStyle','-');
%     quiver(0,0,real(lameda2),imag(lameda2),0,'LineWidth',1,'Color','b','LineStyle','-');
% 
%     quiver(0,0,real(EigenValue_E(1)),imag(EigenValue_E(1)),0,'LineWidth',1,'Color','r','LineStyle','-');
%     quiver(0,0,real(EigenValue_E(2)),imag(EigenValue_E(2)),0,'LineWidth',1,'Color','r','LineStyle','-');
% 
%     quiver(0,0,real(Guesslameda1),imag(Guesslameda1),0,'LineWidth',1,'Color','g','LineStyle','--');
%     quiver(0,0,real(Guesslameda2),imag(Guesslameda2),0,'LineWidth',1,'Color','g','LineStyle','--');
% 
%     axis equal
%     legend('$\overline{\lambda}_{1}$(true)','$\overline{\lambda}_{2}(true)$',...
%     '$E[{\lambda}_{1}]$','$E[{\lambda}_{2}]$',...
%     '$\overline{\lambda}_{1}$(Guess)','$\overline{\lambda}_{2}$(Guess)',...
%     'interpreter','latex');
% 





