clear; clc
% 探究Bias 和Var
coeffList = 1:5 :25;
nb_Loop = 100;
theta_true = [0,pi/3]; 
k = length(theta_true);
P = [1,0.4;0.4,1];
SNR = 2;
sigma2 = 10.^(-SNR/10);

ESPRIT_MSE = zeros(nb_Loop,length(coeffList));
ESPRIT_Bias = zeros(nb_Loop,length(coeffList));
ESPRIT_Var = zeros(nb_Loop,length(coeffList));


GESPRIT_MSE = zeros(nb_Loop,length(coeffList));
GESPRIT_Bias = zeros(nb_Loop,length(coeffList));
GESPRIT_Var = zeros(nb_Loop,length(coeffList));


% 求解g1 g2
[U_P,eigs_P] = eig(P,'vector');
[eigs_P, index] = sort(eigs_P,'descend');
B = U_P(:, index);
eigs_P = eigs_P(1:k);

for coeff_i = 1: length(coeffList)
    disp(coeff_i)
    coeff = coeffList(coeff_i);
    N = 40 * coeff;
    T  = 80* coeff;
    c = N/T;
    k = length(theta_true);
    clear i
    a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
    diffa = @(theta) exp(1i*theta*(0:N-1)') *1i .*(0:N-1).' /sqrt(N);
    A = [];
    diffA = [];
    for tmp_index=1:length(theta_true)
        A = [A a(theta_true(tmp_index))];
        diffA = [diffA diffa(theta_true(tmp_index))];
    end
    D = diffA;
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
        Phi1 = (U_S' * J1' * J1 * U_S) ;
        Phi2 = (U_S' * J1' * J2 * U_S) ;
        Phi =  inv(Phi1)*Phi2;       

        % 传统ESPRIT方法
        [~,eigs_Phi] = eig(Phi,'vector');
        [~, index] = sort(angle(eigs_Phi),'ascend');
        eigs_Phi = eigs_Phi(index);
%         EigenValue_Result(:,jj) = eigs_Phi;
        ESPRIT_DoA = angle(eigs_Phi);
        ESPRIT_DoA = sort(ESPRIT_DoA).';
%         ESPRIT_MSE(jj,SNR_i) = norm(ESPRIT_DoA - theta_true).^2 / k;

        % 修正参数 （ESPRIT-Modified）
        g1 = (1- c * (eigs_P(1)/sigma2)^(-2))/(1 + c * (eigs_P(1)/sigma2)^(-1));
        g2 = (1- c * (eigs_P(2)/sigma2)^(-2))/(1 + c * (eigs_P(2)/sigma2)^(-1));
        Pra1 = [1/g1 1/sqrt(g1*g2);1/sqrt(g1*g2) 1/g2 ];
        Phi2_test = Phi2 .* Pra1; 
        Phi =  inv(Phi1)*Phi2_test;
        [U,eigs_Phi] = eig(Phi,'vector');
        [~, index] = sort(angle(eigs_Phi),'ascend');
        eigs_Phi = eigs_Phi(index);
%         EigenValue_Result2(:,jj) = eigs_Phi;
        GESPRIT_DoA = angle(eigs_Phi);
        GESPRIT_DoA = sort(GESPRIT_DoA).';
%         GESPRIT_MSE(jj,SNR_i) = norm(GESPRIT_DoA - theta_true).^2 / k;


        ESPRIT_MSE(jj,coeff_i) = sum((ESPRIT_DoA - theta_true).^2) / k;
        ESPRIT_Bias(jj,coeff_i)  = sum(abs(ESPRIT_DoA - theta_true)) / k;
        ESPRIT_Var(jj,coeff_i)   = ESPRIT_MSE(jj,coeff_i) - ESPRIT_Bias(jj,coeff_i).^2;

        GESPRIT_MSE(jj,coeff_i) = sum((GESPRIT_DoA - theta_true).^2) / k;
        GESPRIT_Bias(jj,coeff_i)  = sum(abs(GESPRIT_DoA - theta_true)) / k;
        GESPRIT_Var(jj,coeff_i)   = GESPRIT_MSE(jj,coeff_i) - GESPRIT_Bias(jj,coeff_i).^2;

%         MUSIC_MSE(jj,coeff_i) = sum((MUSIC_Theta - theta_true).^2) / k;
%         MUSIC_Bias(jj,coeff_i)  = sum(abs(MUSIC_Theta - theta_true)) / k;
%         MUSIC_Var(jj,coeff_i)   = MUSIC_MSE(jj,coeff_i) - MUSIC_Bias(jj,coeff_i).^2;
% 
%         GMUSIC_MSE(jj,coeff_i) = sum((GMUSIC_Theta - theta_true).^2) / k;
%         GMUSIC_Bias(jj,coeff_i)  = sum(abs(GMUSIC_Theta - theta_true)) / k;
%         GMUSIC_Var(jj,coeff_i)   = GMUSIC_MSE(jj,coeff_i) - GMUSIC_Bias(jj,coeff_i).^2;

    end
   
end

ESPRIT_MSE_E  = mean(ESPRIT_MSE,1);
ESPRIT_Bias_E = mean(ESPRIT_Bias,1);
ESPRIT_Var_E  = mean(ESPRIT_Var,1);

GESPRIT_MSE_E  = mean(GESPRIT_MSE,1);
GESPRIT_Bias_E = mean(GESPRIT_Bias,1);
GESPRIT_Var_E  = mean(GESPRIT_Var,1);



figure;
subplot(1,2,1)
hold on;
plot(coeffList*40,log10(ESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','X','LineWidth',1.5);
plot(coeffList*40,log10(ESPRIT_Bias_E),'LineStyle','-','Color','#0072BD','Marker','s','LineWidth',1.5);
plot(coeffList*40,log10(ESPRIT_Var_E),'LineStyle','-','Color',	'#0072BD','Marker','o','LineWidth',1.5)
legend('MSE','Bias','Var')
title('ESPRIT')

% figure;
subplot(1,2,2)
hold on;
plot(coeffList*40,log10(GESPRIT_MSE_E),'LineStyle','-','Color','#77AC30','Marker','X','LineWidth',1.5);
plot(coeffList*40,log10(GESPRIT_Bias_E),'LineStyle','-','Color','#77AC30','Marker','s','LineWidth',1.5);
plot(coeffList*40,log10(GESPRIT_Var_E),'LineStyle','-','Color',	'#77AC30','Marker','o','LineWidth',1.5)
legend('MSE','Bias','Var')
title('GESPRIT')
% legend('ESPRIT Bias','ESPRIT Var',...
%    'GESPRIT Bias','GESPRIT Var')
% % axis([0 900  -12 2])
% xlabel('N')
% ylabel('log10')
% title('$N/T = 0.5$','Interpreter','latex')