clear; clc
% coeff = 5;

coeffList = [1:5:30];

% coeffList = 20;
nb_Loop =40;
P = [1,0.4;0.4,1];
SNR =8;
sigma2 = 10.^(-SNR/10);

ESPRIT_MSE = zeros(nb_Loop,length(coeffList));
ESPRIT_Bias = zeros(nb_Loop,length(coeffList));
ESPRIT_Var = zeros(nb_Loop,length(coeffList));

MUSIC_MSE = zeros(nb_Loop,length(coeffList));
MUSIC_Bias = zeros(nb_Loop,length(coeffList)); 
MUSIC_Var = zeros(nb_Loop,length(coeffList));

GMUSIC_MSE = zeros(nb_Loop,length(coeffList));
GMUSIC_Bias = zeros(nb_Loop,length(coeffList));
GMUSIC_Var = zeros(nb_Loop,length(coeffList));

for  coeff_i = 1: length(coeffList)
    coeff = coeffList(coeff_i);
    N = 40* coeff ;
    T  = 80* coeff;
    c = N/T;
    theta_true = [0,pi/3];
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
        disp(['coeff_i : '  num2str(coeff_i) ',' 'nb Loop : '  num2str(jj)])
        S = sqrtm(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
        Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
        X = A*S + Z;
        SCM = X*(X')/T;
        [U,eigs_SCM] = eig(SCM,'vector');
        [eigs_SCM, index] = sort(eigs_SCM,'descend');
        U = U(:, index);
        U_S = U(:,1:k);
% 
        J_tmp = eye(N);
         n = N-1;
        J1 = J_tmp(1:n,:);
        J2 = J_tmp(2:n+1,:);
        Phi =  inv(U_S' * J1' * J1 * U_S) * (U_S' * J1' * J2 * U_S);
        Phi1 = U_S' * J1' * J1 * U_S;
        Phi2 = U_S' * J1' * J2* U_S;
%         abs(Phi2)
%         angle(Phi2)
    g1 = (1- c * (1.4)^2)/(1 + c * (1.4)^(-1));
    g2 = (1- c * (0.6)^2)/(1 + c * (0.6)^(-1));
        [~,EigenValues] = eig(Phi2,'vector');
        abs_egin = abs(EigenValues);
        angle_egin = angle(EigenValues);
%         Res = sort(angle(EigenValues)).';


%         MUSIC_Theta = GetMusic(U_S);
%         GMUSIC_Theta = GetGMusic(U_S,eigs_SCM,c);
        ESPRIT_Theta = GetESPRITE(U_S);

        ESPRIT_MSE(jj,coeff_i) = sum((ESPRIT_Theta - theta_true).^2) / k;
        ESPRIT_Bias(jj,coeff_i)  = sum(abs(ESPRIT_Theta - theta_true)) / k;
        ESPRIT_Var(jj,coeff_i)   = ESPRIT_MSE(jj,coeff_i) - ESPRIT_Bias(jj,coeff_i).^2;
% 
%         MUSIC_MSE(jj,coeff_i) = sum((MUSIC_Theta - theta_true).^2) / k;
%         MUSIC_Bias(jj,coeff_i)  = sum(abs(MUSIC_Theta - theta_true)) / k;
%         MUSIC_Var(jj,coeff_i)   = MUSIC_MSE(jj,coeff_i) - MUSIC_Bias(jj,coeff_i).^2;
% % 
%         GMUSIC_MSE(jj,coeff_i) = sum((GMUSIC_Theta - theta_true).^2) / k;
%         GMUSIC_Bias(jj,coeff_i)  = sum(abs(GMUSIC_Theta - theta_true)) / k;
%         GMUSIC_Var(jj,coeff_i)   = GMUSIC_MSE(jj,coeff_i) - GMUSIC_Bias(jj,coeff_i).^2;

    end
   
end

ESPRIT_MSE_E  = mean(ESPRIT_MSE,1);
ESPRIT_Bias_E = mean(ESPRIT_Bias,1);
ESPRIT_Var_E  = mean(ESPRIT_Var,1);

% MUSIC_MSE_E  = mean(MUSIC_MSE,1);
% MUSIC_Bias_E = mean(MUSIC_Bias,1);
% MUSIC_Var_E  = mean(MUSIC_Var,1);
% 
% GMUSIC_MSE_E  = mean(GMUSIC_MSE,1);
% GMUSIC_Bias_E = mean(GMUSIC_Bias,1);
% GMUSIC_Var_E  = mean(GMUSIC_Var,1);

figure;
hold on;
% xline(0,'LineWidth',1,'LineStyle','--');  % condition
% plot(coeffList*40,log10(ESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(coeffList*40,log10(ESPRIT_Bias_E),'LineStyle','-','Color','#0072BD','Marker','s','LineWidth',1.5);
plot(coeffList*40,log10(ESPRIT_Var_E),'LineStyle','-','Color',	'#0072BD','Marker','o','LineWidth',1.5)

% plot(coeffList*40,log10(MUSIC_MSE_E),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
% plot(coeffList*40,log10(GMUSIC_Bias_E),'LineStyle','-','Color','#77AC30','Marker','s','LineWidth',1.5);
% plot(coeffList*40,log10(GMUSIC_Var_E),'LineStyle','-','Color',	'#77AC30','Marker','o','LineWidth',1.5)

% plot(coeffList*40,log10(GMUSIC_MSE_E),'LineStyle','-','Color','#4DBEEE','Marker','x','LineWidth',1.5)
% plot(coeffList*40,log10(GMUSIC_Bias_E),'LineStyle','-','Color','#4DBEEE','Marker','o','LineWidth',1.5);
% plot(coeffList*40,log10(GMUSIC_Var_E),'LineStyle','-','Color',	'#4DBEEE','Marker','s','LineWidth',1.5)

legend('ESPRIT Bias','ESPRIT Var',...
   'GMUSIC Bias','GMUSIC Var')
% axis([0 900  -12 2])
xlabel('N')
ylabel('log10')
% title('$N/T = 0.5$','Interpreter','latex')
%     'GMUSICMSE','GMUSICBiar','GMUSICVar')
% MUSIC_Error_E = mean(MUSIC_Error,1);
% GMUSIC_Error_E = mean(GMUSIC_Error,1);
% ESPRIT_Error_E = mean(ESPRIT_Error,1);
% DESPRIT_Error_E = mean(DESPRIT_Error,1);
% hold on;
% xline(0,'LineWidth',1,'LineStyle','--');  % condition
% plot(SNRList,log10(MUSIC_Error_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(GMUSIC_Error_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5);
% plot(SNRList,log10(ESPRIT_Error_E),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(DESPRIT_Error_E),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(CRB_Theory),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
% legend('threshold','MUSIC','GMUSIC','ESPRIT','DESPRIT','CRB');
% 
% xlabel('SNR','Interpreter','latex')
% ylabel('MSE(log10)','Interpreter','latex')


% g1 = (1- c * (1.4)^2)/(1 + c * (1.4)^(-1));
% g2 = (1- c * (0.6)^2)/(1 + c * (0.6)^(-1));
% a = 1 ;
% b = (g1 + g2)/2 * n/N*(exp(1i*0)+exp(1i*pi/3));
% c = g1*g2 * (n/N)^2 * exp(1i*0)*exp(1i*pi/3);
% 
% 
% deleta = b^2 - 4*a*c;
% 
% ans1 = (-b+sqrt(deleta))/(2*a)
% ans2 = (-b-sqrt(deleta))/(2*a)
% 
% 
% abs(ans1)
% angle(ans1)
% abs(ans2)
% angle(ans2)
