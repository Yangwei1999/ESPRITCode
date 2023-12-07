%% 不同扰动的版本
close all; clear; clc
coeff = 2;
N = 40*coeff;
T  = 80*coeff;
c = N/T;
theta_true = [0,5*2*pi/N]; 
k = length(theta_true);
clear i
rho_noise_List = [0 0.001 0.01 0.1]*2*pi/N;
SNRList = -4:10;
MUSIC_Error_E = zeros(length(rho_noise_List),length(SNRList));
GMUSIC_Error_E = zeros(length(rho_noise_List),length(SNRList));
ESPRIT_Error_E = zeros(length(rho_noise_List),length(SNRList));
DESPRIT_Error_E = zeros(length(rho_noise_List),length(SNRList));

for rho_i = 1 : length(rho_noise_List)

    rho_noise_2 = rho_noise_List(rho_i);
    a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
    a_true = @(theta) exp(1i*(theta+sqrt(rho_noise_2)*randn(1,1))*(0:N-1)')/sqrt(N);
    diffa = @(theta) exp(1i*theta*(0:N-1)') *1i .*(0:N-1).' /sqrt(N);
    A = [];
    A_thoery = [];
    diffA = [];
    for tmp_index=1:length(theta_true)
        A = [A a_true(theta_true(tmp_index))];
        A_thoery = [A_thoery a(theta_true(tmp_index))];
        diffA = [diffA diffa(theta_true(tmp_index))];
    end
    D = diffA;
    nb_Loop = 100;
    CRB_Theory = zeros(1,length(SNRList));
    MUSIC_Error = zeros(nb_Loop,length(SNRList));
    GMUSIC_Error = zeros(nb_Loop,length(SNRList));
    ESPRIT_Error = zeros(nb_Loop,length(SNRList));
    DESPRIT_Error = zeros(nb_Loop,length(SNRList));
    P = [2,0;0,1];
    for SNR_i = 1 : length(SNRList)
        sigma2 =10.^(-SNRList(SNR_i)/10);
        for jj = 1 : nb_Loop
            S = sqrt(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
            Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
            X = A*S + Z;
            SCM = X*(X')/T;
            [U,eigs_SCM] = eig(SCM,'vector');
            [eigs_SCM, index] = sort(eigs_SCM,'descend');
            U = U(:, index);
            U_S = U(:,1:k);
            MUSIC_Theta = GetMusic(U_S);
            GMUSIC_Theta = GetGMusic(U_S,eigs_SCM,c);
            ESPRIT_Theta = GetESPRITE(U_S);
            DESPRIT_Theta = GetDESPRITE(U_S);
           %% error
            MUSIC_Error(jj,SNR_i) = norm(MUSIC_Theta - theta_true).^2 / k;
            GMUSIC_Error(jj,SNR_i) = norm(GMUSIC_Theta - theta_true).^2 /k;
            ESPRIT_Error(jj,SNR_i) =norm(ESPRIT_Theta - theta_true).^2 / k;
            DESPRIT_Error(jj,SNR_i) =norm(DESPRIT_Theta - theta_true).^2 / k;
    
        end
        CRB = sigma2 / (2*T) *inv(real(D'*(eye(N)-A_thoery*inv(A_thoery'*A_thoery)*A_thoery')*D) .*P);
        CRB_Theory(1,SNR_i) = trace(CRB)/k;
    end
    MUSIC_Error_E(rho_i,:)   = mean(MUSIC_Error,1);
    GMUSIC_Error_E(rho_i,:)  = mean(GMUSIC_Error,1);
    ESPRIT_Error_E(rho_i,:)  = mean(ESPRIT_Error,1);
    DESPRIT_Error_E(rho_i,:) = mean(DESPRIT_Error,1);

end
figure;
hold on;    

plot(SNRList,log10(MUSIC_Error_E(1,:)),'LineStyle','-','Color','#0072BD','LineWidth',1.5,'Marker','o')
plot(SNRList,log10(GMUSIC_Error_E(1,:)),'LineStyle','-','Color','#0072BD','LineWidth',1.5,'Marker','x')
plot(SNRList,log10(ESPRIT_Error_E(1,:)),'LineStyle','-','Color','#0072BD','LineWidth',1.5,'Marker','s')

plot(SNRList,log10(MUSIC_Error_E(3,:)),'LineStyle','-','Color','#77AC30','LineWidth',1.5,'Marker','o')
plot(SNRList,log10(GMUSIC_Error_E(3,:)),'LineStyle','-','Color','#77AC30','LineWidth',1.5,'Marker','x')
plot(SNRList,log10(ESPRIT_Error_E(3,:)),'LineStyle','-','Color','#77AC30','LineWidth',1.5,'Marker','s')

plot(SNRList,log10(MUSIC_Error_E(4,:)),'LineStyle','-','Color','#4DBEEE','LineWidth',1.5,'Marker','o')
plot(SNRList,log10(GMUSIC_Error_E(4,:)),'LineStyle','-','Color','#4DBEEE','LineWidth',1.5,'Marker','x')
plot(SNRList,log10(ESPRIT_Error_E(4,:)),'LineStyle','-','Color','#4DBEEE','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(CRB_Theory),'LineStyle','--','Color','cyan','LineWidth',1.5)
% 
legend('$\sigma^2 = 0$(MUSIC)',...
    '$\sigma^2 = 0$(GMUSIC)',...
    '$\sigma^2 = 0$(ESPRIT)',...
    '$\sigma^2 = 0.01$(MUSIC)',...
    '$\sigma^2 = 0.01$(GMUSIC)',...
    '$\sigma^2 = 0.01$(ESPRIT)',...
    '$\sigma^2 = 0.1$(MUSIC)',...
    '$\sigma^2 = 0.1$(GMUSIC)',...
    '$\sigma^2 = 0.1$(ESPRIT)',...
    'CRB','interpreter','latex');

xlabel('SNR')
ylabel('MSE(log10)')




