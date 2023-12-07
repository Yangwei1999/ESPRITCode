close all; clear; clc
% 探究N固定T增长时 CRB的特性
N = 40;
T_List = 100:400:2200;
SNRList = -4:10;
res = zeros(length(T_List),length(SNRList));
theta_true = 0; 
k = length(theta_true);
for TList_i = 1:length(T_List)
    T  = T_List(TList_i);
    c = N/T;
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
    nb_Loop = 20;
    CRB_Theory = zeros(1,length(SNRList));
    MUSIC_Error = zeros(nb_Loop,length(SNRList));
    GMUSIC_Error = zeros(nb_Loop,length(SNRList));
    ESPRIT_Error = zeros(nb_Loop,length(SNRList)); 
    CRB_NbLoop = zeros(1,length(nb_Loop));
    P = 1;
    for SNR_i = 1 : length(SNRList)
        sigma2 = 10.^(-SNRList(SNR_i)/10);
        for jj = 1 : nb_Loop
            S = sqrt(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
            Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
            X = A*S + Z;
            SCM = X*(X')/T;
            [U,eigs_SCM] = eig(SCM,'vector');
            [eigs_SCM, index] = sort(eigs_SCM,'descend');
            U = U(:, index);
            U_S = U(:,1:k);
            CRB_Sum = zeros(k,k);
            for T_i = 1: T
                Xt = diag(S(:,T_i));
                CRB_Sum = CRB_Sum + real(Xt' * D' *(eye(N) - A*inv(A'*A)*A')*D*Xt);
            end
            CRB_Sum = sigma2 *inv(CRB_Sum) /2;
            CRB_NbLoop(1,jj) = CRB_Sum;
        end
        CRB_NbLoop_mean = mean(CRB_NbLoop,2);
        CRB_Theory(1,SNR_i) = trace(CRB_NbLoop_mean)/k;
    end
    res(TList_i,:) = CRB_Theory;
end

figure;
hold on;
plot(SNRList,log10(res(1,:)),'LineStyle','-','Color','#0072BD','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(2,:)),'LineStyle','-','Color','#D95319','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(3,:)),'LineStyle','-','Color','#EDB120','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(4,:)),'LineStyle','-','Color','#7E2F8E','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(5,:)),'LineStyle','-','Color','#77AC30','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(6,:)),'LineStyle','-','Color','#4DBEEE','LineWidth',1.5,'Marker','s')
legend('$T = 100$', ...
    '$T = 500$', ...
    '$T = 900$', ...
    '$T = 1300$', ...
    '$T = 1700$', ...
    '$T = 2100$', ...
    'interpreter','latex');
xlabel('SNR');
ylabel('CRB(log10)')
title('$N = 40$','Interpreter','latex')

