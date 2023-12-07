close all; clear; clc
% 画CRB随着N的变换(T 固定)
T  = 80;
N_List = 100:50:500;
SNRList = 5;
res = zeros(length(N_List),length(SNRList));
theta_true = 0; 
k = length(theta_true);
for NList_i = 1:length(N_List)
    N  = N_List(NList_i);
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
    res(NList_i,:) = CRB_Theory;
end

figure;
hold on;
plot(SNRList,log10(res(1,:)),'LineStyle','-','Color','#0072BD','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(3,:)),'LineStyle','-','Color','#D95319','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(5,:)),'LineStyle','-','Color','#EDB120','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(7,:)),'LineStyle','-','Color','#7E2F8E','LineWidth',1.5,'Marker','s')
plot(SNRList,log10(res(9,:)),'LineStyle','-','Color','#77AC30','LineWidth',1.5,'Marker','s')
legend('$N = 100$', ...
    '$N = 200$', ...
    '$N = 300$', ...
    '$N = 400$', ...
    '$N = 500$', ...
    'interpreter','latex');
xlabel('SNR');
ylabel('CRB')
title('$T = 80$','Interpreter','latex')