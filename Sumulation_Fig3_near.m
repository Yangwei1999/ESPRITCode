close all; clear; clc

N = 40;
T  = 80;
c = N/T;

theta_true = [0,0.4*2*pi/N]; % near angle
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
SNRList = 16:34;
nb_Loop = 400;
CRB_Theory = zeros(1,length(SNRList));
MUSIC_Error = zeros(nb_Loop,length(SNRList));
GMUSIC_Error = zeros(nb_Loop,length(SNRList));
ESPRIT_Error = zeros(nb_Loop,length(SNRList));
for SNR_i = 1 : length(SNRList)
    sigma2 = 10.^(-SNRList(SNR_i)/10);
    for jj = 1 : nb_Loop
        S = sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
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
        %% error
        MUSIC_Error(jj,SNR_i) = (MUSIC_Theta(1) - theta_true(1)).^2;
        GMUSIC_Error(jj,SNR_i) = (GMUSIC_Theta(1) - theta_true(1)).^2;
        ESPRIT_Error(jj,SNR_i) =(ESPRIT_Theta(1) - theta_true(1)).^2;
        ESPRITML_Error(jj,SNR_i) = (ESPRIT_Theta_ML(1) - theta_true(1)).^2;
    end
    CRB = sigma2 / (2*T) *inv(real(D'*(eye(N)-A*inv(A'*A)*A')*D) .*eye(k));
    CRB_Theory(1,SNR_i) = CRB(1,1);
end
figure;
MUSIC_Error_E = mean(MUSIC_Error,1);
GMUSIC_Error_E = mean(GMUSIC_Error,1);
ESPRIT_Error_E = mean(ESPRIT_Error,1);
hold on;
plot(SNRList,log10(MUSIC_Error_E))
Mask1 = SNRList >=0;   % sepeartion condition 
plot(SNRList(Mask1),log10(GMUSIC_Error_E(Mask1)));
plot(SNRList,log10(ESPRIT_Error_E))
plot(SNRList,log10(CRB_Theory))
legend('MUSIC','GMUSIC','ESPRITE','CRB');








