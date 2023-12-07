close all; clear; clc
coeff = 2;
N = 40 * coeff;
T  = 80* coeff;
c = N/T;

theta_true = [0,0.4*2*pi/N];   % closely spaced doa
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
SNRList = 5:30;
nb_Loop = 100;
CRB_Theory = zeros(1,length(SNRList));
MUSIC_Error = zeros(nb_Loop,length(SNRList));
GMUSIC_Error = zeros(nb_Loop,length(SNRList));
ESPRIT_Error = zeros(nb_Loop,length(SNRList));
DESPRIT_Error = zeros(nb_Loop,length(SNRList));
P = [2,0;0,1];

for SNR_i = 1 : length(SNRList)
    sigma2 = 10.^(-SNRList(SNR_i)/10);
    for jj = 1 : nb_Loop
        S = sqrtm(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
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
    CRB = sigma2 / (2*T) *inv(real(D'*(eye(N)-A*inv(A'*A)*A')*D) .*P);
    CRB_Theory(1,SNR_i) = trace(CRB)/k;
end
figure;
MUSIC_Error_E = mean(MUSIC_Error,1);
GMUSIC_Error_E = mean(GMUSIC_Error,1);
ESPRIT_Error_E = mean(ESPRIT_Error,1);
DESPRIT_Error_E = mean(DESPRIT_Error,1);
hold on;
xline(0,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(MUSIC_Error_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GMUSIC_Error_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5);
plot(SNRList,log10(ESPRIT_Error_E),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(DESPRIT_Error_E),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(SNRList,log10(CRB_Theory),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
legend('threshold','MUSIC','GMUSIC','ESPRIT','DESPRIT','CRB');
xlabel('SNR','Interpreter','latex')
ylabel('MSE(log10)','Interpreter','latex')

