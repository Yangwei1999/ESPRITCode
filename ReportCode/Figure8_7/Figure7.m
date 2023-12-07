clear; clc
% 探究Bias 和Var
coeffList = 1:5 :25;
nb_Loop = 20;
P = [1,0.4;0.4,1];
SNR = 3;
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

for coeff_i = 1: length(coeffList)
    coeff = coeffList(coeff_i);
    N = 40 * coeff;
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

        ESPRIT_MSE(jj,coeff_i) = sum((ESPRIT_Theta - theta_true).^2) / k;
        ESPRIT_Bias(jj,coeff_i)  = sum(abs(ESPRIT_Theta - theta_true)) / k;
        ESPRIT_Var(jj,coeff_i)   = ESPRIT_MSE(jj,coeff_i) - ESPRIT_Bias(jj,coeff_i).^2;

        MUSIC_MSE(jj,coeff_i) = sum((MUSIC_Theta - theta_true).^2) / k;
        MUSIC_Bias(jj,coeff_i)  = sum(abs(MUSIC_Theta - theta_true)) / k;
        MUSIC_Var(jj,coeff_i)   = MUSIC_MSE(jj,coeff_i) - MUSIC_Bias(jj,coeff_i).^2;

        GMUSIC_MSE(jj,coeff_i) = sum((GMUSIC_Theta - theta_true).^2) / k;
        GMUSIC_Bias(jj,coeff_i)  = sum(abs(GMUSIC_Theta - theta_true)) / k;
        GMUSIC_Var(jj,coeff_i)   = GMUSIC_MSE(jj,coeff_i) - GMUSIC_Bias(jj,coeff_i).^2;

    end
   
end

ESPRIT_MSE_E  = mean(ESPRIT_MSE,1);
ESPRIT_Bias_E = mean(ESPRIT_Bias,1);
ESPRIT_Var_E  = mean(ESPRIT_Var,1);

MUSIC_MSE_E  = mean(MUSIC_MSE,1);
MUSIC_Bias_E = mean(MUSIC_Bias,1);
MUSIC_Var_E  = mean(MUSIC_Var,1);

GMUSIC_MSE_E  = mean(GMUSIC_MSE,1);
GMUSIC_Bias_E = mean(GMUSIC_Bias,1);
GMUSIC_Var_E  = mean(GMUSIC_Var,1);

figure;
hold on;
plot(coeffList*40,log10(ESPRIT_Bias_E),'LineStyle','-','Color','#0072BD','Marker','s','LineWidth',1.5);
plot(coeffList*40,log10(ESPRIT_Var_E),'LineStyle','-','Color',	'#0072BD','Marker','o','LineWidth',1.5)

% plot(coeffList*40,log10(GMUSIC_Bias_E),'LineStyle','-','Color','#77AC30','Marker','s','LineWidth',1.5);
% plot(coeffList*40,log10(GMUSIC_Var_E),'LineStyle','-','Color',	'#77AC30','Marker','o','LineWidth',1.5)
legend('ESPRIT Bias','ESPRIT Var',...
   'GMUSIC Bias','GMUSIC Var')
% axis([0 900  -12 2])
xlabel('N')
ylabel('log10')
title('$N/T = 0.5$','Interpreter','latex')

