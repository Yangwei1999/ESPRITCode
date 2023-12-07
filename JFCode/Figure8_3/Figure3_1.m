clear; clc
coeffList = 1:5:50;
SNRList = 5;
nb_Loop =20;
MUSIC_Time = zeros(nb_Loop,1);
GMUSIC_Time = zeros(nb_Loop,1);
ESPRIT_Time = zeros(nb_Loop,1);
DESPRIT_Time = zeros(nb_Loop,1);
MUSIC_Time_coeff = zeros(1,length(coeffList));
GMUSIC_Time_coeff = zeros(1,length(coeffList));
ESPRIT_Time_coeff = zeros(1,length(coeffList));
DESPRIT_Time_coeff = zeros(1,length(coeffList));
P = [2,0;0,1];
for coeff_i = 1 : length(coeffList)
    coeff = coeffList(coeff_i);
    N = 40 * coeff;
    T  = 80* coeff;
    c = N/T;
    theta_true = [0,5*2*pi/N]; 
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
    sigma2 = 0.1;
    for jj = 1 : nb_Loop
        S = sqrt(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
        Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
        X = A*S + Z;
        SCM = X*(X')/T;
        [U,eigs_SCM] = eig(SCM,'vector');
        [eigs_SCM, index] = sort(eigs_SCM,'descend');
        U = U(:, index);
        U_S = U(:,1:k);
        tic;
        MUSIC_Theta = GetMusic(U_S);
        MUSIC_Time(jj,1) = toc;
        tic;
        GMUSIC_Theta = GetGMusic(U_S,eigs_SCM,c);
        GMUSIC_Time(jj,1) = toc;
        tic;
        ESPRIT_Theta = GetESPRITE(U_S);
        ESPRIT_Time(jj,1) =toc;
        tic;
        DESPRIT_Theta = GetDESPRITE(U_S);
        DESPRIT_Time(jj,1) =toc;
    end
    MUSIC_Time_coeff(1,coeff_i) =   mean(MUSIC_Time,1);
    GMUSIC_Time_coeff(1,coeff_i) =  mean(GMUSIC_Time,1);
    ESPRIT_Time_coeff(1,coeff_i) =  mean(ESPRIT_Time,1);
    DESPRIT_Time_coeff(1,coeff_i) =  mean(DESPRIT_Time,1);

end
% 
figure;
hold on;
plot(40*coeffList,(MUSIC_Time_coeff),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(40*coeffList,(GMUSIC_Time_coeff),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(40*coeffList,(ESPRIT_Time_coeff),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(40*coeffList,(DESPRIT_Time_coeff),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
legend('MUSIC','GMUSIC','ESPRIT','DESPRIT')
xlabel('N')
ylabel('Time(s)')
title('$N/T = c $' , 'Interpreter','latex')

