clear; clc

%% 探究二阶矩阵 P非对角 时候得          特征值

coeff =20;
nb_Loop =1;
P = [2,0.4;0.4,1];
SNR =5;
sigma2 = 10.^(-SNR/10);

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
    Phi =  inv(U_S' * J1' * J1 * U_S) * (U_S' * J1' * J2 * U_S);
%     Phi1 = U_S' * J1' * J1 * U_S;
    Phi2 = U_S' * J1' * J2* U_S;
    [U,eigs_Phi2] = eig(Phi2,'vector');
    [eigs_Phi2, index] = sort(eigs_Phi2,'descend');
    eigs_Phi2
    eigs_Phi2_abs = abs(eigs_Phi2)
    eigs_Phi2_phase = angle(eigs_Phi2)

%     alpha1 = (g1+g2)/2 * n/N*(exp(1i*theta_true(1))+exp(1i*theta_true(2)));
%     alpha1 = n/N*(g1*exp(1i*theta_true(1)) + g2*exp(1i*theta_true(2)))
%     alpha2 = g1*g2*(n/N)*exp(1i*theta_true(1)) * (n/N) * exp(1i*theta_true(2))
% 
%     sqrtDelta  = sqrt(alpha1.^2 /4 - alpha2)
%     lameda_1 =  alpha1/2 + sqrtDelta
%     lameda_2  =  alpha1/2 - sqrtDelta

    [U_APA,eigs_APA] = eig((A*sqrtm(P))*(A*sqrtm(P))','vector');
    [eigs_APA, index] = sort(eigs_APA,'descend');
    U_APA = U_APA(:, index);
    eigs_APA(1:k);

    [U_P,eigs_P] = eig(P,'vector');
    [eigs_P, index] = sort(eigs_P,'descend');
    B = U_P(:, index);
    eigs_P = eigs_P(1:k)

    g1 = (1- c * (eigs_P(1)/sigma2)^(-2))/(1 + c * (eigs_P(1)/sigma2)^(-1));
    g2 = (1- c * (eigs_P(2)/sigma2)^(-2))/(1 + c * (eigs_P(2)/sigma2)^(-1));


    u1 = U_APA(:,1);
    u2 = U_APA(:,2);

    u1_hat = U_S(:,1);
    u2_hat = U_S(:,2);


    u1_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u1 * g1
    1/2 * (n/N *exp(1i*theta_true(1)) +n/N *exp(1i*theta_true(2))) * g1

    u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat
    u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2

    u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2

    u1_hat' * J1' *J2 * u1_hat * u2_hat' * J1' *J2 * u2_hat - u1_hat' * J1' *J2 * u2_hat * u2_hat' * J1' *J2 * u1_hat
    u1' * J1' *J2 * u1 * u2' * J1' *J2 * u2 *g1 * g2 -u1' * J1' *J2 * u2 * u2' * J1' *J2 * u1 *g1 * g2
    g1*g2 *det(U_APA(:,1:2)'*J1'*J2*U_APA(:,1:2))
    g1*g2 *det(A'*J1'*J2*A)
    g1*g2 * (n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2))

    Alpha1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))+(g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(i*theta_true(2));
    Alpha2 = g1*g2*(n/N)^2 * exp(1i*theta_true(1)) * exp(1i*theta_true(2));

    Delta = Alpha1^2 - 4 * Alpha2;

    lameda1 = (Alpha1 + sqrt(Delta))/2
    lameda2 = (Alpha1 - sqrt(Delta))/2

    (g1*B(1,1)^2 + g2*B(2,1)^2)
    (g1*B(1,2)^2 + g2*B(2,2)^2)
    lameda1Abs = abs(lameda1)
    lameda2Abs = abs(lameda2)
    lameda1Angle = angle(lameda1)
    lameda2Angle = angle(lameda2)

    theta_true

    Guesslameda1 = (g1*B(1,1)^2 + g2*B(2,1)^2)*n/N*exp(1i*theta_true(1))
    Guesslameda2 =  (g1*B(1,2)^2+g2*B(2,2)^2)*n/N*exp(1i*theta_true(2))
    abs(Guesslameda1)
    abs(Guesslameda2)

    g1*g2
    (g1*B(1,1)^2 + g2*B(2,1)^2)*(g1*B(1,2)^2+g2*B(2,2)^2)
    g1^2 * B(1,1)^2*B(1,2)^2 + g1*g2*B(2,1)^2*B(1,2)^2 +g1*g2*B(1,1)^2*B(2,2)^2 + g2^2*B(1,2)^2*B(1,1)^2
    g1*g2+B(1,1)^2 * B(1,2)^2 *(g1 - g2)^2
   
    g1*g2
    B(1,1)^2 * B(1,2)^2 *(g1 - g2)^2
    figure;
    hold on ;

    b = [real(sqrt(Delta)) imag(sqrt(Delta))].';

    quiver(0,0,real(Alpha1)/2,imag(Alpha1)/2,0,'LineWidth',1,'Color','r','DisplayName','Alpha1');
    quiver(0,0,real(sqrt(Delta))/2,imag(sqrt(Delta))/2,0,'LineWidth',1,'Color','r','DisplayName','Sqrt(Delta)');

    quiver(0,0,real(Alpha1+sqrt(Delta))/2,imag(Alpha1+sqrt(Delta))/2,0,'LineWidth',1,'Color','b','LineStyle','-'...
        ,'DisplayName','lameda_1_true');
    quiver(0,0,real(Alpha1-sqrt(Delta))/2,imag(Alpha1-sqrt(Delta))/2,0,'LineWidth',1,'Color','b','LineStyle','-'...
        ,'DisplayName','lameda_2_true');


    quiver(0,0,real(Guesslameda1),imag(Guesslameda1),0,'LineWidth',1,'Color','g','LineStyle','--',...
        'DisplayName','lameda_1_guess');
    quiver(0,0,real(Guesslameda2),imag(Guesslameda2),0,'LineWidth',1,'Color','g','LineStyle','--',...
         'DisplayName','lameda_2_guess');

    axis equal
    legend();











%     [~,EigenValues] = eig(Phi2,'vector');

   


end

% delta = alpha1.^2 /4 - alpha2
% abs(delta)
% angle(delta)
% c = sqrt(delta)
% abs_C = abs(c)
% angle_c = angle(c)

% figure;
% hold on;
% % xline(0,'LineWidth',1,'LineStyle','--');  % condition
% % plot(coeffList*40,log10(ESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(coeffList*40,log10(ESPRIT_Bias_E),'LineStyle','-','Color','#0072BD','Marker','s','LineWidth',1.5);
% plot(coeffList*40,log10(ESPRIT_Var_E),'LineStyle','-','Color',	'#0072BD','Marker','o','LineWidth',1.5)
% 
% % plot(coeffList*40,log10(MUSIC_MSE_E),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
% plot(coeffList*40,log10(GMUSIC_Bias_E),'LineStyle','-','Color','#77AC30','Marker','s','LineWidth',1.5);
% plot(coeffList*40,log10(GMUSIC_Var_E),'LineStyle','-','Color',	'#77AC30','Marker','o','LineWidth',1.5)
% 
% % plot(coeffList*40,log10(GMUSIC_MSE_E),'LineStyle','-','Color','#4DBEEE','Marker','x','LineWidth',1.5)
% % plot(coeffList*40,log10(GMUSIC_Bias_E),'LineStyle','-','Color','#4DBEEE','Marker','o','LineWidth',1.5);
% % plot(coeffList*40,log10(GMUSIC_Var_E),'LineStyle','-','Color',	'#4DBEEE','Marker','s','LineWidth',1.5)
% 
% legend('ESPRIT Bias','ESPRIT Var',...
%    'GMUSIC Bias','GMUSIC Var')
% % axis([0 900  -12 2])
% xlabel('N')
% ylabel('log10')
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

