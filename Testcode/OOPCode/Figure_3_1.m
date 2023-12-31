
clear ;
clc;
% 初始settings
coeff =10;
N = 40 * coeff;
T = 80 * coeff;
% theta_true = [0,5*2*pi/N];
theta_true = [0,pi/3];
P = [1 0.4; 0.4 1];
SNRList = 2;
nbLoop = 100;

% 实例化所有对象
ArrayObject = [];
for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
end

% 跟Loop 有关的变量
ESPRITDoA_Nb = zeros(nbLoop,2);
ESPRITMSE_Nb =  zeros(nbLoop,1);
ESPRITEiValue_Nb= zeros(nbLoop,2);
GESPRITDoA_Nb = zeros(nbLoop,2);
GESPRITMSE_Nb =  zeros(nbLoop,1);
GESPRITEiValue_Nb = zeros(nbLoop,2);


% 跟自变量SNRList有关的变量
ESPRITMSE =  zeros(1,length(SNRList));
ESPRITVar =  zeros(1,length(SNRList));
ESPRITBias =  zeros(1,length(SNRList));
GESPRITMSE = zeros(1,length(SNRList));
GESPRITVar =  zeros(1,length(SNRList));
GESPRITBias =  zeros(1,length(SNRList));
CRB_Res         = zeros(1,length(SNRList));

for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        % 跟loop 相关的写在此处
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [ESPRITDoA_Nb(Loop_i,:),ESPRITMSE_Nb(Loop_i,1),ESPRITEiValue_Nb(Loop_i,:)]  = ObjectNow.GetESPRIT();
        [GESPRITDoA_Nb(Loop_i,:),GESPRITMSE_Nb(Loop_i,1),GESPRITEiValue_Nb(Loop_i,:)]  = ObjectNow.GetGESPRIT('Theory');
%         [GESPRIT2DoA,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-1');
    end
    % 跟obejct相关的写在此处
    [ESPRITMSE(1,object_i),ESPRITVar(1,object_i),ESPRITBias(1,object_i)] = ObjectNow.GetStatNum(ESPRITDoA_Nb,ESPRITMSE_Nb);
    [GESPRITMSE(1,object_i),GESPRITVar(1,object_i),GESPRITBias(1,object_i)] = ObjectNow.GetStatNum(GESPRITDoA_Nb,GESPRITMSE_Nb);
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end


ObjectNow = ArrayObject(1);

%% 实验部分
% 实验DOA
% MSE
[ESPRITMSE(1,1),ESPRITVar(1,1),ESPRITBias(1,1)]
[GESPRITMSE(1,1),GESPRITVar(1,1),GESPRITBias(1,1)]
% 实验特征值
ESPRITEigenValue_E = (mean(ESPRITEiValue_Nb,1))
GESPRITEigenValue_E= (mean(GESPRITEiValue_Nb,1))

%% 理论部分
% 获取对象信息
U_APA = ObjectNow.UsTrue;
g = (1- ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-2))./...
    (1 + ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-1));
J_tmp = eye(ObjectNow.N);
n = ObjectNow.N-1;
J1 = J_tmp(1:n,:);
J2 = J_tmp(2:end,:);

% ESPRIT算法理论特征值
u1 = U_APA(:,1);
u2 = U_APA(:,2);
Alpha1 = g(1)  *  u1'*J1'*J2*u1 + g(2) * u2'*J1'*J2*u2;
Alpha2 = g(1)  *  g(2) *(n/N).^2 * exp(1i * theta_true(1)) * exp(1i * theta_true(2));
Delta = Alpha1^2 - 4 * Alpha2;
Lambda_Lit_ESPRIT = [(Alpha1 + sqrt(Delta))/2 *N/n (Alpha1 - sqrt(Delta))/2 *N/n]
% GESPRIT算法理论特征值
Lambda_Lit_GESPRIT = (exp(1i*ObjectNow.ThetaTrue))

figure;
hold on ;
% 
quiver(0,0,real(Lambda_Lit_ESPRIT(1)),imag(Lambda_Lit_ESPRIT(1)),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
quiver(0,0,real(Lambda_Lit_ESPRIT(2)),imag(Lambda_Lit_ESPRIT(2)),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
% 
quiver(0,0,real(ESPRITEigenValue_E(1)),imag(ESPRITEigenValue_E(1)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
quiver(0,0,real(ESPRITEigenValue_E(2)),imag(ESPRITEigenValue_E(2)),0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');

quiver(0,0,real(Lambda_Lit_GESPRIT(1)),imag(Lambda_Lit_GESPRIT(1)),0,'LineWidth',1,'Color','#EDB120','LineStyle','--');
quiver(0,0,real(Lambda_Lit_GESPRIT(2)),imag(Lambda_Lit_GESPRIT(2)),0,'LineWidth',1,'Color','#EDB120','LineStyle','--');

quiver(0,0,real(GESPRITEigenValue_E(1)),imag(GESPRITEigenValue_E(1)),0,'LineWidth',1,'Color',	'#77AC30','LineStyle','-');
quiver(0,0,real(GESPRITEigenValue_E(2)),imag(GESPRITEigenValue_E(2)),0,'LineWidth',1,'Color',	'#77AC30','LineStyle','-');

% annotation('textarrow',[.3,.6],[.7,.4],'String','ABC');
axis equal
legend('$\overline{\lambda}_{1}$(Theory-ESPRIT)','$\overline{\lambda}_{2}$(Theory-ESPRIT)',...
    '$\overline{\lambda}_{1}$(Empirical-ESPRIT)','$\overline{\lambda}_{2}$(Empirical-ESPRIT)',...
    '$\overline{\lambda}_{1}$(Theory-GESPRIT)','$\overline{\lambda}_{2}$(Theory-GESPRIT)',...
    '$\overline{\lambda}_{1}$(Empirical-GESPRIT)','$\overline{\lambda}_{2}$(Empirical-GESPRIT)',...
    'interpreter','latex');



% Lambda2_Lit_GESPRIT = (Alpha1 - sqrt(Delta))/2
% GESPRIT2_MSE_E = mean(GESPRIT2_MSE_Res,1);
% 
% CRB_Res_E     = CRB_Res;
% GESPRIT2_MSE_E = mean(GESPRIT2_MSE_Res,1);
% figure;
% hold on ;
% % 
% % xline(2,'LineWidth',1,'LineStyle','--');  % condition
% plot(SNRList,log10(ESPRIT_MSE_E),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(GESPRIT_MSE_E),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(GESPRIT2_MSE_E),'LineStyle','-','Color','#77AC30','Marker','d','LineWidth',1.5)
% % plot(SNRList,log10(MUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% % plot(SNRList,log10(GMUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(CRB_Res_E),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
% % legend('threshold','ESPRIT','DESPRIT','MUSIC','GMUSIC','CRB');
% % % legend('threshold','ESPRIT','GESPRIT','CRB')
