
clear ;
clc;
coeff =1:5:40;
N = 40 ;
T = 80 ;
% theta_true = [0,5*2*pi/N];
theta_true = [0,pi/3];
P = [1 0.4; 0.4 1];

SNR = 2;

% 实例化所有对象
ArrayObject = [];
for ii = 1:length(coeff)
    ArrayObject = [ArrayObject ArraySignalModel(N*coeff(ii),T*coeff(ii),theta_true,P,SNR)];
end

nbLoop = 100;
% 跟Loop 有关的变量
ESPRITDoA_Nb = zeros(nbLoop,2);
ESPRITMSE_Nb =  zeros(nbLoop,1);
ESPRITEiValue_Nb= zeros(nbLoop,2);
GESPRITDoA_Nb = zeros(nbLoop,2);
GESPRITMSE_Nb =  zeros(nbLoop,1);
GESPRITEiValue_Nb = zeros(nbLoop,2);

% 跟自变量SNRList有关的变量
ESPRITMSE =  zeros(1,length(coeff));
ESPRITVar =  zeros(1,length(coeff));
ESPRITBias =  zeros(1,length(coeff));
GESPRITMSE = zeros(1,length(coeff));
GESPRITVar =  zeros(1,length(coeff));
GESPRITBias =  zeros(1,length(coeff));
CRB_Res         = zeros(1,length(coeff));
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [ESPRITDoA_Nb(Loop_i,:),ESPRITMSE_Nb(Loop_i,1),ESPRITEiValue_Nb(Loop_i,:)]  = ObjectNow.GetESPRIT();
        [GESPRITDoA_Nb(Loop_i,:),GESPRITMSE_Nb(Loop_i,1),GESPRITEiValue_Nb(Loop_i,:)]  = ObjectNow.GetGESPRIT('Theory');

    end
    [ESPRITMSE(1,object_i),ESPRITVar(1,object_i),ESPRITBias(1,object_i)] = ObjectNow.GetStatNum(ESPRITDoA_Nb,ESPRITMSE_Nb);
    [GESPRITMSE(1,object_i),GESPRITVar(1,object_i),GESPRITBias(1,object_i)] = ObjectNow.GetStatNum(GESPRITDoA_Nb,GESPRITMSE_Nb);
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end
% theta_true
% ESPRITDoA
% ObjectNow = ArrayObject(1);
%% 实验部分
% 实验DOA

figure;
hold on ;
subplot(1,2,1)
hold on ;
% xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(40*coeff,log10(ESPRITMSE),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(40*coeff,log10(ESPRITVar),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(40*coeff,log10(ESPRITBias),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
legend('MSE','Var','Bias')
title('ESPRIT')
xlabel('N')

subplot(1,2,2)
hold on ;
% xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(40*coeff,log10(GESPRITMSE),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(40*coeff,log10(GESPRITVar),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(40*coeff,log10(GESPRITBias),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
legend('MSE','Var','Bias')
title('GESPRIT')
xlabel('N')
