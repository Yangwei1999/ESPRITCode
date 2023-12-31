
clear ;
clc;
coeff =2;
N = 40 * coeff;
T = 80 * coeff;
theta_true = [0,5*2*pi/N];
P = [1 0.4; 0.4 1];

SNRList = 2:15;
ScanArea = [-pi/4 pi/4];
ScanPrec = 5000;

% 实例化所有对象
ArrayObject = [];
for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
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
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [~,ESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetESPRIT();
        [~,GESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetGESPRIT('Empirical-1');
%         [~,MUSIC_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetMusic(ScanArea,ScanPrec);
%         [~,GMUSIC_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGMusic(ScanArea,ScanPrec);
        % [~,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-2');
    end
    [ESPRITMSE(1,object_i),ESPRITVar(1,object_i),ESPRITBias(1,object_i)] = ObjectNow.GetStatNum(ESPRITDoA_Nb,ESPRITMSE_Nb);
    [GESPRITMSE(1,object_i),GESPRITVar(1,object_i),GESPRITBias(1,object_i)] = ObjectNow.GetStatNum(GESPRITDoA_Nb,GESPRITMSE_Nb);
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end

% ESPRIT_MSE_E = mean(ESPRIT_MSE_Res,1);
% GESPRIT_MSE_E = mean(GESPRIT_MSE_Res,1);
% MUSIC_MSE_E = mean(MUSIC_MSE_Res,1);
% GMUSIC_MSE_E = mean(GMUSIC_MSE_Res,1);
% CRB_Res_E     = CRB_Res;
% GESPRIT2_MSE_E = mean(GESPRIT2_MSE_Res,1);
figure;
hold on ;

xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(ESPRITMSE),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GESPRITMSE),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(MUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(GMUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(SNRList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
% legend('threshold','ESPRIT','DESPRIT','MUSIC','GMUSIC','CRB');
legend('threshold','ESPRIT','GESPRIT','CRB')
