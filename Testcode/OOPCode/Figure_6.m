% 当 n = round(N)/2 时 情况
clear ;
clc;
coeff =2;
N = 40 * coeff;
T = 80 * coeff;
theta_true = [5*2*pi/N];
k = length(theta_true);;

% P = [1 0; 0 1];
P = 1;
SNRList =0:14;
ScanArea = [-pi/4 pi/4];
ScanPrec = 2000;

% 实例化所有对象
ArrayObject = [];
for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
end

nbLoop = 100;
% 跟Loop 有关的变量
ESPRITDoA_Nb = zeros(nbLoop,k);
ESPRITMSE_Nb =  zeros(nbLoop,1);
ESPRITEiValue_Nb= zeros(nbLoop,k);


ESPRITDoAsub_Nb = zeros(nbLoop,k);
ESPRITMSEsub_Nb =  zeros(nbLoop,1);
ESPRITEiValuesub_Nb= zeros(nbLoop,k);

GESPRITDoA_Nb = zeros(nbLoop,k);
GESPRITMSE_Nb =  zeros(nbLoop,1);
GESPRITEiValue_Nb = zeros(nbLoop,k);

% MUSICDoA_Nb = zeros(nbLoop,2);
% MUSICMSE_Nb =  zeros(nbLoop,1);
% MUSICEiValue_Nb= zeros(nbLoop,2);
% 
% GMUSICDoA_Nb = zeros(nbLoop,2);
% GMUSICMSE_Nb =  zeros(nbLoop,1);
% GMUSICEiValue_Nb= zeros(nbLoop,2);
% 跟自变量SNRList有关的变量
ESPRITMSE =  zeros(1,length(SNRList)); 
ESPRITVar =  zeros(1,length(SNRList));
ESPRITBias =  zeros(1,length(SNRList));

ESPRITMSESub =  zeros(1,length(SNRList)); 
ESPRITVarSub =  zeros(1,length(SNRList));
ESPRITBiasSub =  zeros(1,length(SNRList));

GESPRITMSE = zeros(1,length(SNRList));
GESPRITVar =  zeros(1,length(SNRList));
GESPRITBias =  zeros(1,length(SNRList));

MUSICMSE =  zeros(1,length(SNRList)); 
MUSICVar =  zeros(1,length(SNRList));
MUSICBias =  zeros(1,length(SNRList));

GMUSICMSE =  zeros(1,length(SNRList)); 
GMUSICVar =  zeros(1,length(SNRList));
GMUSICBias =  zeros(1,length(SNRList));

CRB_Res         = zeros(1,length(SNRList));
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [ESPRITDoA_Nb(Loop_i,:),ESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetESPRIT();
%         [GESPRITDoA_Nb(Loop_i,:),GESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetGESPRIT('Empirical-1');
        [ESPRITDoAsub_Nb4(Loop_i,:),ESPRITMSEsub_Nb4(Loop_i,1),~]  = ObjectNow.GetESPRITsub(1,ObjectNow.N-4,4);

        [ESPRITDoAsub_Nb(Loop_i,:),ESPRITMSEsub_Nb(Loop_i,1),~,angle]  = ObjectNow.GetESPRITsub(1,ObjectNow.N/2,ObjectNow.N/2);
%          espritdoa()
        % 找到角度中不为0的那个点 Angle  = delta * theta
        a = espritdoa(ObjectNow.SCMHat,1);
        pi* sind(a);
        % 修正周期2pi 角度
%         [~,index1] = max(abs(angle));
        % 扫描2*pi 周期中与真实角度最近的点
        k = -20:20;
        kList = 2*pi*k;
        ScanAngle = (angle(1) + kList)/(ObjectNow.N/2);
        % 找到距离真实角度最近的点
        [~,index] = min(abs(ScanAngle - theta_true));
        ScanAngle(index);

        RepairAngle = [ScanAngle(index)];

%         theta_true
        ESPRITMSEsub_Nb(Loop_i,1) = (RepairAngle - theta_true).^2;
        ESPRITDoAsub_Nb(Loop_i,:) = RepairAngle;
%         [MUSICDoA_Nb(Loop_i,:),MUSICMSE_Nb(Loop_i,1),~]  = ObjectNow.GetMusic(ScanArea,ScanPrec);
%         [GMUSICDoA_Nb(Loop_i,:),GMUSICMSE_Nb(Loop_i,1),~]  = ObjectNow.GetGMusic(ScanArea,ScanPrec);
        % [~,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-2');
    end
    [ESPRITMSE(1,object_i),ESPRITVar(1,object_i),ESPRITBias(1,object_i)] = ObjectNow.GetStatNum(ESPRITDoA_Nb,ESPRITMSE_Nb);
    [ESPRITMSESub4(1,object_i),ESPRITVarSub4(1,object_i),ESPRITBiasSub4(1,object_i)] =...
    ObjectNow.GetStatNum(ESPRITDoAsub_Nb4,ESPRITMSEsub_Nb4);
    [ESPRITMSESub(1,object_i),ESPRITVarSub(1,object_i),ESPRITBiasSub(1,object_i)] =...
        ObjectNow.GetStatNum(ESPRITDoAsub_Nb,ESPRITMSEsub_Nb);
%     [GESPRITMSE(1,object_i),GESPRITVar(1,object_i),GESPRITBias(1,object_i)] = ObjectNow.GetStatNum(GESPRITDoA_Nb,GESPRITMSE_Nb);
%     [MUSICMSE(1,object_i),MUSICVar(1,object_i),MUSICBias(1,object_i)] = ObjectNow.GetStatNum(MUSICDoA_Nb,MUSICMSE_Nb);
%     [GMUSICMSE(1,object_i),GMUSICVar(1,object_i),GMUSICBias(1,object_i)] = ObjectNow.GetStatNum(GMUSICDoA_Nb,GMUSICMSE_Nb);
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end



figure;
hold on ;

xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(ESPRITMSE),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(ESPRITMSESub4),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
plot(SNRList,log10(ESPRITMSESub),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(MUSICMSE),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(GMUSICMSE),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(SNRList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
title('MSE')
% legend('threshold','ESPRIT','GESPRIT','ESPRIT-7','CRB')
legend('threshold','ESPRIT-1','ESPRIT-4','ESPRIT-N*0.5','CRB');
% legend('threshold','ESPRIT','GESPRIT','CRB')


figure

hold on;

plot(SNRList,log10(ESPRITBias),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(ESPRITBiasSub4),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
plot(SNRList,log10(ESPRITBiasSub),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
title('bias\^2')
legend('ESPRIT-1','ESPRIT-4','ESPRIT-N*0.5');
% plot(SNRList,log10(MUSICMSE),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(SNRList,log10(GMUSICMSE),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
% plot(SNRList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)


figure;

hold on;
plot(SNRList,log10(ESPRITVar),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(ESPRITVarSub4),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
plot(SNRList,log10(ESPRITVarSub),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
title('var')

legend('ESPRIT-1','ESPRIT-4','ESPRIT-N*0.5');



