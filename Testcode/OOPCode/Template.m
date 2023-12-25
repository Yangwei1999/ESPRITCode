
clear ;
clc;
coeff =2;
N = 40 * coeff;
T = 80 * coeff;
theta_true = [0,5*2*pi/N];
P = [1 0.4; 0.4 1];

SNRList = 2:15;
ScanArea = [-pi/2 pi/2];
ScanPrec = 4000;

% 变量
VariableList = SNRList;

% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,VariableList(ii))];
end

nbLoop = 20;
% 跟Loop 有关的变量
ESPRITDoA_Nb = zeros(nbLoop,2);
ESPRITMSE_Nb =  zeros(nbLoop,1);
ESPRITEiValue_Nb= zeros(nbLoop,2);

GESPRITDoA_Nb = zeros(nbLoop,2);
GESPRITMSE_Nb =  zeros(nbLoop,1);
GESPRITEiValue_Nb = zeros(nbLoop,2);

MUSICDoA_Nb = zeros(nbLoop,2);
MUSICMSE_Nb =  zeros(nbLoop,1);
MUSICEiValue_Nb= zeros(nbLoop,2);

GMUSICDoA_Nb = zeros(nbLoop,2);
GMUSICMSE_Nb =  zeros(nbLoop,1);
GMUSICEiValue_Nb= zeros(nbLoop,2);
% 跟自变量SNRList有关的变量
ESPRITMSE =  zeros(1,length(VariableList)); 
ESPRITVar =  zeros(1,length(VariableList));
ESPRITBias =  zeros(1,length(VariableList));

GESPRITMSE = zeros(1,length(VariableList));
GESPRITVar =  zeros(1,length(VariableList));
GESPRITBias =  zeros(1,length(VariableList));

MUSICMSE =  zeros(1,length(VariableList)); 
MUSICVar =  zeros(1,length(VariableList));
MUSICBias =  zeros(1,length(VariableList));

GMUSICMSE =  zeros(1,length(VariableList)); 
GMUSICVar =  zeros(1,length(VariableList));
GMUSICBias =  zeros(1,length(VariableList));

CRB_Res         = zeros(1,length(VariableList));
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [ESPRITDoA_Nb(Loop_i,:),ESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetESPRIT();
        [GESPRITDoA_Nb(Loop_i,:),GESPRITMSE_Nb(Loop_i,1),~]  = ObjectNow.GetGESPRIT('Empirical-1');
        [MUSICDoA_Nb(Loop_i,:),MUSICMSE_Nb(Loop_i,1),~]  = ObjectNow.GetMusic(ScanArea,ScanPrec);
        [GMUSICDoA_Nb(Loop_i,:),GMUSICMSE_Nb(Loop_i,1),~]  = ObjectNow.GetGMusic(ScanArea,ScanPrec);
        % [~,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-2');
    end
    [ESPRITMSE(1,object_i),ESPRITVar(1,object_i),ESPRITBias(1,object_i)] = ObjectNow.GetStatNum(ESPRITDoA_Nb,ESPRITMSE_Nb);
    [GESPRITMSE(1,object_i),GESPRITVar(1,object_i),GESPRITBias(1,object_i)] = ObjectNow.GetStatNum(GESPRITDoA_Nb,GESPRITMSE_Nb);
    [MUSICMSE(1,object_i),MUSICVar(1,object_i),MUSICBias(1,object_i)] = ObjectNow.GetStatNum(MUSICDoA_Nb,MUSICMSE_Nb);
    [GMUSICMSE(1,object_i),GMUSICVar(1,object_i),GMUSICBias(1,object_i)] = ObjectNow.GetStatNum(GMUSICDoA_Nb,GMUSICMSE_Nb);
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end


figure;
hold on ;
xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(VariableList,log10(ESPRITMSE),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList,log10(GESPRITMSE),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList,log10(MUSICMSE),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(VariableList,log10(GMUSICMSE),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(VariableList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
legend('threshold','ESPRIT','DESPRIT','MUSIC','GMUSIC','CRB');
% legend('threshold','ESPRIT','GESPRIT','CRB')


if(1)
    figure;
    subplot(1,2,1)
    hold on ;
    plot(VariableList,log10(ESPRITMSE),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
    plot(VariableList,log10(ESPRITVar),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
    plot(VariableList,log10(ESPRITBias),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
    legend('MSE','Var','Bias')
    title('ESPRIT')
    xlabel('N')
%     
    subplot(1,2,2)
    hold on ;
    plot(VariableList,log10(GESPRITMSE),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
    plot(VariableList,log10(GESPRITVar),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
    plot(VariableList,log10(GESPRITBias),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
    legend('MSE','Var','Bias')
    title('GESPRIT')
    xlabel('N')
end

