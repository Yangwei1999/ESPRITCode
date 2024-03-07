
clear ;
clc;
%% 初始化参数
coeff =2;
N = 40 * coeff;
T = 80 * coeff;
% theta_true = [0,5*2*pi/N];
theta_true = [0,pi/4];
k = length(theta_true);
P = [2 0; 0 1];

SNRList = -4:10;
ScanArea = [-pi/2 pi/2];
ScanPrec = 4000;

% 变量
VariableList = SNRList;
VariableLabel = 'SNR';

% legend 
ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};
% ShowLegend ={'Threshold','ESPRIT','GESPRIT','CRB'};
%% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,VariableList(ii))];
end

nbLoop = 100;
% 跟Loop 有关的变量  
ReceivedNum1 = 4;
DoA_Nb = zeros(ReceivedNum1,nbLoop,k);
MSE_Nb =  zeros(ReceivedNum1,nbLoop);
EiValue_Nb= zeros(ReceivedNum1,nbLoop,k);

% 跟自变量SNRList有关的变量 
ReceivedNum2 = ReceivedNum1;
MSE_VList = zeros(ReceivedNum2,length(VariableList));
Var_VList = zeros(ReceivedNum2,length(VariableList));
Bias_VList = zeros(ReceivedNum2,length(VariableList));
CRB_Res = zeros(1,length(VariableList));

%%代码部分

for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [DoA_Nb(1,Loop_i,:),MSE_Nb(1,Loop_i),EiValue_Nb(1,Loop_i,:)]  = ObjectNow.GetESPRIT();             
        [DoA_Nb(2,Loop_i,:),MSE_Nb(2,Loop_i),EiValue_Nb(2,Loop_i,:)]  = ObjectNow.GetGESPRIT('Empirical-2');   
        [DoA_Nb(3,Loop_i,:),MSE_Nb(3,Loop_i),...
         DoA_Nb(4,Loop_i,:),MSE_Nb(4,Loop_i)] = ObjectNow.GetMusicType(ScanArea,ScanPrec,'Empirical-2');
    end

    for kk = 1:ReceivedNum2
        [MSE_VList(kk,object_i),Var_VList(kk,object_i),Bias_VList(kk,object_i)] = ObjectNow.GetStatNum(squeeze(DoA_Nb(kk,:,:)),MSE_Nb(kk,:));
    end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end


figure;
hold on ;
xline(ObjectNow.SepCondition,'LineWidth',1,'LineStyle','--' ,'Label', ...
    num2str(round(ObjectNow.SepCondition,-1)) ,'LabelVerticalAlignment','bottom',...
    'LabelOrientation','horizontal');  % condition
plot(VariableList,log10(MSE_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList,log10(MSE_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList,log10(MSE_VList(3,:)),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(VariableList,log10(MSE_VList(4,:)),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(VariableList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
legend(ShowLegend())


