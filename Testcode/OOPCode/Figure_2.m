
clear ;
clc;
%% 初始化参数
coeff =10;
N = 40 * coeff;
T = 100 * coeff;
theta_true = [0,5*2*pi/N];
% theta_true = [0,pi/4];
k = length(theta_true);
P = [1 0.4; 0.4 1];

SNRList = 0:2:20;
ScanArea = [-pi/2 pi/2];
ScanPrec = 8000;

% 变量
VariableList = SNRList;
VariableLabel = 'SNR';

% legend 
ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};

%% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,VariableList(ii))];
end

nbLoop = 20;
% 跟Loop 有关的变量  
ReceivedNum1 = 4;
DoA_Nb = zeros(ReceivedNum1,nbLoop,k);
MSE_Nb =  zeros(ReceivedNum1,nbLoop);
EiValue_Nb= zeros(ReceivedNum1,nbLoop,k);

% 跟自变量VariableList有关的变量 
ReceivedNum2 = ReceivedNum1;
MSE_VList = zeros(ReceivedNum2,length(VariableList));
Var_VList = zeros(ReceivedNum2,length(VariableList));
Bias_VList = zeros(ReceivedNum2,length(VariableList));
CRB_Res = zeros(1,length(VariableList));

%%算法部分
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [DoA_Nb(1,Loop_i,:),MSE_Nb(1,Loop_i),EiValue_Nb(1,Loop_i,:)]  = ObjectNow.GetESPRIT();                
        [DoA_Nb(2,Loop_i,:),MSE_Nb(2,Loop_i),EiValue_Nb(2,Loop_i,:)]  = ObjectNow.GetGESPRIT('Empirical-2');   
%         [DoA_Nb(3,Loop_i,:),MSE_Nb(3,Loop_i)]  = ObjectNow.GetMusic(ScanArea,ScanPrec);
        if(ReceivedNum1==4)
            [DoA_Nb(4,Loop_i,:),MSE_Nb(4,Loop_i)]  = ObjectNow.GetGMusic(ScanArea,ScanPrec); 
            [DoA_Nb(3,Loop_i,:),MSE_Nb(3,Loop_i),...
             DoA_Nb(4,Loop_i,:),MSE_Nb(4,Loop_i)] = ObjectNow.GetMusicType(ScanArea,ScanPrec,'Empirical-2');
        end
    end

    for kk = 1:ReceivedNum2
        [MSE_VList(kk,object_i),Var_VList(kk,object_i),Bias_VList(kk,object_i)] = ObjectNow.GetStatNum(squeeze(DoA_Nb(kk,:,:)),MSE_Nb(kk,:));
    end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end


%% 画图部分
figure;
hold on ;
xline(ObjectNow.SepCondition,'LineWidth',1,'LineStyle','--' ,'Label', ...
    num2str(round(ObjectNow.SepCondition,-1)) ,'LabelVerticalAlignment','bottom',...
    'LabelOrientation','horizontal');  % condition
plot(VariableList,log10(MSE_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList,log10(MSE_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
if(ReceivedNum1==4)
    plot(VariableList,log10(MSE_VList(3,:)),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
    plot(VariableList,log10(MSE_VList(4,:)),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
    
    ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};
else
    ShowLegend ={'Threshold','ESPRIT','GESPRIT','CRB'};
end
plot(VariableList,log10(CRB_Res),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
legend(ShowLegend())


% figure;
% hold on;
% plot(VariableList,log10(MSE_VList(3,:)),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(VariableList,log10(MSE_VList(4,:)),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
% 
% figure;
% hold on;
% plot(VariableList,log10(MSE_VList(5,:)),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
% plot(VariableList,log10(MSE_VList(6,:)),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
