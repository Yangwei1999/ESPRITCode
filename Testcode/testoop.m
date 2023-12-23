
clear ;
clc;
coeff =4;
N = 40 * coeff;
T = 80 * coeff;
theta_true = [0*2*pi/N,5*2*pi/N];
P = [1 0.4 ; 0.4 1];

ArrayObject = [];

SNRList = 0:14;

for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
end
nbLoop = 100;
% h = waitbar(0,'sad');
ESPRIT_MSE_Res = zeros(nbLoop,length(SNRList));
GESPRIT_MSE_Res = zeros(nbLoop,length(SNRList));
CRB_Res         = zeros(1,length(SNRList));
% GESPRIT2_MSE_Res = zeros(nbLoop,length(SNRList));
for object_i = 1:length(ArrayObject)
    object_i
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop       
        ObjectNow.GenerateGuass();
        [~,ESPRIT_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetESPRIT();
        [~,GESPRIT_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-1');
        % [~,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-2');
    end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end

ESPRIT_MSE_E = mean(ESPRIT_MSE_Res,1);
GESPRIT_MSE_E = mean(GESPRIT_MSE_Res,1);
CRB_Res_E     = CRB_Res;
% GESPRIT2_MSE_E = mean(GESPRIT2_MSE_Res,1);
figure;
hold on ;

xline(2,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(ESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GESPRIT_MSE_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(SNRList,log10(CRB_Res_E),'LineStyle','--','Color','#0072BD','Marker','none','LineWidth',1.5)
% xline(2)
legend('threshold','ESPRIT','GESPRIT','CRB')
