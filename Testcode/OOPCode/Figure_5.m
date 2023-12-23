
clear ;
clc;
coeff =2;
N = 40 * coeff;
T = 80 * coeff;
theta_true = [0,0.4*2*pi/N]; % closely spaced doar
P = [1 0; 0 1];

ArrayObject = [];

SNRList = 15:34;
ScanArea = [-pi/4 pi/4];
ScanPrec = 5000;
for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
end
nbLoop = 100;
ESPRIT_MSE_Res = zeros(nbLoop,length(SNRList));
GESPRIT_MSE_Res = zeros(nbLoop,length(SNRList));
MUSIC_MSE_Res = zeros(nbLoop,length(SNRList));
GMUSIC_MSE_Res = zeros(nbLoop,length(SNRList));
CRB_Res         = zeros(1,length(SNRList));
% GESPRIT2_MSE_Res = zeros(nbLoop,length(SNRList));
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [~,ESPRIT_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetESPRIT();
        [~,GESPRIT_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-1');
        [~,MUSIC_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetMusic(ScanArea,ScanPrec);
        [~,GMUSIC_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGMusic(ScanArea,ScanPrec);
        % [~,GESPRIT2_MSE_Res(Loop_i,object_i),~]  = ObjectNow.GetGESPRIT('Empirical-2');
    end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end

ESPRIT_MSE_E = mean(ESPRIT_MSE_Res,1);
GESPRIT_MSE_E = mean(GESPRIT_MSE_Res,1);
MUSIC_MSE_E = mean(MUSIC_MSE_Res,1);
GMUSIC_MSE_E = mean(GMUSIC_MSE_Res,1);
CRB_Res_E     = CRB_Res;
% GESPRIT2_MSE_E = mean(GESPRIT2_MSE_Res,1);
figure;
hold on ;

% xline(9.2,'LineWidth',1,'LineStyle','--');  % condition
plot(SNRList,log10(ESPRIT_MSE_E),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GESPRIT_MSE_E),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(SNRList,log10(MUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','x','LineWidth',1.5)
plot(SNRList,log10(GMUSIC_MSE_E),'LineStyle','-','Color','#0072BD','Marker','o','LineWidth',1.5)
plot(SNRList,log10(CRB_Res_E),'LineStyle','--','Color','#4DBEEE','LineWidth',1.5)
legend('ESPRIT','DESPRIT','MUSIC','GMUSIC','CRB');
% legend('threshold','ESPRIT','GESPRIT','CRB')
