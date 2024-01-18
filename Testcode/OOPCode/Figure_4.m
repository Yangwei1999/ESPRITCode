% 我想知道不同的子阵列对MSE的影响  把BIAS  也画出来
clear ;
clc;
coeff =2;
N = 40 *coeff;
T = 80 *coeff;
theta_true = [0,5*2*pi/N];

% theta_true = [5*2*pi/N]
k = length(theta_true);
% theta_true = [0,pi/4];
% P = [1 0; 0 1];
P = eye(k);

ArrayObject = [];

SNRList = 0:14;
% ScanArea = [-pi/4 pi/4];
% ScanPrec = 5000;
for ii = 1:length(SNRList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,SNRList(ii))];
end
nbLoop =100;
TestNum = 7;
ESPRIT_Doa_Res = zeros(nbLoop,length(SNRList),TestNum,k);
ESPRIT_MSE_Res = zeros(nbLoop,length(SNRList),TestNum);
ESPRIT_Angle_Res = zeros(nbLoop,length(SNRList),TestNum,k);


CRB_Res         = zeros(1,length(SNRList));

for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        for Test_i  = 1: TestNum
                [ESPRIT_Doa_Res(Loop_i,object_i,Test_i,:),ESPRIT_MSE_Res(Loop_i,object_i,Test_i),~,ESPRIT_Angle_Res(Loop_i,object_i,Test_i,:)]  = ...
                    ObjectNow.GetESPRITsub(1,ObjectNow.N-Test_i,Test_i);
        end
    end

%     for Test_i  = 1: TestNum
%            [ESPRITMSE(Test_i,object_i),ESPRITVar(Test_i,object_i),ESPRITBias(Test_i,object_i)] = ObjectNow.GetStatNum(ESPRIT_Doa_Res(:,:,Test_i,:),ESPRIT_MSE_Res(:,:,Test_i));
%     end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end
label1 ={};
for ii = 1:TestNum
    label1{ii} = ['Delta ', num2str(ii) ];
end

% 求 \Delta theta 的方差
% Var1 = squeeze(var(ESPRIT_Angle_Res,0,1));
% Mean1  = mean(Var1,3);
% 
% para1 = 1./(1:TestNum);
% para1 = repmat(para1,[length(SNRList),1]);
% 
% Mean2 = Mean1 .* para1;


% figure ;
% hold on;
% plot(SNRList,log10(Mean1.'))
% legend(label1)
% title('Var(Delta*theta)')

% figure ;
% hold on;
% plot(SNRList,log10(Mean2.'))
% legend(label1)


% ESPRIT_MSE_E = squeeze(mean(ESPRIT_MSE_Res,1));

if(1)
    ESPRIT_MSE_E = squeeze(mean(ESPRIT_MSE_Res,1));

    figure;
    hold on ;
    plot(SNRList,log10(ESPRIT_MSE_E.'),'LineWidth',0.5)
%     legend('1','2','3', ...
%         '4','5','6','7','8','9','10')
    legend(label1)
    title('MSE(theta)')
% legend('2')
end


