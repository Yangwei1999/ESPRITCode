
clear ;
clc;
% 初始settings  三个角度的情况
coeff =1:5:30;
N = 40 * 1;
T = 80 * 1;
% theta_true = [0,5*2*pi/N];
theta_true = [0,8/pi,pi/3];
k = length(theta_true);
% P = [1 0.4; 0.4 1];
P = [1 0.4 0;0.4 1 0; 0  0 1]
SNR = 2;

ScanArea = [-pi/2 pi/2];
ScanPrec = 4000;

% 变量
VariableList = coeff;
VariableLabel = 'N';
espritdoa
% legend 
ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};

%% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N*VariableList(ii),T*VariableList(ii),theta_true,P,SNR)];
end

nbLoop = 100;
% 跟Loop 有关的变量  
ReceivedNum1 = 2;
DoA_Nb = zeros(ReceivedNum1,nbLoop,k);
MSE_Nb =  zeros(ReceivedNum1,nbLoop);
EiValue_Nb= zeros(ReceivedNum1,nbLoop,k);

% 跟自变量VariableList有关的变量 
ReceivedNum2 = ReceivedNum1;
MSE_VList = zeros(ReceivedNum2,length(VariableList));
Var_VList = zeros(ReceivedNum2,length(VariableList));
Bias_VList = zeros(ReceivedNum2,length(VariableList));
EiValue_VList = zeros(2,length(VariableList),k);  % Only ESPRIT-Type methods
EiValue_VListTheory = zeros(2,length(VariableList),k);  % Only ESPRIT-Type methods
CRB_Res = zeros(1,length(VariableList));

%%代码部分

for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [DoA_Nb(1,Loop_i,:),MSE_Nb(1,Loop_i),EiValue_Nb(1,Loop_i,:)]  = ObjectNow.GetESPRIT();                
        [DoA_Nb(2,Loop_i,:),MSE_Nb(2,Loop_i),EiValue_Nb(2,Loop_i,:)]  = ObjectNow.GetGESPRIT('theory');   
%         [DoA_Nb(3,Loop_i,:),MSE_Nb(3,Loop_i)]                         = ObjectNow.GetMusic(ScanArea,ScanPrec);
%         [DoA_Nb(4,Loop_i,:),MSE_Nb(4,Loop_i)]                         = ObjectNow.GetGMusic(ScanArea,ScanPrec); 
    end

    for kk = 1:ReceivedNum2
        [MSE_VList(kk,object_i),Var_VList(kk,object_i),Bias_VList(kk,object_i)] = ObjectNow.GetStatNum(squeeze(DoA_Nb(kk,:,:)),MSE_Nb(kk,:));
    end
    % Emperical EiegnValue
    for kk = 1:2
        EiValue_VList(kk,object_i,:) = mean(squeeze(EiValue_Nb(kk,:,:)),1);
    end
%     % Theory   EigenValue
%     % 获取对象信息
%     U_APA = ObjectNow.UsTrue;
%     g = (1- ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-2))./...
%         (1 + ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-1));
%     J_tmp = eye(ObjectNow.N);
%     n = ObjectNow.N-1;
%     J1 = J_tmp(1:n,:);
%     J2 = J_tmp(2:end,:);
%     
%     % ESPRIT算法理论特征值
%     u1 = U_APA(:,1);
%     u2 = U_APA(:,2);
%     Alpha1 = g(1)  *  u1'*J1'*J2*u1 + g(2) * u2'*J1'*J2*u2;
%     Alpha2 = g(1)  *  g(2) *(n/ObjectNow.N).^2 * exp(1i * theta_true(1)) * exp(1i * theta_true(2));
%     Delta = Alpha1^2 - 4 * Alpha2;
%     EiValue_VListTheory(1,object_i,:) = [(Alpha1 + sqrt(Delta))/2 *ObjectNow.N/n (Alpha1 - sqrt(Delta))/2 *ObjectNow.N/n];
%     % GESPRIT算法理论特征值
%     EiValue_VListTheory(2,object_i,:) = (exp(1i*ObjectNow.ThetaTrue));
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end



% 
% figure;
% for ii = 1 : 6
%     subplot(2,3,ii)
%     hold on ;
%     % 
%     index1 = ii;
%     quiver(0,0,real(EiValue_VListTheory(1,index1,1)),imag(EiValue_VListTheory(1,index1,1)),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
%     quiver(0,0,real(EiValue_VListTheory(1,index1,2)),imag(EiValue_VListTheory(1,index1,2)),0,'LineWidth',1,'Color','#0072BD','LineStyle','--');
%     % 
%     quiver(0,0,real(EiValue_VList(1,index1,1)),imag(EiValue_VList(1,index1,1))            ,0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
%     quiver(0,0,real(EiValue_VList(1,index1,2)),imag(EiValue_VList(1,index1,2))            ,0,'LineWidth',1,'Color',	'#D95319','LineStyle','-');
%     
%     quiver(0,0,real(EiValue_VListTheory(2,index1,1)),imag(EiValue_VListTheory(2,index1,1)),0,'LineWidth',1,'Color','#EDB120','LineStyle','--');
%     quiver(0,0,real(EiValue_VListTheory(2,index1,2)),imag(EiValue_VListTheory(2,index1,2)),0,'LineWidth',1,'Color','#EDB120','LineStyle','--');
%     
%     quiver(0,0,real(EiValue_VList(2,index1,1)),imag(EiValue_VList(2,index1,1))            ,0,'LineWidth',1,'Color',	'#77AC30','LineStyle','-');
%     quiver(0,0,real(EiValue_VList(2,index1,2)),imag(EiValue_VList(2,index1,2))            ,0,'LineWidth',1,'Color',	'#77AC30','LineStyle','-');
% 
%     legend('$\overline{\lambda}_{1}$','$\overline{\lambda}_{2}$',...
%     '$\overline{\lambda}_{1}$','$\overline{\lambda}_{2}$',...
%     '$\overline{\lambda}_{1}$','$\overline{\lambda}_{2}$',...
%     '$\overline{\lambda}_{1}$','$\overline{\lambda}_{2}$',...
%     'interpreter','latex');
%     title(['N = ' num2str(N*VariableList(ii))])
% end

if(1)
    figure;
    subplot(1,2,1)
    hold on ;
    plot(VariableList*40,log10(MSE_VList(1,:)),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
    plot(VariableList*40,log10(Var_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
    plot(VariableList*40,log10(Bias_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
    legend('MSE','Var','Bias')
    title('ESPRIT')
    xlabel(VariableLabel)
%     
    subplot(1,2,2)
    hold on ;
    plot(VariableList*40,log10(MSE_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
    plot(VariableList*40,log10(Var_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
    plot(VariableList*40,log10(Bias_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
    legend('MSE','Var','Bias')
    title('GESPRIT')
    xlabel(VariableLabel)
end

if(1)
    figure;
    hold on ;
    plot(VariableList*40,log10(Var_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
    plot(VariableList*40,log10(Var_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
    legend('Var\_ESPRIT','Var\_GESPRIT')
    title('ESPRIT')
    xlabel(VariableLabel)
end




