%% Figure 1
clear;
load('3.mat')
figure;
hold on ;
subplot(2,1,1)
hold on;
%thoery
xline(Angle_Lit_ESPRIT(1),'LineWidth',1,'Color','#D95319','LineStyle','--')
xline(Angle_Lit_ESPRIT(2),'LineWidth',1,'Color','#D95319','LineStyle','--')
%Emperical
xline(Angle_Emp_ESPRIT(1),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
xline(Angle_Emp_ESPRIT(2),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
%true
xline(Angle_Lit_GESPRIT(1),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
xline(Angle_Lit_GESPRIT(2),'LineWidth',2,'Color','#A2142F','LineStyle','-')
legend('theory-1','theory-2','Emp-1','Emp-2','True-1','True-2')
% annotation('doublearrow',Angle_Lit_GESPRIT,[0.5,0.5])
title('Tradition ESPRIT')
axis([-0.2 1.2 0 2])

subplot(2,1,2)
hold on;
%thoery
xline(Angle_Lit_GESPRIT(1),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
xline(Angle_Lit_GESPRIT(2),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
%Emperical
xline(Angle_Emp_GESPRIT(1),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
xline(Angle_Emp_GESPRIT(2),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
axis([-0.2 1.2 0 2])
legend('theory(True)-1','theory(True)-2','Emp-1','Emp-2')
title('Improved ESPRIT(GESPRIT)')


%%Figure 2
clear;
load("4.mat")

figure;
subplot(1,2,1)
hold on ;
plot(VariableList*N,log10(MSE_VList(1,:)),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList*N,log10(Var_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList*N,log10(Bias_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
legend('MSE','Var','Bias')
title('ESPRIT')
xlabel(VariableLabel)
axis([min(VariableList*N) max(VariableList*N) -7 0])
%     
subplot(1,2,2)

hold on ;
plot(VariableList*N,log10(MSE_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList*N,log10(Var_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList*N,log10(Bias_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
legend('MSE','Var','Bias')
title('GESPRIT')
xlabel(VariableLabel)
axis([min(VariableList*N) max(VariableList*N) -7 0])



