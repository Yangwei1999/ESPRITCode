

c = 0.5;

l = (sqrt(c):0.1:10);

x = l./1;

y1 = (1-c*x.^(-2))./(1+c*x.^(-1));

y2 = (c*(x.^2+2.*x+c))./(x.^2+c.*x).^2;


figure;


hold on;

plot(l,y1,'LineWidth',1);

plot(l,y2,'LineWidth',1);

xlabel('l')
xline(sqrt(c),'LineStyle','--')
yline(1,'LineStyle','--')
legend('y1','y2','sqrt(c)','1')