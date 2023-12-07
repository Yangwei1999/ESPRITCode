%% MUSIC、distributed MUSIC算法MSE随子阵元个数变化关系
clear; close all;
derad = pi/180;      %角度->弧度
M = 2;               % 信源数目
theta = [-10.2,35.5];  % 待估计角度
snr = -20;            % 信噪比
sigma2 = 1/db2pow(snr);
K = 512;             % 快拍数
 
dd = 0.5;            % 阵元间距 
p = 4;               %分布式节点个数
qq = 5:1:11;
iter = 200;
var_MU = zeros(length(qq),iter);
var_dMU = zeros(length(qq),iter);
mse_MU = zeros(1,length(qq));
mse_dMU = zeros(1,length(qq));
for it = 1:1:length(qq)
    q = qq(it);
    N = p*q;               % 阵元个数    
    d=0:dd:(N-1)*dd;
    ddd = 0:dd:(q-1)*dd;
    A=exp(-1i*2*pi*d.'*sin(theta*derad));  %方向矢量
    x = zeros(q,K,p);
    RXX = zeros(q,q,p);
    %%%%构建信号模型%%%%%
    P = eye(2,2);

    for cnt = 1:1:iter
        S=sqrtm(P)*randn(M,K);             %信源信号，入射信号
        W = complex(randn(N,K), randn(N,K));
        X = A*S + sqrt(sigma2/2)*W;
        for i =1:1:p
            x(:,:,i) = A(1+(i-1)*q:i*q,:)*S+ sqrt(sigma2/2)*complex(randn(q,K), randn(q,K));
        end
        X1=x;
        SCM = X*X'/K;
        [U,eigs_SCM] = eig(SCM,'vector');
        [eigs_SCM, index] = sort(eigs_SCM,'descend');
        U = U(:, index);
        U_S = U(:,1:3);

        for j = 1:1:p
                RXX(:,:,j) = squeeze(X1(:,:,j))*squeeze(X1(:,:,j))'/K;
        end
        e  = zeros(q,p);
        e2  = zeros(q,p);
        e3  = zeros(q,p);
        e4  = zeros(q,p);
        for i = 1:1:p
            e(:,i)  = randn(q,1);
            for k = 1:1:100
                e(:,i) = squeeze(RXX(:,:,i))*e(:,i);
                e(:,i) = e(:,i) ./ vecnorm(e(:,i));
            end
        end
        for i = 1:1:p
            e2(:,i)  = randn(q,1);
            for k = 1:1:100
                e2(:,i) = (diag(ones(1,q))-e(:,i)*e(:,i)')*squeeze(RXX(:,:,i))*e2(:,i);
                e2(:,i) = e2(:,i) ./ vecnorm(e2(:,i));
            end
        end
        for i = 1:1:p
            e3(:,i)  = randn(q,1);
            for k = 1:1:100
                e3(:,i) = (diag(ones(1,q))-e(:,i)*e(:,i)'-e2(:,i)*e2(:,i)')*squeeze(RXX(:,:,i))*e3(:,i);
                e3(:,i) = e3(:,i) ./ vecnorm(e3(:,i));
            end
        end
        Pmusic = zeros(361,2);
        for iang = 1:361
            angle(iang)=(iang-181)/2;
            phim=derad*angle(iang);
            a=exp(-1i*2*pi*ddd*sin(phim)).';
            aa = @(theta) exp(-1i*2*pi*d*sin(theta)).';
            aver = 0;
            for i =1:1:p
                Es = e(:,i);
                Es2 = e2(:,i); 
                Es3 = e3(:,i); 
                aver = a'*Es*Es'*a + a'*Es2*Es2'*a +aver; 
            end
            Pmusic(iang,1) = real(N - aver);
            Pmusic(iang,2) =  real(N-aa(phim)'*U_S*(U_S')*aa(phim));
            
            Pmusic1=abs(Pmusic(:,1));
            Pmmax1=max(Pmusic1);
            Pmusic1=10*log10(Pmusic1/Pmmax1);            % 归一化处理
            Pmusic2=abs(Pmusic(:,2));
            Pmmax2=max(Pmusic2);
            Pmusic2=10*log10(Pmusic2/Pmmax2);            % 归一化处理           
        end
        [pks,locs] = findpeaks(-Pmusic1);
        [~, index] = sort(pks,'descend');
        locs = locs(index);
%         if length(locs)==2
%             locs = locs(1:2);
%             locs = [locs;1];
%         else
            locs = locs(1:M);
%         end
        MUSIC_estim = sort(derad*angle(locs));
        [pks,locs] = findpeaks(-Pmusic2);
        [~, index] = sort(pks,'descend');
        locs = locs(index);
        locs = locs(1:M);
        dMUSIC_estim = sort(derad*angle(locs));
        var_MU(it,cnt) = sum((MUSIC_estim-deg2rad(theta)).^2)/2;
        var_dMU(it,cnt) = sum((dMUSIC_estim-deg2rad(theta)).^2)/2;
    end
    mse_MU(it)= sum(var_MU(it,:))/iter;
    mse_dMU(it)= sum(var_dMU(it,:))/iter;
end
    
figure()
plot(qq,log10(mse_MU),'+-','Color',[rand(1,3)],'linewidth',2);
hold on;
plot(qq,log10(mse_dMU),'+-','Color',[rand(1,3)],'linewidth',2);
hold off;
legend('MUSIC','Distributed MUSIC');
xlabel('单个传感器的子阵元数目');
ylabel('log(MSE)');

%% MUSIC、distributed MUSIC算法DOA估计结果
%%%%%%%%%%%二、多个目标仅取均值
clear; close all;
%%%%%%%% MUSIC for Uniform Linear Array%%%%%%%%
derad = pi/180;      %角度->弧度
p = 4;               %分布式节点个数
q = 10;               %每个节点SENSOR个数
N = p*q;               % 阵元个数        
M = 3;               % 信源数目
theta = [-10,35,37];  % 待估计角度
snr = 0;            % 信噪比
sigma2 = 1/db2pow(snr);
K = 512;             % 快拍数
 
dd = 0.5;            % 阵元间距 
d=0:dd:(N-1)*dd;
ddd = 0:dd:(q-1)*dd;
A=exp(-1i*2*pi*d.'*sin(theta*derad));  %方向矢量
x = zeros(q,K,p);
RXX = zeros(q,q,p);
%%%%构建信号模型%%%%%
P = eye(3,3);
S=sqrtm(P)*randn(M,K);             %信源信号，入射信号
W = complex(randn(N,K), randn(N,K));
X = A*S + sqrt(sigma2/2)*W;
for i =1:1:p
    x(:,:,i) = A(1+(i-1)*q:i*q,:)*S+ sqrt(sigma2/2)*complex(randn(q,K), randn(q,K));
end
% X1=awgn(x,snr,'measured'); %将白色高斯噪声添加到信号中
X1=x;
% 计算协方差矩阵
% Rxx=X1*X1'/K;
% % 特征值分解
% [EV,D]=eig(Rxx);                   %特征值分解
% EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
% [EVA,I]=sort(EVA);                 %将特征值排序 从小到大
% EV=fliplr(EV(:,I));                % 对应特征矢量排序
SCM = X*X'/K;
[U,eigs_SCM] = eig(SCM,'vector');
[eigs_SCM, index] = sort(eigs_SCM,'descend');
U = U(:, index);
U_S = U(:,1:3);
            
for j = 1:1:p
        RXX(:,:,j) = squeeze(X1(:,:,j))*squeeze(X1(:,:,j))'/K;
end
e  = zeros(q,p);
e2  = zeros(q,p);
e3  = zeros(q,p);
e4  = zeros(q,p);
for i = 1:1:p
    e(:,i)  = randn(q,1);
    for k = 1:1:100
        e(:,i) = squeeze(RXX(:,:,i))*e(:,i);
        e(:,i) = e(:,i) ./ vecnorm(e(:,i));
    end
end
for i = 1:1:p
    e2(:,i)  = randn(q,1);
    for k = 1:1:100
        e2(:,i) = (diag(ones(1,q))-e(:,i)*e(:,i)')*squeeze(RXX(:,:,i))*e2(:,i);
        e2(:,i) = e2(:,i) ./ vecnorm(e2(:,i));
    end
end
for i = 1:1:p
    e3(:,i)  = randn(q,1);
    for k = 1:1:100
        e3(:,i) = (diag(ones(1,q))-e(:,i)*e(:,i)'-e2(:,i)*e2(:,i)')*squeeze(RXX(:,:,i))*e3(:,i);
        e3(:,i) = e3(:,i) ./ vecnorm(e3(:,i));
    end
end
Pmusic = zeros(361,2);
% 遍历每个角度，计算空间谱
for iang = 1:361
    angle(iang)=(iang-181)/2;
    phim=derad*angle(iang);
    a=exp(-1i*2*pi*ddd*sin(phim)).';
    aa = @(theta) exp(-1i*2*pi*d*sin(theta)).';
    aver = 0;
    for i =1:1:p
        Es = e(:,i);
        Es2 = e2(:,i); 
        Es3 = e3(:,i); 
        aver = a'*Es*Es'*a + a'*Es2*Es2'*a + a'*Es3*Es3'*a +aver; 
    end
    Pmusic(iang,1) = real(N - aver);
    Pmusic(iang,2) =  real(N-aa(phim)'*U_S*(U_S')*aa(phim));
end


figure()
Pmusic1=abs(Pmusic(:,1));
Pmmax1=max(Pmusic1)
Pmusic1=10*log10(Pmusic1/Pmmax1);            % 归一化处理
Pmusic2=abs(Pmusic(:,2));
Pmmax2=max(Pmusic2)
Pmusic2=10*log10(Pmusic2/Pmmax2);            % 归一化处理
h=plot(angle,Pmusic1,'Linewidth',1.5);
hold on;
plot(angle,Pmusic2,'Linewidth',1.5);
hold off;
legend('Distributed MUSIC','MUSIC');
xlabel('入射角/(degree)');
ylabel('空间谱/(dB)');
set(gca, 'XTick',[-90:30:90]);
grid on;
