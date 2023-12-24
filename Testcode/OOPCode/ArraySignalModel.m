
classdef ArraySignalModel < handle
%%
    properties
        A    % 导向矢量矩阵
        D    % 导向矢量求导矩阵
        S    % 信号接收矩阵
        Z    % 噪声接收矩阵
        k    % 信号个数
        ThetaTrue    % 实际信号角度  频率域
        P    % 信号功率矩阵
        N    % Array Num
        T    % Sanpshot Num
        c
        sigma2 % -log10(SNR)
        UsHat   %  Signal Space
        EigsHat      %  Emperical eigenvalues
        UsTrue       %  Gorund Signal Space
        EigsTrue     %  Gorund Signal eigenvalus
%         GridGap    
%         Gridinterval
        
    end

    methods
        function obj = ArraySignalModel(N,T,Theta_true,P,SNR)

            a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
            diffa = @(theta) exp(1i*theta*(0:N-1)') *1i .*(0:N-1).' /sqrt(N);

            obj.k = length(Theta_true) ;
            obj.N = N;
            obj.T = T;
            obj.c = N/T;
            obj.ThetaTrue = sort(Theta_true);
            % 生成A
            A = [];
            diffA = [];
            for tmp_index=1:length(Theta_true)
                A = [A a(Theta_true(tmp_index))];
                diffA = [diffA diffa(Theta_true(tmp_index))];
            end
            obj.D = diffA;
            obj.A = A;
            obj.P = P;

            % 求解真实信号子空间  与信号不关  不用考虑相关 不相关  也不考虑 widely or closely
            [U_APA,eigs_APA] = eig((obj.A*sqrtm(obj.P))*(obj.A*sqrtm(obj.P))','vector');
            [eigs_APA, index] = sort(eigs_APA,'descend');
            U_APA = U_APA(:, index);
            obj.UsTrue = U_APA(:,1:obj.k);
            eigs_APA = eigs_APA(1:obj.k);
            obj.EigsTrue = eigs_APA;

            sigma2 = 10.^(-SNR/10);
            obj.sigma2 =sigma2;

%             obj.GridGap = GridGap;
% 
%             intervalabs = pi/4;
%             obj.Gridinterval = [-intervalabs,+intervalabs];

            % GenerateGuass(obj);

            % % 生成S
            % obj.S = sqrtm(P)*sqrt(1/2) *(randn(k,T) + 1i *randn(k,T));
            % % 生成Z noise 
            % sigma2 = 10.^(-SNR/10);
            % Z = sqrt(sigma2/2) * (randn(N,T) + 1i* randn(N,T));
            % obj.Z = Z;
            % obj.sigma2 =sigma2;


        end

        function [DoA,MSE,Bias,EigenValues] = GetESPRITsub(obj,index1,index2,Bias)
            %% Origin esprit methods
                J_tmp = eye(obj.N);
                n =  index2 - index1;
                J1 = J_tmp(index1:index2,:);
                J2 = J_tmp(index1+Bias:index2+Bias,:);
                Phi1 = obj.UsHat' * (J1' * J1) * obj.UsHat;
                Phi2 = obj.UsHat' * (J1' * J2) * obj.UsHat;
                Phi =  inv(Phi1)* Phi2;
                [~,EigenValues] = eig(Phi,'vector');
                [DoA,index] = sort(angle(EigenValues));
                DoA = DoA.'/Bias  ;
                EigenValues = EigenValues(index);
                MSE = sum((DoA - obj.ThetaTrue).^2) / obj.k;
                Bias =  sum(abs(DoA - obj.ThetaTrue))  / obj.k;
        end

        
        function [DoA,MSE,EigenValues] = GetESPRIT(obj)
            %% Origin esprit methods
                J_tmp = eye(obj.N);
                n =  obj.N-1;
                J1 = J_tmp(1:n,:);
                J2 = J_tmp(2:n+1,:);
                Phi1 = obj.UsHat' * (J1' * J1) * obj.UsHat;
                Phi2 = obj.UsHat' * (J1' * J2) * obj.UsHat;
                Phi =  inv(Phi1)* Phi2;
                [~,EigenValues] = eig(Phi,'vector');
                [DoA,index] = sort(angle(EigenValues));
                DoA = DoA.' ;
                EigenValues = EigenValues(index);
                MSE = sum((DoA - obj.ThetaTrue).^2) / obj.k;
%                 Bias =  sum(abs(DoA - obj.ThetaTrue))  / obj.k;
        end
        function [DoA,MSE,EigenValues] = GetGESPRIT(obj,type)
            %% Gereral esprit methods

                switch type
                    case 'Theory'
                        % case 1 : 采用理论修正构造G 
                        g = (1- obj.c .* (obj.EigsTrue./obj.sigma2).^(-2))./(1 + obj.c .* (obj.EigsTrue./obj.sigma2).^(-1));
            
                    case 'Empirical-1'
                        % case 2 : 采用估计修正 估计信噪比 估计l !真实情况
                        sigma2_estim = mean(obj.EigsHat(obj.k+1:end));
                        ell_estim = zeros(obj.k,1);
                        for l = 1:obj.k
                            lambda = obj.EigsHat(l)/sigma2_estim;
                            if lambda>=(1+sqrt(obj.c))^2
                                ell_estim(l) = (lambda-(1+obj.c))/2 + sqrt( (lambda-(1+obj.c))^2 - 4*obj.c)/2;
                            end
                        end
                        g = ((ell_estim).^(2)- obj.c)./((ell_estim).^(2) + obj.c .* (ell_estim).^(1));
                    case 'Empirical-2'
                        % case 2 : 采用估计修正 估计信噪比 估计l !真实情况
                        sigma2_estim = mean(obj.EigsHat(obj.k+1:end));
                        ell_estim = zeros(obj.k,1);
                        for l = 1:obj.k
                            lambda = obj.EigsHat(l)/sigma2_estim;
                            if lambda>=(1+sqrt(obj.c))^2
                                ell_estim(l) = (lambda-(1+obj.c))/2 + sqrt( (lambda-(1+obj.c))^2 - 4*obj.c)/2;
                            end
                        end
                        g = ((ell_estim).^(2)- obj.c)./((ell_estim).^(2) + obj.c .* (ell_estim).^(1));
                        mask1 = isinf(g);
                        g(mask1) = 1;
                    otherwise
                        g = (1- obj.c .* (obj.EigsTrue./obj.sigma2).^(-2))./(1 + obj.c .* (obj.EigsTrue./obj.sigma2).^(-1));

                end
                G = sqrt(1./(g*g.'));
                J_tmp = eye(obj.N);
                n = obj.N-1;
                J1 = J_tmp(1:n,:);
                J2 = J_tmp(2:end,:);
                Phi1 = obj.UsHat' * (J1' * J1) * obj.UsHat;
                Phi2 = obj.UsHat' * (J1' * J2) * obj.UsHat;

                Phi2_Repair = Phi2.*G;
                Phi =  inv(Phi1)* Phi2_Repair;
                [~,EigenValues] = eig(Phi,'vector');
                [DoA,index] = sort(angle(EigenValues));
                DoA = DoA.';
                EigenValues = EigenValues(index);
                MSE = sum((DoA - obj.ThetaTrue).^2) / obj.k;
%                 Bias =  sum(abs(DoA - obj.ThetaTrue))  / obj.k;
        end

        
        function [DoA,MSE,Bias] = GetMusic(obj,GridInterval,GridGap)
            % MUSIC可以和GMUSIC同时计算
            w_theta = linspace(GridInterval(1),GridInterval(2),GridGap);
            a = @(theta) exp(1i*theta*(0:obj.N-1)')/sqrt(obj.N);
            store_output = zeros(length(w_theta),1);
            for j = 1:length(w_theta)
                w_theta_i = w_theta(j);
                %MUSIC
                store_output(j) = (real(a(w_theta_i)'*obj.UsHat*(obj.UsHat')*a(w_theta_i)));
            end
            [pks,locs] = findpeaks(store_output);
            [~, index] = sort(pks,'descend');
            locs = locs(index);
            locs = locs(1:obj.k);
            DoA = sort(w_theta(locs),'ascend');
            MSE = sum((DoA - obj.ThetaTrue).^2) / obj.k;
            Bias =  sum(abs(DoA - obj.ThetaTrue))  / obj.k;
        end

        function [DoA,MSE,Bias] = GetGMusic(obj,GridInterval,GridGap)
            % MUSIC可以和GMUSIC同时计算
            w_theta = linspace(GridInterval(1),GridInterval(2),GridGap);
            a = @(theta) exp(1i*theta*(0:obj.N-1)')/sqrt(obj.N);
            store_output = zeros(length(w_theta),1);
            for j = 1:length(w_theta)
                %GMUSIC 方法
                theta = w_theta(j);
                sigma2_estim = mean(obj.EigsHat(obj.k+1:end));
                DGMUSIC = zeros(obj.k,obj.k);
                for l = 1:obj.k
                    lambda = obj.EigsHat(l)/sigma2_estim;
                    if lambda>=(1+sqrt(obj.c))^2
                        ell_estim = (lambda-(1+obj.c))/2 + sqrt( (lambda-(1+obj.c))^2 - 4*obj.c)/2;
                        DGMUSIC(l,l) = (ell_estim^2+obj.c*ell_estim)/(ell_estim^2-obj.c);
                    end
                end
                store_output(j,1) = (real(( a(theta)'*obj.UsHat*DGMUSIC*(obj.UsHat')*a(theta))));
            end
            [pks,locs] = findpeaks(store_output);
            [~, index] = sort(pks,'descend');
            locs = locs(index);
            locs = locs(1:obj.k);
            DoA = sort(w_theta(locs),'ascend');
            MSE = sum((DoA - obj.ThetaTrue).^2) / obj.k;
            Bias =  sum(abs(DoA - obj.ThetaTrue))  / obj.k;
        end

        function GenerateGuass(obj)
            % 此函数为了多次重复实验
            % 生成S
            obj.S = sqrtm(obj.P)*sqrt(1/2) *(randn(obj.k,obj.T) + 1i *randn(obj.k,obj.T));
            % 生成Z noise 
            obj.Z = sqrt(obj.sigma2/2) * (randn(obj.N,obj.T) + 1i* randn(obj.N,obj.T));

            % 构造协方差矩阵  生成新的子空间
            X = obj.A*obj.S + obj.Z;
            SCM = X*(X')/obj.T;
            [U,eigs_SCM] = eig(SCM,'vector');
            [eigs_SCM, index] = sort(eigs_SCM,'descend');
            U = U(:, index);
            obj.UsHat= U(:,1:obj.k);
            obj.EigsHat = eigs_SCM;
        end


        function CRB = GetCRB(obj)
            CRB = obj.sigma2 / (2*obj.T) *inv(real(obj.D'*(eye(obj.N)-obj.A*inv(obj.A'*obj.A)*obj.A')*obj.D) .*eye(obj.k));
            
        end

        function [MSE,Var,Bias] = GetStatNum(obj,DOA_Nb,MSE_Nb)
            MSE     =   mean(MSE_Nb,1);
            DOA_E   =   mean(DOA_Nb,1);
            
            Var   =   var(DOA_Nb,0,1);
            Var   =   mean(Var);

            Bias    =  (DOA_E - obj.ThetaTrue).^2;
            Bias    =  mean(Bias);
        end
    end
end