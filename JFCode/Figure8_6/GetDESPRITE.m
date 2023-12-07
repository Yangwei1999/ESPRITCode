function Res= GetDESPRITE(U_S)
%GETESPRITE 此处显示有关此函数的摘要
%   此处显示详细说明

    [N,k] = size(U_S);
    J_tmp = eye(N);
    n = N-1;
    J1 = J_tmp(1:n,:);
    J2 = J_tmp(2:n+1,:);
    Phi =  diag(U_S' * J1' * J2 * U_S);
%     Phi = zeros(1,k);
%     for ii = 1 : k
%         Phi(1,ii) = U_S(:,ii)' * J1' * J2 * U_S(:,ii);
%     end
%     Phi =  diag(Phi);
%     [~,EigenValues] = eig(Phi,'vector');
    Res = sort((angle(Phi))).';
end

