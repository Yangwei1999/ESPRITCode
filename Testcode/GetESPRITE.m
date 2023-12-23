function Res= GetESPRITE(U_S)
%GETESPRITE 此处显示有关此函数的摘要
%   此处显示详细说明

    [N,~] = size(U_S);
    J_tmp = eye(N);
%     n = round(N/2);
    n = N-1;
    J1 = J_tmp(1:n,:);
    J2 = J_tmp(2:end,:);
    Phi =  inv(U_S' * J1' * J1 * U_S) * (U_S' * J1' * J2 * U_S);
    [~,EigenValues] = eig(Phi,'vector');
    Res = sort(angle(EigenValues)).';
end

