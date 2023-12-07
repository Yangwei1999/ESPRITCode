function res = GetWellDIS(a,b,N,T)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明  生产wellbull分布
    tes  = wblrnd(a,b,[N,T]);
    res = normalize(tes,2,'center','mean','scale','std');
end