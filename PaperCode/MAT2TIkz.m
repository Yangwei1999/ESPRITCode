%% matlab mat to tikz picture

[index1,index2] = size(MSE_VList);
for ii =  1 : index1
    for jj = 1: index2
        fprintf('(%d,%e)',(VariableList(jj)),MSE_VList(ii,jj))
    end
    disp('\n');
end

for jj = 1: length(CRB_Res)
    fprintf('(%d,%e)',(VariableList(jj)),CRB_Res(1,jj))
end
