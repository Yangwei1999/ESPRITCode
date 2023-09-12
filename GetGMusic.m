function MUSIC_estim = GetGMusic(U_S,eigs_SCM,c)
%GETGMUSIC 此处显示有关此函数的摘要
%   此处显示详细说明

    [N,k] = size(U_S);
    w_theta = linspace(-pi,pi,10000);
    a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
 
    store_output = zeros(length(w_theta),1);

    for j = 1:length(w_theta)
        theta = w_theta(j);
        sigma2_estim = mean(eigs_SCM(k+1:end));
        D = zeros(k,k);
        for l = 1:k
            lambda = eigs_SCM(l)/sigma2_estim;
            if lambda>=(1+sqrt(c))^2
                ell_estim = (lambda-(1+c))/2 + sqrt( (lambda-(1+c))^2 - 4*c)/2;
                D(l,l) = (ell_estim^2+c*ell_estim)/(ell_estim^2-c);
            else
                D(l,l) = 1;
            end
        end
        store_output(j,1) = (real(( a(theta)'*U_S*D*(U_S')*a(theta))));
    end


    [pks,locs] = findpeaks(store_output);
    [~, index] = sort(pks,'descend');
    locs = locs(index);
    locs = locs(1:k);
    MUSIC_estim = sort(w_theta(locs));

end

