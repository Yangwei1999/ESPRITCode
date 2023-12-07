function MUSIC_estim = GetMusic(U_S)

    [N,k] = size(U_S);
    w_theta = linspace(-pi/2,pi/2,10000);
    a = @(theta) exp(1i*theta*(0:N-1)')/sqrt(N);
    store_output = zeros(length(w_theta),1);
    for j = 1:length(w_theta)
        w_theta_i = w_theta(j);
        %MUSIC
        store_output(j) = (real(a(w_theta_i)'*U_S*(U_S')*a(w_theta_i)));
    end

    [pks,locs] = findpeaks(store_output);
    [~, index] = sort(pks,'descend');
    locs = locs(index);
    locs = locs(1:k);
    MUSIC_estim = sort(w_theta(locs));
end

