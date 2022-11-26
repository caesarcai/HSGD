function [si,x,err,timer] = HSGD(z_K,n,r,K,alpha,eta,gamma,gamma_decay_factot,proj,tol,max_iter)if mod(n,2)    p = (n+1)/2;  %not sample rate here    DD = [1:p p-1:-1:1].';else    p = n/2;    DD = [1:p p p-1:-1:1].';endq = n+1-p;norm_obs = norm(z_K);m = numel(z_K); % number of observed samplesstep_size = n/m;     % inversed sample rate% other variables to be pre-allocatedHV = zeros(p,r);HtU = zeros(q,r);timer = zeros(max_iter,1);err = zeros(max_iter,1);if proj    %% Estimate mu via one-step Cadzow    mu = cadzow_inco_partial_obs(z_K,K,n,r);    c_s = max(n/p,n/q);end%%%%%%%%%%%%%%%%%%%%For fast fft svd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%t1 = tic; [val,ind] = sort_frac(z_K,alpha);   s_K = zeros(m,1);s_K(ind) = val;x = zeros(n,1);x(K)= step_size.*(z_K-s_K);L = 2^nextpow2(n); % next power of 2 for faster fft% indeces for fhmvmultiply to useind1 = 1:q; ind2 = q:n; ind3 = 1:p; ind4 = p:n; Yforward = @(y) fhmvmultiply_1D(x,y);Ytranspose = @(y) fhmvmultiply_1D(conj(x),y);opts = []; opts.eta = 1e-16;[U,Sigma,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);eta = eta/Sigma(1,1); % constant stepsizesigma = sqrt(diag(Sigma)).';  U = bsxfun(@times,U,sigma);V = bsxfun(@times,V,sigma);if proj    incoh_U = (2*mu*r*c_s/n)*norm(U,2)^2;    incoh_V = (2*mu*r*c_s/n)*norm(V,2)^2;    [U,V] = proj_oper( U, V, incoh_U, incoh_V );endx = zeros(n,1);for i = 1:r    ui = U(:,i);    vi = conj(V(:,i));    ui = fft(ui,L);    vi = fft(vi,L);    ui = ui.*vi;    ui = ifft(ui);    ui = ui(1:n);    x = x+ui;endx = x./DD;     [val,ind] = sort_frac(z_K-x(K),gamma*alpha); s_K = zeros(m,1);s_K(ind) = val;    init_err = norm(x(K)+s_K-z_K)/norm_obs;init_timer = toc(t1);si = 0;% Successive Iterationsfor iter = 1:max_iter    tic;    UtU = U'*U;    VtV = V'*V;    x(K) = step_size*(z_K-s_K)+(1-step_size)*x(K);        fft_x = fft(x,L);       for i = 1:r        vi = V(:,i);        HV(:,i) = fhmvmultiply(fft_x,vi,n,L,ind1,ind2);    end            S = VtV;    Ug = HV-U*S;        fft_conj_x = fft(conj(x),L);      for i = 1:r        ui = U(:,i);        HtU(:,i) = fhmvmultiply(fft_conj_x,ui,n,L,ind3,ind4);    end    S = UtU;    Vg = HtU-V*S;        U = U+eta*Ug-1/16*eta*U*(U'*U-V'*V);    V = V+eta*Vg+1/16*eta*V*(U'*U-V'*V);    if proj        [U,V] = proj_oper( U, V, incoh_U, incoh_V );    end        x = zeros(n,1);    for i = 1:r        ui = U(:,i);        vi = conj(V(:,i));        ui = fft(ui,L);        vi = fft(vi,L);        ui = ui.*vi;        ui = ifft(ui);        ui = ui(1:n);        x = x+ui;    end    x = x./DD;        gamma = (gamma-1.05)*gamma_decay_factot+1.05;    [val,ind] = sort_frac(z_K-x(K),gamma*alpha);        s_K = zeros(m,1);    s_K(ind) = val;        err(iter) = norm(x(K)+s_K-z_K)/norm_obs;    timer(iter) = toc;        if err(iter) < tol      si = 1;      err(2:iter+1) = err(1:iter);      err(1) = init_err;      err = err(1:iter+1);      timer(2:iter+1) = timer(1:iter);      timer(1) = init_timer;      timer = timer(1:iter+1);      fprintf('Total %d iteration, final error: %e, total time without init: %f , with init: %f\n======================================\n', iter, err(iter), sum(timer(timer>0)),sum(timer(timer>0))+init_timer);      return    elseif err(iter) > 10 % blow up      err(2:iter+1) = err(1:iter);      err(1) = init_err;      err = err(1:iter+1);      timer(2:iter+1) = timer(1:iter);      timer(1) = init_timer;      timer = timer(1:iter+1);      fprintf('Total %d iteration, final error: %e, diverged\n======================================\n', iter, err(iter));      return    else       fprintf('Iteration %d: error: %e, timer: %f \n', iter, err(iter), timer(iter));    end endend % main function ends herefunction [val,ind] = sort_frac(obs,alpha)m = length(obs);k = round(alpha*m);[val,ind] = maxk(obs,k,'ComparisonMethod','abs');endfunction [U, V] = proj_oper( U, V, incoh_U, incoh_V )row_norm_square_U = sum(U.^2,2);big_rows_U = row_norm_square_U > incoh_U;U(big_rows_U,:) = bsxfun(@times,U(big_rows_U,:),(incoh_U ./ sqrt(row_norm_square_U(big_rows_U)))); row_norm_square_V = sum(V.^2,2);big_rows_V = row_norm_square_V > incoh_V;V(big_rows_V,:) = bsxfun(@times,V(big_rows_V,:),(incoh_V ./ sqrt(row_norm_square_V(big_rows_V))));  end