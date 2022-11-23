function mu = cadzow_inco_partial_obs(x_K,K,n,r)
% cadzow denoising for 1D signal

x=zeros(n,1);
x(K)=x_K;

%n = size(x,1);
if mod(n,2)
    p = (n+1)/2;
    DD = [1:p p-1:-1:1]';
else
    p = n/2;
    DD = [1:p p p-1:-1:1]';
end
q = n+1-p;

c_s = max(n/p,n/q);

opts = []; opts.eta = 1e-15;
xx = x;

Yforward = @(z) fhmvmultiply_1D(xx,z);
Ytranspose = @(z) fhmvmultiply_1D(conj(xx),z);

try
    [U,SS,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
ss = diag(SS(1:r,1:r)); %sigma_D = ss(1);
%U = U(:,1:r);
%V = V(:,1:r);

for i = 1:r
    ui = U(:,i);
    vi = V(:,i);
    x = x+ss(i)*conv_fft(ui,conj(vi));
end
x = x./DD;

Yforward = @(y) fhmvmultiply_1D(x,y);
Ytranspose = @(y) fhmvmultiply_1D(conj(x),y);
try
    [U,SS,V] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
catch
    fprintf('SVD did not converge.\n');
    return;
end
%ss = diag(SS(1:r,1:r)); 

%sigma_L = ss(1);

%U = U(:,1:r);
%V = V(:,1:r);

row_norm = zeros(p,1);
for i = 1:p
    row_norm(i) = norm(U(i,:))^2;
end

col_norm = zeros(q,1);
for j = 1:q
    col_norm(j) = norm(V(j,:))^2;
end

mu = max(max(row_norm),max(col_norm))*n/(c_s*r);

end