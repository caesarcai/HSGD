clear all
addpath('PROPACK')
n = 2^15;                 % dimension
r = 10;                   % rank
sample_rate = 0.1;        % sample rate
m = round(sample_rate*n); 
alpha = 0.1;              % 10% outlier
k = round(alpha*m);

[M,ox,~] = generate_signal(n,r,m,'true','false');

K = randsample(m,k);

obs = ox(M);
mean_of_ox = mean(abs(ox));
obs(K) = obs(K) + 20*mean_of_ox*2*(rand(k,1)-0.5+1i*(rand(k,1)-0.5));

eta = 0.6;
gamma_init = max(min(1.5, (1-2*log(n)*r/n)*(1/alpha)),1.2);
gamma_decay = 0.95;
tol = 1e-6;
max_iter = 1000;
proj = true;
[x,err,timer] = HSGD(obs,n,r,M,alpha,eta,gamma_init,gamma_decay,proj,tol,max_iter);

%SNR = snr(ox,ox-x)
norm(ox-x)/norm(ox)
