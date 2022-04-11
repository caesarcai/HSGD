function y = fhmvmultiply(fft_h,x,lh,L,ind1,ind2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fast (multi-level) Hankel matrix vector mulitplication.
%
%Inputs
% fft_h: h generates the multi-level Hankel matrix H, fft_h is the 1D fft
%        of h of length L. 
%     x: column vector to be left multiplied by H.
%    lh: length of h.
%     L: the next number of power 2 greater than or equal to lh.
%  ind1: indeces to pad reversed x.
%  ind2: indeces to extract the vector after multiplication.
%
%Output
%     y: H*x.  
%
%Reference: Lu L, Xu W, Qiao S. A fast SVD for multilevel block Hankel
%matrices with minimal memory storage[J]. Numerical Algorithms, 2015, 
%69(4): 875-891.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lx = length(x);

xrev = x(lx:-1:1); % reverse x

xx = zeros(lh,1);

xx(ind1) = xrev; % put reversed x in the correct locations

fft_xx = fft(xx,L);

yy = ifft(fft_h.*fft_xx);

y = yy(ind2); % extract the vector after mulitplication 