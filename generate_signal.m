function [K,ox,f] = generate_signal(N,r,m,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate simulated (1D, 2D, 3D) spectrally sparse signals, with or without 
%separation between the frequencies, with or without damping.
%
%Inputs
%        N: a vector specify lengths of the signal in different dimensions, 
%           the number of dimensions <= 3.           
%        r: model order of the signal.
%        m: number or percentage of observed entries. If less than 1, it's
%           considered as a percentage.
% varargin: (optional) first specify whether to generate frequencies with 
%           separations. Input should be without (false,'false',0) or with 
%           (default, true,'true',1). Also specify whether the signal is 
%           damped. Either undamped (default, false,'false',0) or damped 
%           (true,'true',1).
%
%Outputs
%        K: indices for the observed entries.
%       ox: generated signal as a vector.
%        f: frequencies in different dimensions, stored column by column.
%
%< Example 1 > 
%[K,ox,f] = generate_signal(128,5,40) will generate a 1D undamped signal of
%length 128 with model order 5, with separation between frequencies and it
%is observed at 40 locations.  
%< Example 2 >
%[K,ox,f] = generate_signal([128 64 32],10,0.04,false,true) will generate a
%3D damped signal of size 128*64*32 with model order 10, without separation 
%between the frequencies and we have 4 percent samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = nargin;

if nargs < 3
    error('Not enough input parameters!');
end

switch nargs
    case 3
        separation = true;
        damping = false;
    case 4
        separation = varargin{1};
        damping = false;
    case 5
        separation = varargin{1};
        damping = varargin{2};
    otherwise
        error('Too many input parameters!')
end

% Generate complex coefficients of the signal
dynamic_range = 10;
c = exp(1i*2*pi*rand(r,1)).*(1+10.^(rand(r,1).*(dynamic_range/20)));

dim = numel(N); % dimension of the signal

switch dim
    
    case 1
        
        N1 = N; s1 = 0:N1-1;
        
        switch separation     
            case {false,'false',0} % without separation  
                f1 = rand(r,1);        
            case {true,'true',1} % with separation
                d1 = 1.5/N1;
                E1 = 1-r*d1; % excess space
                if E1 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f1 = rand(r+1,1); % wrap around distance
                    fs = E1*f1(1:r)/sum(f1);
                    fs = d1*ones(r,1)+fs;
                    f1 = cumsum(fs);
                end
            otherwise
                error('Separation should either be true or false!')
        end
        
        switch damping    
            case {false,'false',0} % without damping
                ox = exp(s1'*(1i*2*pi*f1'))*c;              
            case {true,'true',1} % with damping
                d1 = round(N1/4)+round(N1/4)*rand(r,1);
                d1 = -1./d1;
                ox = exp(s1'*(d1'+1i*2*pi*f1'))*c;
            otherwise
                error('Damping should either be true or false!')
        end
        
        f = f1;
        
        if m < 1
            m = round(N1*m);
        end
        
        K = randsample(N1,m); % sample without replacement
        
    case 2

        N1 = N(1); s1 = 0:N1-1;
        N2 = N(2); s2 = 0:N2-1;

        switch separation
            case {false,'false',0} % without separation 
                f1 = rand(r,1);
                f2 = rand(r,1);
            case {true,'true',1} % with separation
                d1 = 1.5/N1;
                E1 = 1-r*d1; % excess space
                if E1 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f1 = rand(r+1,1); % wrap around distance
                    fs = E1*f1(1:r)/sum(f1);
                    fs = d1*ones(r,1)+fs;
                    f1 = cumsum(fs);
                end
                d2 = 1.5/N2;
                E2 = 1-r*d2; % excess space
                if E2 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f2 = rand(r+1,1); % wrap around distance
                    fs = E2*f2(1:r)/sum(f2);
                    fs = d2*ones(r,1)+fs;
                    f2 = cumsum(fs);
                end
            otherwise
                error('Separation should either be true or false!')
        end
        
        switch damping        
            case {false,'false',0} % without damping
                ox = exp(kron(ones(N2,1),s1')*(1i*2*pi*f1')+kron(s2',ones(N1,1))...
                    *(1i*2*pi*f2'))*c;    
            case {true,'true',1} % with damping
                d1 = round(N1/4)+round(N1/4)*rand(r,1);
                d1 = -1./d1;
                d2 = round(N2/4)+round(N2/4)*rand(r,1);
                d2 = -1./d2;
                ox = exp(kron(ones(N2,1),s1')*(d1'+1i*2*pi*f1')+kron(s2',ones(N1,1))...
                    *(d2'+1i*2*pi*f2'))*c;
            otherwise
                error('Damping should either be true or false!')
        end
        
        f = [f1 f2];
        
        if m < 1
            m = round(N1*N2*m);
        end
        
        K = randsample(N1*N2,m); % sample without replacement
        
    case 3

        N1 = N(1); s1 = 0:N1-1;
        N2 = N(2); s2 = 0:N2-1;
        N3 = N(3); s3 = 0:N3-1;

        switch separation
            case {false,'false',0} % without separation 
                f1 = rand(r,1);
                f2 = rand(r,1);
                f3 = rand(r,1);
            case {true,'true',1} % with separation
                d1 = 1.5/N1;
                E1 = 1-r*d1; % excess space
                if E1 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f1 = rand(r+1,1); % wrap around distance
                    fs = E1*f1(1:r)/sum(f1);
                    fs = d1*ones(r,1)+fs;
                    f1 = cumsum(fs);
                end
                d2 = 1.5/N2;
                E2 = 1-r*d2; % excess space
                if E2 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f2 = rand(r+1,1); % wrap around distance
                    fs = E2*f2(1:r)/sum(f2);
                    fs = d2*ones(r,1)+fs;
                    f2 = cumsum(fs);
                end
                d3 = 1.5/N3;
                E3 = 1-r*d3; % excess space
                if E3 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f3 = rand(r+1,1); % wrap around distance
                    fs = E3*f3(1:r)/sum(f3);
                    fs = d3*ones(r,1)+fs;
                    f3 = cumsum(fs);
                end
            otherwise
                error('Separation should either be true or false!')
        end
        
        switch damping       
            case {false,'false',0} % without damping
                ox = exp(kron(ones(N3,1),kron(ones(N2,1),s1'))*(1i*2*pi*f1')...
                    +kron(ones(N3,1),kron(s2',ones(N1,1)))*(1i*2*pi*f2')...
                    +kron(s3',kron(ones(N2,1),ones(N1,1)))*(1i*2*pi*f3'))*c;     
            case {true,'true',1} % with damping
                d1 = round(N1/4)+round(N1/4)*rand(r,1);
                d1 = -1./d1;
                d2 = round(N2/4)+round(N2/4)*rand(r,1);
                d2 = -1./d2;
                d3 = round(N3/4)+round(N3/4)*rand(r,1);
                d3 = -1./d3;
                ox = exp(kron(ones(N3,1),kron(ones(N2,1),s1'))*(d1'+1i*2*pi*f1')...
                    +kron(ones(N3,1),kron(s2',ones(N1,1)))*(d2'+1i*2*pi*f2')...
                    +kron(s3',kron(ones(N2,1),ones(N1,1)))*(d3'+1i*2*pi*f3'))*c;
            otherwise
                error('Damping should either be true or false!')
        end

        f = [f1 f2 f3];

        if m < 1
            m = round(N1*N2*N3*m);
        end
        
        K = randsample(N1*N2*N3,m); % sample without replacement
        
    otherwise
        
        error('Dimension of the signal should be <= 3!')
       
end