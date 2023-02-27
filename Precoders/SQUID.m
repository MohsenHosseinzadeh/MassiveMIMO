function [x, beta] = SQUID(s,H,N0)
% =========================================================================
% squared infinity-norm relaxation with Douglas-Rachford splitting (SQUID)
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 precoded vector
%       - beta: precoding factor (scalar)
% =========================================================================

    % dimensions
    [U, B] = size(H);

    % number of iterations
    iter = 50; 
    
    % gain: affects the convergence of SQUID 
    % - set to 1 for large problems (e.g., 128 BS antennas) or low SNR
    % - set to small value for small problems or high SNR
    gain = 1; % default value (must be optimized)
    
    % relxation parameter: affects the convergence of SQUID 
    rho = 1; % default value (must be optimized)

    % convert to real-valued channel
    HR = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [ real(s) ; imag(s) ];
  
    % initialize
    b = zeros(2*B,1);
    c = zeros(2*B,1);
    
    % pre-processing
    Q = HR'/((0.5/gain)*eye(2*U) + HR*HR');
    sMF = HR'*sR;
    sREG = (2*gain)*(sMF - Q*(HR*sMF));
   
    for t=1:iter % SQUID loop
        z = 2*b -c;
        a = sREG + z - Q*(HR*z);
        b = prox_infinity_norm_squared(c+a-b,2*U*B*N0);        
        c = c + rho*(a - b);
    end
    
    % extract binary solution
    x = sign(b);     
    x = 1/sqrt(2*B)*(x(1:B,1)+1i*x(B+1:2*B,1));

    % compute beta
    Hx = H*x; 
    beta = real(Hx'*s)/(norm(Hx,2)^2+U*N0);
    
    % flip negative beta
    if beta < 0
        x = -x;
        beta = -beta;
    end

end

