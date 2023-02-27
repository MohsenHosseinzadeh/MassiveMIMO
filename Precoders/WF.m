function [x, beta, P] = WF(s, H, N0)
% =========================================================================
% Wiener filter (WF) precoder
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 precoded vector
%       - beta: precoding factor (scalar)
%       - P: BxU precoding matrix
% =========================================================================

    % number of UEs
    [U, ~] = size(H);
    
    % precoding matrix (before normalization)
    T = H' / (H*H' + U*N0*eye(U));

    % precoding factor
    beta = sqrt(real(trace((T*T'))));
    
    % precoding matrix
    P = 1/beta*T;
    
    % precoded vector
    x = P*s;

end
