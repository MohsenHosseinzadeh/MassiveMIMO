function [x, beta] = EXS(s, H, N0)
% =========================================================================
% exhaustive search (EXS)
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
    
    % 1-bit alphabet and quantizer
    alphabet = [-1-1i; 1-1i; -1+1i; 1+1i] / sqrt(2*B);
    quantizer = @(z) (sign(real(z)) + 1i*sign(imag(z)))/sqrt(2*B);

    % augmented channel matrix
    Hb = [H; sqrt(U*N0)*eye(B)];
 
    % QR decomposition
    [~,R] = qr(Hb); order = 1:B;

    % MRT vector (scaled) before/after quantization
    zMRT = H(:,order)'*s;
    xMRT = quantizer(zMRT);

    % initialization
    PA = zeros(B,1); % path
    ST = zeros(B,length(alphabet)); % stack
    radius = inf;
    
    % preprocessing
    num_next = nan(4,B); % numerator (next step)
    num_potential = cell(1,B); % numerator (potential steps)
    den_present = nan(4,B); % denominator (present)
    den_future = nan(1,B); % denominator (future)
    for l = 1:B
        num_potential{:,l} = R(1:l-1,l)*alphabet.';
        num_next(:,l) = R(l,l)*alphabet;
        den_present(:,l) = real(zMRT(l)'*alphabet);
        den_future(l) = real(zMRT(1:l-1)'*xMRT(1:l-1));
    end
    
    % root node
    level = B;
    
    % numerator (lower bound)
    numerator = abs(num_next(:,l)).^2; % present only
    
    % denominator (upper bound)
    denominator = (abs(den_present(:,level)) + abs(den_future(level))).^2; % present and future
    
    % add root node to stack
    ST(level,:) = numerator./denominator;
    
    % sphere decoding
    while level <= B
        
        % find smallest PED in boundary
        [minPED,idx] = min(ST(level,:));
        
        % proceed only if list is not empty
        if minPED<inf
            
            ST(level,idx) = inf; % mark child as tested
            NewPath = [idx; PA(level+1:end,1)]; % new best path
            
            % search child
            if minPED < radius
                
                % valid candidate found
                if level>1 
                    
                    % expand this best node
                    PA(level:end,1) = NewPath;
                    level = level - 1; % downstep
                    
                    % numerator (lower bound)
                    num_past = norm(R(level+1:end,level+1:end)*alphabet(PA(level+1:end,1)),2)^2;
                    num_present = abs(num_next(:, level) + R(level,level+1:end) * alphabet(PA(level+1:end,1))).^2;
                    num_future = 0; % no future metric
                    numerator =  num_past + num_present + num_future;
                    
                    % denominator (upper bound)
                    den_past = real(zMRT(level+1:end,1)'*alphabet(PA(level+1:end,1)));
                    denominator = ( abs(den_present(:,level) +  den_past) +  abs(den_future(level)) ).^2;
                    
                    % add to stack
                    ST(level,:) = numerator./denominator;
                    
                else
                    
                    % valid leaf found
                    idxhat = NewPath;
                    x = alphabet(idxhat);
                    
                    % update radius (radius reduction)
                    radius = minPED;
                    
                end
                
            end
            
        else
            
            % no more child nodes to be checked
            level=level+1;
            
        end
        
    end
    
    % re-order vector and compute beta
    x(order,1) = x;
    beta = real(x'*H'*s)/(norm(H*x,2)^2+U*N0);
    
    % flip negative beta
    if beta < 0 
        x = -x;
        beta = -beta;
    end
    
end
