function [ model ] = cjw_pls(X, Y, a)
%CJW_PLS Function to perform PLS on a given input
%
% Synopsis:
%  [ model ] = cjw_pls(X, Y, a, verb)
%
% Description:
%  Perform Partial Least Squares Regression, using a factors, on a given 
%  independent block (X) and a given dependent block (Y).
%
% Inputs:
%  X    (matrix)    (n x m)    Matrix of independant variables.
%                              Each row of X is a sample.
%                              Each column of X is a variable.
%
%  Y    (matrix)    (n x p)    Matrix of dependant variables.
%                              Each row of Y is a sample.
%                              Each column of Y is a variable.
%
%  a    (int)       (scalar)   Number of factors (aka latent components)
%                              to use when performing PLS.
%
%  verb (-)         (-)        Flag which, if exists, activates verbose
%                              mode.
%
% Outputs:
%  model            (struct)   Structure containing the model of the PLS.  
%
% Examples:
%  [ model ] = cjw_pls(X, Y, 5, 'v')
%
% Requirements:
%  None
% 
% References:
%  [a] Geladi, P. (1986). Partial least-squares regression: a tutorial. 
%      Analytica Chimica Acta, 185, 1?17. doi:10.1016/0003-2670(86)80028-9
%  [b] Abdi, H. (2003). Partial Least Squares ( PLS ) Regression. In 
%      Encyclopedia of Social Sciences Research Methods (pp. 1-7). 
%      doi:10.1097/EDE.0b013e3181d74bf5
%  [c] Rosipal, R., & Krämer, N. (2006). Overview and Recent Advances in 
%      Partial Least Squares. Subspace Latent Structure and Feature 
%      Selection, 3940, 34-51. doi:10.1007/11752790_2
%  [d] Pedersen, M. S., et al. (2008). The Matrix Cookbook. 
%      Matrix (Vol. M, pp. 1?71). doi:10.1111/j.1365-294X.2006.03161.x
%  
% Author:
%  Chris Winchurch <cjw111(at)imperial.ac.uk>
%
% License:
%  The program is free for non-commercial academic use. Please 
%  contact the authors if you are interested in using the software
%  for commercial purposes. The software must not modified or
%  re-distributed without prior permission of the author.
%
% Notes on implementation:
%  -Algorithm is implemented based on Geladi & Kowalski's PLS tutorial[1].
%  -The parenthesised numbers refer to the steps for the full PLS algorithm
%   as written in the appendix of [1].
%
% Changes:
%  21/02/2015  First edition
%  11/06/2015  Added verbose mode
%
%%

% Normalise vector x to have magnitude of 1.
normalise = @(x) x./norm(x);

v=exist('verb', 'var');
   
%% Pre-processing and Checks
% Deduce variables from inputs
n = size(X,1);
m = size(X,2);
p = size(Y,2);

% Check inputs
assert(size(X,1)==size(Y,1), ...
    ['Independant feature matrix X has %d samples, '  ...
     'dependant feature matrix Y has %d samples. '    ...
     'Number of samples should be equal for X & Y.'], ...
     size(X,1), size(Y,1));
assert(a<=size(Y,1), ...
    ['Dependant feature matrix Y has %d samples, ' ...
     'number of factors (a) specified as %d. '    ...
     'Number of factors must be less than samples in Y a<=size(Y,1).'], ...
     size(Y,1), a);

% Summarise inputs
if v
fprintf('Summary of inputs:\n');
fprintf(['Independant feature matrix X (%dx%d) ' ...
        'has %d samples, %d dimensions.\n'], n, m, n, m);
fprintf(['Dependant feature matrix Y (%dx%d) '   ...
        'has %d samples, %d dimensions.\n'], n, p, n, p);
fprintf('%d factors will be used in the PLS regression.\n', a);
end
% Retain scaling information about independant and dependant matrices
mu_X = mean(X);
mu_Y = mean(Y);
sigma_X = std(X);
sigma_Y = std(Y);

% Standardise independant and dependant variable matrices
X = zscore(X);
Y = zscore(Y);

%% Algorithm

% Settings 
tol = eps;
max_iter = 100;

E{1} = X;
F{1} = Y;

model = struct();

% For each component:
for h = 1:a
    y_idx = 1;
    u = Y(:,y_idx)/norm(Y(:,y_idx));           % [a](1) (with u normalised)
    t=u;
    
    converged = false;
    iter = 0;
    while ~converged
        
        % In the X block
        w = normalise(E{h}.'*u);                                   % [b](1)
        t = normalise(E{h}*w);                                     % [b](2)

        % In the Y block
        if p>1 % Only perform if Y block has more than one variable ...
            q = normalise(F{h}.'*t);                               % [b](3)
            u = F{h}*q;                                            % [b](4)
        else   % ... else set to q=1.
            q=1;
            u = F{h}*q;                                            % [b](4)
            break;
        end
        
        % Only check for convergence after 2nd iteration.
        if iter>0 
            % Check convergence                                      [a](8)
            converged = norm(t-t_last) < tol || iter > max_iter;
        end
        t_last = t;
        iter = iter + 1;
        
    end 

    % Compute the factor loadings for X                               % [b]
    pp = E{h}.'*t;

    % Compute the value of b
    b = t.'*u;                                                        % [b]
    
    % Calculate residual for use in next iteration:
    % (E corresponds to X, F corresponds to Y)
    E{h+1} = E{h} - t*pp.';
    F{h+1} = F{h} - b*t*q.';
    
    % Calculate R^2 (R-squared) for this factor
    R2X = (t'*t)*(pp'*pp)/(sum(sum(X.^2)));
    R2Y = (b^2)*(t'*t)*(q'*q)/(sum(sum(Y.^2)));
   
    % Assign the variables into the return model struct array:
    model(h).E = E{h};
    model(h).F = F{h};
    model(h).p = pp;
    model(h).q = q;
    model(h).w = w;
    model(h).t = t;
    model(h).u = u;
    model(h).b = b;
    model(h).R2X = R2X;
    model(h).R2Y = R2Y;

end

[model(1:a).B] = deal(diag([model.b]));

% Calcuate the regression coefficients:                      [b](Section 3)
% b_coeff = (P^(T+)BQ') (where P^(T+) is  Moore-Penrose pseudo-inverse) 
%         = (W*(P'*W)^(-1))BQ'                                   [d](3.6.3)
b_coeff = [model.w]*([model.p].'*[model.w])^(-1)*[model(1).B]*[model.q]';

% Scale the coefficients according to the original scaling of the data:
% Create scale matrix for input mean and variance:
tmp = bsxfun(@rdivide, [-mu_X; eye(m)], sigma_X);
% Use matrix to scale coefficients:
tmp2 = tmp*b_coeff;
% Scale the coefficients according to output variance:
tmp3 = bsxfun(@times, tmp2, sigma_Y);
% Scale the coefficients according to output mean:
b_coeff_star = tmp3 + [mu_Y; zeros(m,p)];

% Assign the variables into the return model struct array:
model(1).b_coeff = b_coeff;
model(1).b_coeff_star = b_coeff_star;

end
