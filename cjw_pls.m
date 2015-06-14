function [ PRESS, RESS, c_rate ] = cjw_plscv( X, Y, DA )
%CJW_PLSCV Function to perform PLS cross-validation on a given input.
%
% Synopsis:
%  [ PRESS, RESS, Q2, c_rate ] = cjw_plscv( X, Y, 'DA' )
%  [ PRESS, RESS, Q2 ] = cjw_plscv( X, Y )
%
% Description:
%  Perform Partial Least Squares cross-validation on a given 
%  independent block (X) and a given dependent block (Y).
%
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
%  DA               (flag)     If exists, enable classification accuracy.
%
% Outputs:
%  PRESS  (vector)   (a-1)     Vector containing PRESS score for each LV. 
%  RESS   (vector)   (a-1)     Vector containing RESS score for each LV. 
%  c_rate (vector)   (a-1)     Vector containing classification rate score 
%                              for a model using h LVs.
%
% Examples: 
%  [ PRESS, RESS, Q2, c_rate ] = cjw_plscv( X, Y )
%
% Requirements:
%  cjw_pls
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

[I, K] = size(Y);

a = rank(X);

leaveout = eye(size(X,1));

RESS(1) = K*(I-1);


for h=1:a-1
    Yh = zeros(size(Y));
    for j=1:size(X,1)
        l = leaveout(:,j);
        this_x = X(~l,:);
        this_y = Y(~l,:);
        this_model = cjw_pls2(this_x, this_y, h);
        pred = [1 X(j,:)]*this_model(1).b_coeff_star;
        Yh(j,:) = pred;
        
        if exist('DA', 'var')
            [~,guess]=min(abs(pred-1),[],2);
            actual = find(Y(j,:));
            correct(j) = (actual==guess);
        end
    end
    c_rate(h) = sum(correct)/length(correct);
    PRESS(h) = norm(Y-Yh).^2;
    
    this_model = cjw_pls2(X, Y, h);
    pred =  [ones(size(X,1),1) X]*this_model(1).b_coeff_star;
    RESS(h+1) = norm(Y-pred).^2;
    
end

end
