function Z = minus(X,Y)
% WAVEFUN/MINUS Minus function for wave function class
%    Z = MINUS(X,Y) returns the wave function of X - Y.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isa(X,'SpwWavefun') && isa(Y,'SpwWavefun')
    Z = X;
    Z.gpsi = X.gpsi - Y.gpsi;
    return;
end

if isa(X,'SpwWavefun')
    Z = X;
    Z.gpsi = X.gpsi - Y;
    return;
end

if isa(Y,'SpwWavefun')
    Z = Y;
    Z.gpsi = X - Y.gpsi;
    return;
end

end