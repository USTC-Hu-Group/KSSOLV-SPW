function Z = mtimes(X,Y)
% WAVEFUN/MTIMES Mtimes function for wave function class
%    M = MTIMES(X,Y) returns a matrix as the multiplication of two wave
%    functions.
%
%    Z = MTIMES(M,Y) or Z = MTIMES(X,M) returns a wave function as the
%    multiplication of a matrix or scalar with a wave function.   
%
%    Note: the left multiplication of a non-scalar matrix with a wave
%    function returns a regular matrix.  
%
%    See also Wavefun.
%    Y1*Y2 = Y; 
%    Y1*M1 = Y;
%    M1*Y1 = M;

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isa(X,'SpwWavefun') && isa(Y,'SpwWavefun')
    if X.trans
        Xpsi = X.gpsi';
    else
        Xpsi = X.gpsi;
    end
    if Y.trans
        Ypsi = Y.gpsi';
    else
        Ypsi = Y.gpsi;
    end
    Z = Xpsi*Ypsi;
    return;
end

if isa(X,'SpwWavefun')
    if X.trans
        Z = X.gpsi'*Y;
    else
        Z = X;
        Z.gpsi = X.gpsi*Y;
    end
    return;
end

if isa(Y,'SpwWavefun')
    if Y.trans
        Z = Y;
        Z.gpsi = Y.gpsi*X';
    else
        if isscalar(X)
            Z = Y;
            Z.gpsi = X * Y.gpsi;
        else
            Z = X*Y.gpsi;
        end
    end
    return;
end

end
