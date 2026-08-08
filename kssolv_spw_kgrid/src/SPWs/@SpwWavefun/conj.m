function Xconj = conj(X)
% WAVEFUN/CONJ Conj function for wave function class
%    Xconj = CONJ(X) returns the conjugate value of the wave function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% I don't think that is the conjugte of the wavefunction, because conjugate
% will influence both the coefs of PWs and PWs itself, so here we can only
% say it is a conjugate of the coefs of a wavefunction. 
% 2025/12/01 by liangyu Wang
Xconj = X;
Xconj.gpsi = conj(X.gpsi);

end