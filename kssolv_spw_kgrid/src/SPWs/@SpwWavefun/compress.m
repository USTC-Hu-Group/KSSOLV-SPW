function X = compress(X,idxnz)
% WAVEFUN/COMPRESS Compress function for wave function class
%    X = COMPRESS(X,idxnz) returns the compact wave function with the
%    non-zero index idxnz.
%
%    See also Wavefun.
% This step may not exist in GpaWavefun, but pay attention to iscompact.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% idxnz will not appear in SPWs subspace, so I keep the psi unchanged.
X.idxnz = idxnz;
X.psi = X.psi(idxnz,:);
X.iscompact = 1;

end