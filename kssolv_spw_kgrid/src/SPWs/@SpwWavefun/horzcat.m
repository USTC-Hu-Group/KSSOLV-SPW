function Y = horzcat(varargin)
% WAVEFUN/HORZCAT Horzcat function for wave function class
%    Y = HORZCAT(X1,X2,...) returns the horizontal cat of wave functions.
%
%    Y = [X1,X2,...] returns the horizontal cat of wave functions. 
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

X = varargin{1};

nXcols = cellfun(@ncols,varargin);

if isa(varargin{1}.gpsi, 'gpuArray')&&isa(varargin{2}.gpsi, 'gpuArray')
    Ypsi = gpuArray.zeros(nrows(X),sum(nXcols));
else
    Ypsi = zeros(nrows(X),sum(nXcols));
end

idx = 0;
for it = 1:length(varargin)
    idx = idx(end)+(1:nXcols(it));
    Ypsi(:,idx) = varargin{it}.gpsi;
end
if X.iscompact
    Y = SpwWavefun(Ypsi,X.n1,X.n2,X.n3,X.idxnz);
else
    Y = SpwWavefun(Ypsi,X.n1,X.n2,X.n3);
end
Y.trans = X.trans;
Y.ispin = X.ispin;

% because horzcat do not get the gmat value, we need to add the info:
Y.gmat = X.gmat;
end
