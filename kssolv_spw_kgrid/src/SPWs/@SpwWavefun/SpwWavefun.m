classdef SpwWavefun
    % GPAWAVEFUN KSSOLV class for GPA wave function
    %    Most functions is a copy of Wavefun Class.
    %
    %    X = Gpa2WF(Y) transform a  gap-wavefunction in subspace to a
    %    wavefuntion in the full space.
    %    Y = genGpaInfo(Y,X) delivers the transformation matrix from X to Y.   
    %      ----------------------------------------------------------
% case 0    Y = GpaWavefun()                 ?   empty Y   
%
% case 1    Y = GpaWavefun(psi)              F   Y = psi(n1,n2,n3)
%                                                ndims(psi) = 3 in (n1,n2,n3) format. 
%
% case 3    Y = GpaWavefun(n1,n2,n3)         F   empty Y 
%
% case 4-1  Y = GpaWavefun(n1,n2,n3,ncols)   ?   zeros(ng,ncols) Y 
%
% case 4-2  Y = GpaWavefun(gpsi,n1,n2,n3)    T   Y = gpsi
%
% case 5    Y = GpaWavefun(gpsi,n1,n2,n3,idxnz) T  Y = psi(n1,n2,n3)
%           maybe this form will be used in most cases.
%
% case 6    Y = GpaWavefun(gpsi,n1,n2,n3,idxnz,ispin)  F 
    %
    %    The wavefun class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        n1,n2,n3    Number of discretization points in each dimension
    %        nrows       Number of rows in the wave function
    %        ncols       Number of columns in the wave function
    %        idxnz       Indices for non-zero entries in the n1,n2,n3
    %        iscompact   Indicator for compact format storage
    %        trans       Indicator for transpose
    %        psi         2D data matrix storing each wave function in each
    %                    column
    %        occ         Occupation rate
    %      ----------------------------GPA------------------------------
    %        gorder      Order of symmetry operation used by GPA
    %        geigval     Eigenvalue of G in this subspace
    %        gmat        Eigenvectors of G (GPA Basis set of this subspace)                     
    %        gconj       Label of complex conjugate eigenpairs 
    %        gpsi        Coordination of Wave functions on GPA basis set
    %        igpa        index of GPA blocks
    %    See also Atom, Wavefun, Gpaprec.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    %  -----------------------------------------------
    %  GPA part by Liangyu Wang 2025
    properties (SetAccess = protected)
        n1 = 0
        n2 = 0
        n3 = 0
        idxnz
        iscompact = 0
        trans = 0
        psi
	ik
	wks
    %----------GPA----------%
    irep
    gmat          % Eigenvector of G/Transformation matrix
    gpsi          % <--Coordination of Wave functions on GPA basis set--> %

    %----------GPA----------%
    end
    properties (SetAccess = public)
        occ
        ispin = 1
    end
    methods
        function X = SpwWavefun(varargin)
            switch (nargin)
                case 0
                    return;
                case 1
                    Xin = varargin{1};
                    if iscell(Xin)
                        [X.n1,X.n2,X.n3] = size(Xin{1});
                        nrows = numel(Xin{1});
                        ncols = length(Xin);
                        X.gpsi = zeros(nrows,ncols);
                        for it = 1:length(Xin)
                            X.gpsi(:,it) = Xin{it}(:);
                        end
                    else
                        if ndims(Xin) == 3
                            [X.n1,X.n2,X.n3]=size(Xin);
                            X.gpsi  = Xin(:);
                        elseif ndims(Xin) == 4
                            [X.n1,X.n2,X.n3,ncols]=size(Xin);
                            nrows = X.n1*X.n2*X.n3;
                            X.gpsi  = reshape(Xin,nrows,ncols);
                        end
                    end
                case 3
                     X.n1  = varargin{1};
                     X.n2  = varargin{2};
                     X.n3  = varargin{3};
                case 4
                    if size(varargin{1},1) > 1
                    %if ndims(varargin{1}) > 1
                        X.n1  = varargin{2};
                        X.n2  = varargin{3};
                        X.n3  = varargin{4};
                        X.gpsi = varargin{1};
                    else
                        X.n1  = varargin{1};
                        X.n2  = varargin{2};
                        X.n3  = varargin{3};
                        nrows  = X.n1*X.n2*X.n3;
                        ncols  = varargin{4};
                        X.gpsi = zeros(nrows,ncols);
                    end
                case 5
                    Xin = varargin{1};
                    X.n1 = varargin{2};
                    X.n2 = varargin{3};
                    X.n3 = varargin{4};
                    X.idxnz = varargin{5};
                    X.iscompact = 1;
                    if iscell(Xin)
                        X.gpsi = zeros(numel(Xin{1}),length(Xin));
                        for it = 1:length(Xin)
                            X.gpsi(:,it) = Xin{it}(:);
                        end
                    else
                        X.gpsi  = Xin;
                    end
                case 6
                    Xin = varargin{1};
                    X.n1 = varargin{2};
                    X.n2 = varargin{3};
                    X.n3 = varargin{4};
                    X.idxnz = varargin{5};
                    X.ispin = varargin{6};
                    X.iscompact = 1;
                    if iscell(Xin)
                        X.gpsi = zeros(numel(Xin{1}),length(Xin));
                        for it = 1:length(Xin)
                            X.gpsi(:,it) = Xin{it}(:);
                        end
                    else
                        X.gpsi  = Xin;
                    end
                otherwise
                    error('Wrong number of arguments');
            end
        end
        function ind = end(X,k,n)
            szd = size(X.gpsi);
            if k < n
                ind = szd(k);
            else
                ind = prod(szd(k:end));
            end
            
        end

        function X = Gpa2WF(Y)
            X = Wavefun();
            X.n1 = Y.n1;
            X.n2 = Y.n2;
            X.n3 = Y.n3;
            X.idxnz = Y.idxnz;
            X.iscompact = Y.iscompact;
            X.trans = Y.trans;
            X.psi = Y.gmat*Y.gpsi;
            X.ispin = Y.ispin;
        end

        function Y = genGpaInfo(Y,X)
            % Here, X is GpaWavefun
            Y.gmat = X.gmat;
        end
    end
end
