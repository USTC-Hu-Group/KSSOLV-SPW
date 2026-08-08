classdef SpwWavefunSet
    % GpaWavefunSet KSSOLV class for Bloch gpa wave function
    % BY = genGpaInfo(BX,Gprec)
    % 
    %    Most functions is a copy of BlochWavefun Class.
    % 
    %    Remark: This BlochWavefun serves as a data container for wave
    %    functions. 
    %    It is the same with GpaWavefunSet and GpaWavefun.
    %    But BlochWavefun contains Wavefuns of all k-points, while GpaWavefunSet
    %    contains the GpaWavefuns of all gpa subspace for ONE certain
    %    k-point.
    %        
    %
    %    The Blochwavefun class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        nspin       Number of density components
    %        nkpts       Number of k-points
    %        wks         The weight for each k-potins
    %        wavefuncell Cell of wave functions
    %      ----------------------------GPA------------------------------
    %        gorder      Order of symmetry operation used by GPA
    %        geigval     Eigenvalue of G in this subspace
    %        gmat        Eigenvectors of G (GPA Basis set of this subspace)                     
    %        gconj       Label of complex conjugate eigenpairs 
    %        gpsi        Coordination of Wave functions on GPA basis set
    %        igpa        index of GPA blocks
    %    See also Atom, Wavefun, GpaWavefun, Gpaprec.
   
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    %  -----------------------------------------------
    %  GPA part by Liangyu Wang 2025
    properties (SetAccess = protected)
        ik 
        nspin = 1
        nkpts = 0
        wks
        spwwavefuncell = {}
    %----------GPA----------%
        dim
        gmat          % Eigenvector of G/Transformation matrix
    end
    methods
        function BX = SpwWavefunSet(varargin)
            switch (nargin)
                case 0
                    return;
                case 1
                    Xin = varargin{1};
                    BX.nkpts = numel(Xin);
                    BX.wks = ones(BX.nkpts,1);
                    BX.spwwavefuncell = cellfun(@SpwWavefun,Xin, ...
                        'UniformOutput',0);
                case 2
                    Xin = varargin{1};
                    BX.nkpts = numel(Xin);
                    BX.wks = wks;
                    BX.spwwavefuncell = cellfun(@SpwWavefun,Xin, ...
                        'UniformOutput',0);
                case 4
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    BX.nkpts = numel(Xin);
                    BX.wks = ones(BX.nkpts,1);
                    BX.spwwavefuncell = ...
                        cellfun(@(X)SpwWavefun(X,n1,n2,n3),Xin, ...
                        'UniformOutput',0);
                case 5
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    vec = varargin{5};
                    BX.nkpts = numel(Xin);
                    if numel(vec) == numel(Xin)
                        BX.wks = vec;
                        BX.spwwavefuncell = ...
                            cellfun(@(X)SpwWavefun(X,n1,n2,n3),Xin, ...
                            'UniformOutput',0);
                    else
                        BX.wks = ones(BX.nkpts,1);
                        BX.spwwavefuncell = ...
                            cellfun(@(X)SpwWavefun(X,n1,n2,n3,vec),Xin, ...
                            'UniformOutput',0);
                    end
                case 6
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    idxnz = varargin{5};
                    BX.wks = varargin{6};
                    BX.nkpts = numel(Xin);
                    BX.spwwavefuncell = ...
                        cellfun(@(X)SpwWavefun(X,n1,n2,n3,idxnz),Xin, ...
                        'UniformOutput',0);
                case 7
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    idxnz = varargin{5};
                    BX.wks = varargin{6};
                    BX.nspin = varargin{7};
                    if BX.nspin ~= 2
                        BX.nkpts = numel(Xin);
                        BX.spwwavefuncell = ...
                        cellfun(@(X)SpwWavefun(X,n1,n2,n3,idxnz),Xin, ...
                        'UniformOutput',0);
                    else
                        nkpts = numel(Xin)/2;
                        BX.nkpts = nkpts;
                        X_up = cell(nkpts,1);
                        X_dw = cell(nkpts,1);
                        for ik = 1:nkpts
                            X_up{ik} = Xin{ik};
                            X_dw{ik} = Xin{ik+nkpts};
                        end
                        wavefuncell_up = ...
                        cellfun(@(X)SpwWavefun(X,n1,n2,n3,idxnz,1),X_up, ...
                        'UniformOutput',0);
                        wavefuncell_dw = ...
                        cellfun(@(X)SpwWavefun(X,n1,n2,n3,idxnz,2),X_dw, ...
                        'UniformOutput',0);
                        BX.spwwavefuncell = [wavefuncell_up;wavefuncell_dw];
                    end
                otherwise
                    error('Wrong number of arguements');
            end
        end
        function BY = genGpaInfo(BX,Gprec)
            BY = BX;
            BY.gmat = Gprec.gmat;
            BY.dim = Gprec.dim;
        end
    end
end
