classdef Ggrid
    % GGRID KSSOLV class for the reciprocal space grid
    %    Grid = GGRID returns an empty frequency mask.
    %
    %    Grid = GGRID(mol) returns a frequency mask for the mol. 
    %    Shift is not added.
    %
    %    Grid = GGRID(mol, ecut) returns a frequency mask for the mol
    %    with the given energy cut. Shift is not added.
    %
    %    Grid = GGRID(ik, mol) or GGRID(ik, mol, ecut) returns a frequency 
    %    mask for the mol with a k-point shift specified by index ik.
    %
    %    Grid = GGRID(gx, gy, gz, ecut) returns a frequency mask for
    %    locations (gx, gy, gz) with the given energy cut. Shift is not added.
    %    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    %  Add k-grid by Liangyu Wang, USTC, 260323
    properties (SetAccess = protected)
        ecut
        ng
        idxnz
        gkk
        gkx
        gky
        gkz
        kkxyz
    end
    
    methods
        function Grid = Ggrid(varargin)
            
            %===========================================================
            % Input process
            
            % Return empty mask
            if nargin == 0
                return;
            end
            
            ik = 0;  
            mol = []; 
            
            % Check for signature: Ggrid(ik, mol, [ecut])  -> Backward compatibility for Ggridk
            if nargin >= 2 && isnumeric(varargin{1}) && isscalar(varargin{1}) && isa(varargin{2}, 'Molecule')
                ik = varargin{1};
                mol = varargin{2};
                if nargin >= 3 && isnumeric(varargin{3})
                    Gecut = varargin{3};
                else
                    Gecut = mol.ecut;
                end
                
            % Check for signature: Ggrid(mol, [ecut]) -> Standard Ggrid
            elseif isa(varargin{1}, 'Molecule')
                mol = varargin{1};
                if nargin >= 2 && isnumeric(varargin{2})
                    Gecut = varargin{2};
                else
                    Gecut = mol.ecut;
                end
                
            % Check for signature: Ggrid(gx, gy, gz, ecut) -> Location input
            elseif nargin >= 4 && isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                    && isnumeric(varargin{3}) && isnumeric(varargin{4})
                pregkx = varargin{1};
                pregky = varargin{2};
                pregkz = varargin{3};
                C  = eye(3);
                Gecut = varargin{4};
            else
                error('Invalid input.');
            end
            
            if ~isempty(mol)
                n1 = mol.n1;
                n2 = mol.n2;
                n3 = mol.n3;
                [I,J,K] = ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
                    (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
                    (0:n3-1)-((0:n3-1) >= n3/2)*n3);
                pregkx = I(:);
                pregky = J(:);
                pregkz = K(:);
                C  = mol.supercell;
            end
            
            %===========================================================
            % Construction of Ggrid
            %
            Grid.ecut = Gecut;
            xyz = [pregkx pregky pregkz];
            
            % NOTE: the transpose is very important.
            Creci = 2*pi*inv(C)';
            
            kkxyz = xyz*Creci;
            
            % Apply K-point shift if ik is specified and molecule is present
            if ik > 0 && ~isempty(mol)
                kkxyz_shift = kkxyz + mol.kpts(ik,:);
                kk = sum(kkxyz_shift.*kkxyz_shift, 2);
            else
                kk = sum(kkxyz.*kkxyz, 2);
            end
            
            Grid.idxnz = find(kk <= Grid.ecut * 2 * meDef());
                   
            Grid.ng    = length(Grid.idxnz);
            Grid.gkk   = kk(Grid.idxnz);
            Grid.gkx   = kkxyz(Grid.idxnz, 1);
            Grid.gky   = kkxyz(Grid.idxnz, 2);
            Grid.gkz   = kkxyz(Grid.idxnz, 3);
            
            Grid.kkxyz = kkxyz;
        end
    end
end