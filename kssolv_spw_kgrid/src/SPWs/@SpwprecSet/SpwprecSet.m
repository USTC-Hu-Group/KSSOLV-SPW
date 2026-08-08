classdef SpwprecSet 
    % SpwprecSet is the Symmetrized Plane-waves basis for each k-point.
    % However, the irreducible representation table for all k-points is
    % hard to generate. This code Only Works for Gamma point by now!
    %    The gpaprecset class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %    nkpts          Number of k-points
    %    kidx           Index of k-points
    %    spwpreccell    Cell of gpaprec matrix of each k-groups
    %   
    %    See also GpaprecSet.

% Created by Sheng Chen and Liangyu Wang, 2025
% University of Science and Technology of China
    properties
        nkpts         % number of k-points 
        spwpreccell = {}
    end
    
    methods
         
        function obj = SpwprecSet(cry,sgidx,varargin)
                nk = cry.nkpts;
                obj.nkpts = nk;
                obj.spwpreccell = cell(nk, 1);
                for ik = 1:nk
                    Gprec = Spwprec(cry,ik,sgidx,varargin{:});
                    Gprec.kset = ik;
                    obj.spwpreccell{ik} = Gprec;
                end
        end

        function X = subsref(BX,S)
                switch S(1).type
                    case '{}'
                        ik = S(1).subs{1};
                        if numel(ik) > 1
                            error('Wrong sub index.')
                        end
                        if numel(S) > 1
                            X = builtin('subsref',BX.spwpreccell{ik},S(2:end));
                        else
                            X = BX.spwpreccell{ik};           
                        end
                    otherwise
                        X = builtin('subsref',BX,S);
                end
       end
        
    end
end
