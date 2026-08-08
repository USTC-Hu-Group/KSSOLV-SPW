classdef Spwprec
    % GPAPREC KSSOLV class for GPA pre-transformation matrix for a certain k-point.
    %    Gprec = GpaprecSet(cry,sym_input_value) returns the transformation matrix P
    %    for Gamma point (the first k-point), normally using the highest order of transformation G.
    %
    %    Gprec = Gpaprec(cry,sym_input_value,n) returns the transformation matrix P 
    %    of the n_th group of k-points, using the same transformation G_n. 
    %
    %    The gpaprec class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %    nkpts          Number of k-points
    %    kgroup         Group index of ik
    %                   (k-points with same symmetry operation is stored in one group )
    %    kset           Index of k-points in kgroup               
    %    gmat           Eigenvector of G/Transformation matrix P
    %    dim            Dimension of each Irrep
    %    
    %    See also GpaprecSet.
    
    %  Copyright (c) 2025  Liangyu Wang,
    %                      University of Science and Technology of China
    %  This file is distributed under the terms of the MIT License.
    properties
        %--------------------------------%
        nkpts         % Number of k-points
        kgroup        % Group index of ik
        kset            % Index of k-point
        gmat          % Eigenvector of G/Transformation matrix
        dim           % dimension of each irrep
        % gnbnd         % Number of eigenval to calculate in subspace
    end
    
    methods
        
        function obj = Spwprec(cry,ik,sgidx,varargin)
            obj.nkpts = cry.nkpts;
            smat = genSPWs_allirep_k(cry,ik,sgidx,varargin{:});
            obj.gmat = smat; 
            dim = sum(~cellfun(@isempty,smat),2);
            obj.dim = dim;
        end

        % function GY = genGprecnbnd(GX,X)
        %     nbnd = genGpanbnd(GX,X)';
        %     GY = GX;
        %     GY.gnbnd = mat2cell(nbnd,ones(GY.gorder,1));
        % end
        % the func above is used to set the band to calculate without the
        % info from the first step, it is not a recommend.
   
    end
end
