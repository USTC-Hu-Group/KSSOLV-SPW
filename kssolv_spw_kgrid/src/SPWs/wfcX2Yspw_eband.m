function [Y, eband_spw] = wfcX2Yspw_eband(spwprec, H, X)
% wfcX2Y transform Wavefun Class X to SpwWavefunSet Y and outputs the 
% corresponding energy distributions.
%
% Outputs:
%   Y         - SpwWavefunSet class containing the projected wavefunctions
%   eband_spw - Cell array (nirep x 1) storing the energies of the new basis
% wfcX2Y transform Wavefun Class X to GpaWavefunSet Y or cell of array Y0.
% It returns GpaWavefunSet Y with Y{i} being its subset.
% Y{i} is the coordinate of GpaWavefun in subspace Gi, which is projected
% and orthogonalized from X. 
% Here H & X should be that of 1 k-point, not BH or BX.
%-------------------------------------------------------------------------%
% Mind: Because of band degeneracy, if the highest band is degenerate and
% its invariant subspace in not fully included, there will be some extra
% bands distributed equally on GPA subspaces. 
% So there will be several extra bands after projection than the input Wavefun X. 
% The treatment of projection and orthogonalization   
%-------------------------------------------------------------------------%
% see block_orthogonalize & analyze_array for details about orth-steps.

% Created by Liangyu Wang, 2025
% University of Science and Technology of China
    if ~isa(spwprec,'Spwprec')
        error('Input should be Spwprec Class, use Spwprec(cry,idx,P).');
    end

    smat = spwprec.gmat(:,1);
    nirep = size(smat,1);
    
    Y0 = cell(nirep,1);
    eband_spw = cell(nirep,1); 
    
    tol_vec = 1e-1;
    tol_deg = 1e-3;
    tol_mgs = 2e-1;

    if isa(X,'Wavefun') && isa(H,'Ham')
        for i = 1:nirep
            if isempty(smat{i,1})
                continue;
            end
            
            Y_temp = smat{i,1}' * X.psi;
            MMX = smat{i,1} * Y_temp;
            
            normvecc = vecnorm(MMX);
            idx = normvecc > tol_vec;

            Y0{i} = Y_temp(:, idx) ./ normvecc(:, idx);
            
            ebnd = H.eband(idx);
            
            [unique_values, ~, e_index] = analyze_array(ebnd, tol_deg);
            
            [Y0{i}, eband_spw{i}] = block_orthogonalize(Y0{i}, e_index, unique_values,tol_mgs);
        end
    else
        error('X should be Wavefun & H should be Ham');
    end
    
    % transform Y0 cell to SpwWavefunSet class
    Y = SpwWavefunSet(Y0, X.n1, X.n2, X.n3, X.idxnz);
    Y.gmat = spwprec.gmat;
    Y.dim = spwprec.dim;

end

function [orthogonalized_matrix, block_energies] = block_orthogonalize(input_matrix, e_index, unique_values,tol_mgs)
% block_orth: compute the orthogonalization of input column vectors block by block.
% Outputs both the orthogonalized matrix and the updated energy vector for this irep.

    orthogonalized_matrix = [];
    block_energies = []; 
    
    for i = 1:length(e_index)
        current_idx = e_index{i};
        current_energy = unique_values(i); 
        
        if length(current_idx) < 2
            orthogonalized_block = input_matrix(:, current_idx);
        else
            
            orthogonalized_block = gram_schmidt_mgs(input_matrix(:, current_idx),tol_mgs);
        end
        
        num_survived = size(orthogonalized_block, 2);
        
        if num_survived > 0
            orthogonalized_matrix = [orthogonalized_matrix, orthogonalized_block];
            block_energies = [block_energies; repmat(current_energy, num_survived, 1)];
        end
    end
end

function independent_vectors = gram_schmidt_mgs(Y,tol_mgs)
% Modified Gram-Schmidt (MGS) method to extract linearly independent column vectors.
% MGS is numerically much more stable than Classical Gram-Schmidt.

    [n, m] = size(Y);
    independent_vectors = zeros(n, m);  
    count = 0;  

    for i = 1:m
        v = Y(:, i);  

        for j = 1:count
            q = independent_vectors(:, j);
            v = v - q * (q' * v);
        end
        
        v_norm = norm(v);
        % If the vector is not zero, add it to the set of linearly independent vectors
        if v_norm > tol_mgs
            count = count + 1;  
            independent_vectors(:, count) = v / v_norm;  % Normalize and store
        end
    end

    % Trim to the effective linearly independent vectors
    independent_vectors = independent_vectors(:, 1:count);
end

function [unique_values, counts, indices] = analyze_array(arr,tol_deg)

    tol = tol_deg;
    n = length(arr);  

    unique_values = [];  
    counts = [];         
    indices = {};        

    if n == 0
        return; 
    end

    current_sum = arr(1);  
    count = 1;                
    index_list = 1;          

    for i = 2:n

        if abs(arr(i) - arr(i-1)) < tol
            count = count + 1;  
            index_list(end + 1) = i;  
            current_sum = current_sum + arr(i); 
        else
            unique_values = [unique_values; current_sum / count]; 
            counts = [counts; count];
            indices{end + 1} = index_list;  

            count = 1;  
            index_list = i;  
            current_sum = arr(i);
        end
    end

    unique_values = [unique_values; current_sum / count];
    counts = [counts; count];
    indices{end + 1} = index_list;
end