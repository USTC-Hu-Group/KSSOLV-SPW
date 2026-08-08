function [gen_cell, relations] = find_minimal_generators(R_cell)
% Input: R_cell - including 3x3 matrix in cell format.
% Output: 
%   gen_cell - the minimal generator cell 
%   relations - 1xN cell, the generation relation. 

    N = numel(R_cell);
    tol = 1e-6;
    
    all_traces = cellfun(@(x) trace(x), R_cell);
    all_dets = cellfun(@(x) det(x), R_cell);
    
    is_identity = (abs(all_traces - 3) < tol) & (abs(all_dets - 1) < tol);
    id_idx = find(is_identity, 1);
    if isempty(id_idx)
        error('No identity matrix in input.');
    end
    
    gen_cell = {R_cell{id_idx}}; 

    relations = cell(N, 1);
    relations{id_idx} = 1; 
    

    known_indices = id_idx;

    while numel(known_indices) < N
        missing = setdiff(1:N, known_indices);
        if isempty(missing), break; end
        
        new_gen_idx_in_R = missing(1);
        gen_cell{end+1} = R_cell{new_gen_idx_in_R};
        
        current_gen_num = numel(gen_cell);
        relations{new_gen_idx_in_R} = current_gen_num; 
        
        new_this_round = new_gen_idx_in_R;
        while ~isempty(new_this_round)
            next_round = [];
            for i = new_this_round(:).'
                for g = 2:current_gen_num
                    test_mat = R_cell{i} * gen_cell{g};
                    
                    t_tr = trace(test_mat);
                    t_det = det(test_mat);
                    candidates = find(abs(all_traces - t_tr) < tol & abs(all_dets - t_det) < tol);
                    
                    found_idx = [];
                    for k = candidates(:).'
                        if norm(R_cell{k} - test_mat, 'fro') < tol
                            found_idx = k;
                            break;
                        end
                    end
                    
                    if ~isempty(found_idx) && isempty(relations{found_idx})
                        if isequal(relations{i}, 1)
                            relations{found_idx} = g;
                        else
                            relations{found_idx} = [relations{i}, g];
                        end
                        next_round(end+1) = found_idx;
                    end
                end
            end 
            known_indices = unique([known_indices, next_round, new_this_round]);
            new_this_round = unique(next_round);
        end
    end
    gen_cell = gen_cell';
end