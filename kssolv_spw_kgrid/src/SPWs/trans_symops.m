function Xtc_transform = trans_symops(Xtc, P, p)

    if size(p,1) < 3
        p = p.';
    end
    M = [P, p; 0, 0, 0, 1];
    invM = inv(M);


    Xtc_transform = cellfun(@(D) invM * D * M, Xtc, 'UniformOutput', false);


    Xtc_transform = filter_rot_trans(Xtc_transform);
end

function sym_out = filter_rot_trans(sym_input)

    tol = 1e-5;

    cleaned_ops = cell(size(sym_input));
    for i = 1:length(sym_input)
        R = round(sym_input{i}(1:3, 1:3)); 
        t = sym_input{i}(1:3, 4);
        

        t = mod(t, 1);
        t(abs(t) < tol) = 0;          
        t(abs(t - 1) < tol) = 0;      
        
        cleaned_ops{i} = [R, t; 0, 0, 0, 1];
    end

    num_ops = length(cleaned_ops);
    flat_ops = zeros(num_ops, 16);
    for i = 1:num_ops
        flat_ops(i, :) = reshape(cleaned_ops{i}, 1, 16);
    end

    [~, ia, ~] = unique(round(flat_ops, 6), 'rows', 'stable');
    
    sym_out = cleaned_ops(ia);
end
