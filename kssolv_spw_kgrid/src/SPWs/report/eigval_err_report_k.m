function eigval_err_report_k(cry,info_store)
    fprintf('%s\n',repelem('-',20));
    fprintf('Total Energy Error (in Hartree):\n');
    
    if length(info_store) >= 2 && ~isempty(info_store{2})
        fprintf('Extra Bands:\n');
        err_etot = abs(info_store{1}.Etot - info_store{2}.Etot);
        fprintf('  Absolute Error: %2.2e\n',err_etot);
        fprintf('  Absolute Error per Atom: %2.2e\n',err_etot/sum(cry.natoms));
    end
    
    if length(info_store) >= 3 && ~isempty(info_store{3})
        fprintf('SPWs:\n');
        err_etot = abs(info_store{1}.Etot - info_store{3}.Etot);
        fprintf('  Absolute Error: %2.2e\n',err_etot);
        fprintf('  Absolute Error per Atom: %2.2e\n',err_etot/sum(cry.natoms));
    end

    fprintf('%s\n',repelem('-',20));
    fprintf('Eigval Error Report (in Hartree):\n')
    
    % Identify number of bands and K-points from the cell array
    nbnd0 = size(info_store{1}.Eigvals{1},1);
    nk = length(info_store{1}.Eigvals);

    % e00_mat will be a [nbnd0 x nk] matrix, where each column is a K-point
    e00_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{1}.Eigvals,'UniformOutput',false).');
    
    if length(info_store) >= 2 && ~isempty(info_store{2})
        e10_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{2}.Eigvals,'UniformOutput',false).');
        dife_c = e10_mat - e00_mat;
        fprintf('Extra Bands (Overall):\n  Maximum Absolute Error: %.2e \n',max(abs(dife_c(:))));
    end

    if length(info_store) >= 3 && ~isempty(info_store{3})
        e11_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{3}.Eigvals,'UniformOutput',false).');
        
        % dife_spw is now size [nbnd0 x nk]
        dife_spw = e11_mat - e00_mat; 
        
        fprintf('SPWs (Overall):\n')
        dife_spw_all = dife_spw(:);
        
        fprintf('  Maximum Absolute Error: %2.2e \n', max(abs(dife_spw_all)));
        
        % Fixed normalizations: using mean() to average over all elements
        e_L1 = mean(abs(dife_spw_all));
        fprintf('  Mean Absolute Error(L1 Norm): %2.2e \n', e_L1);

        e_RMSE = sqrt(mean(dife_spw_all.^2));
        fprintf('  Root Mean Square Error: %2.2e \n', e_RMSE);

        e_MBE = mean(dife_spw_all);
        fprintf('  Mean Bias Error: %2.2e \n', e_MBE); % Fixed copy-paste bug here

        % =================================================================
        % Per-K-Point Eigval Error Breakdown
        % =================================================================
        if nk > 1
            fprintf('\n====== Per-K-Point Eigval Error Breakdown ======\n');
            for ik = 1:nk
                % Extract the specific column for the current K-point
                diff_k = dife_spw(:, ik);
                
                fprintf(' [K-Point #%d]\n', ik);
                fprintf('   ├─ Maximum Absolute Error : %2.2e\n', max(abs(diff_k)));
                fprintf('   ├─ Mean Absolute Error    : %2.2e\n', mean(abs(diff_k)));
                fprintf('   ├─ Root Mean Square Error : %2.2e\n', sqrt(mean(diff_k.^2)));
                fprintf('   └─ Mean Bias Error        : %2.2e\n', mean(diff_k));
            end
            fprintf('================================================\n');
        end
    end
end