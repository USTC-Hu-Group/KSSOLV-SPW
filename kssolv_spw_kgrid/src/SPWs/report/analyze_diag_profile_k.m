function [total_stats, step_stats, kpt_stats, irrep_stats_kpt] = analyze_diag_profile_k(prof_history)

    if ~iscell(prof_history)
        error('Input info.prof_history should be cell.');
    end

    [nscfiter, nk] = size(prof_history);

    fields = {'t_total', 't_hpsi', 't_proj', 't_eig', 't_upd', ...
              'n_hpsi_ops', 'n_restarts', 'dav_iter'};

    total_stats = struct();
    for f = 1:length(fields)
        total_stats.(fields{f}) = 0;
    end
    total_stats.restart_dims = [];

    empty_stat = total_stats;
    step_stats = repmat(empty_stat, nscfiter, 1);
    
    kpt_stats = repmat(empty_stat, nk, 1);

    irrep_stats_kpt = cell(nk, 1);
    for ik = 1:nk
        irrep_stats_kpt{ik} = repmat(empty_stat, 0, 1);
    end

    algo_type = 'Standard Global SCF';
    if ~isempty(prof_history) && iscell(prof_history{1,1})
        algo_type = 'SPW SCF (Irrep-Decoupled)';
    end

    for iscf = 1:nscfiter
        for ik = 1:nk
            item = prof_history{iscf, ik};
            
            if isstruct(item)
                items_to_process = {item}; 
            elseif iscell(item)
                items_to_process = item;   
            else
                continue; 
            end

            for idx = 1:length(items_to_process)
                pd = items_to_process{idx};
                if isempty(pd), continue; end

                if length(irrep_stats_kpt{ik}) < idx
                    num_padding = idx - length(irrep_stats_kpt{ik});
                    irrep_stats_kpt{ik} = [irrep_stats_kpt{ik}; repmat(empty_stat, num_padding, 1)]; 
                end

                for f = 1:length(fields)
                    fn = fields{f};
                    if isfield(pd, fn)
                        total_stats.(fn) = total_stats.(fn) + pd.(fn);
                        step_stats(iscf).(fn) = step_stats(iscf).(fn) + pd.(fn);
                        kpt_stats(ik).(fn) = kpt_stats(ik).(fn) + pd.(fn);
                                            
                        irrep_stats_kpt{ik}(idx).(fn) = irrep_stats_kpt{ik}(idx).(fn) + pd.(fn);
                    end
                end


                if isfield(pd, 'restart_dims') && ~isempty(pd.restart_dims)
                    total_stats.restart_dims = [total_stats.restart_dims, pd.restart_dims];
                    step_stats(iscf).restart_dims = [step_stats(iscf).restart_dims, pd.restart_dims];
                    kpt_stats(ik).restart_dims = [kpt_stats(ik).restart_dims, pd.restart_dims];
                    
                    irrep_stats_kpt{ik}(idx).restart_dims = [irrep_stats_kpt{ik}(idx).restart_dims, pd.restart_dims];
                end
            end
        end
    end

    fprintf('\n========================================================\n');
    fprintf('  Profiling Analysis \n');
    fprintf('  [%s]\n', algo_type);
    fprintf('========================================================\n');
    fprintf('SCF Iterations Evaluated : %d \n', nscfiter);
    fprintf('K-Points Evaluated       : %d \n', nk);
    fprintf('--------------------------------------------------------\n');
    fprintf('Total Diagonalization Time (t_total) : %10.4f s\n', total_stats.t_total);
    fprintf('  ├─ H*psi Time (t_hpsi)             : %10.4f s\n', total_stats.t_hpsi);
    fprintf('  ├─ Projection Time (t_proj)        : %10.4f s\n', total_stats.t_proj);
    fprintf('  ├─ Diagonalization Time (t_eig)    : %10.4f s\n', total_stats.t_eig);
    fprintf('  └─ Wavefun Update Time (t_upd)     : %10.4f s\n', total_stats.t_upd);
    
    overhead = total_stats.t_total - (total_stats.t_hpsi + total_stats.t_proj + total_stats.t_eig + total_stats.t_upd);
    fprintf('Framework Overhead (t_total - sum)   : %10.4f s\n', overhead);
    fprintf('--------------------------------------------------------\n');
    fprintf('Total Inner Iterations (dav_iter)    : %d \n', total_stats.dav_iter);
    fprintf('Total H*psi operations (n_hpsi_ops)  : %d \n', total_stats.n_hpsi_ops);
    fprintf('Total Subspace Restarts (n_restarts) : %d \n', total_stats.n_restarts);
    if total_stats.n_restarts > 0
        u_dims = unique(total_stats.restart_dims);
        dim_str = sprintf('%d, ', u_dims);
        fprintf('  └─ Restart dimension : [%s]\n', dim_str(1:end-2));
    end
    
 
    if nk > 0
        fprintf('\n====== Per-K-Point & Irrep Breakdown ======\n');
        for ik = 1:nk
            kstats = kpt_stats(ik);
            if kstats.dav_iter > 0 || kstats.t_total > 0
                fprintf(' [K-Point #%d] Total Evaluation\n', ik);
                fprintf('   ├─ Time  : %7.4f s (H*psi: %.4fs, proj: %.4fs, eig: %.4fs, upd: %.4fs)\n', ...
                    kstats.t_total, kstats.t_hpsi, kstats.t_proj, kstats.t_eig, kstats.t_upd);
                fprintf('   ├─ Ops   : %5d Iters, %5d H*psi ops\n', kstats.dav_iter, kstats.n_hpsi_ops);
                fprintf('   └─ Memory: %5d Restarts', kstats.n_restarts);
                if kstats.n_restarts > 0
                    u_dims = unique(kstats.restart_dims);
                    dim_str = sprintf('%d, ', u_dims);
                    fprintf(' (Hit limits at dim: [%s])\n', dim_str(1:end-2));
                else
                    fprintf(' (No memory limits hit)\n');
                end
                
                curr_irrep_stats = irrep_stats_kpt{ik};
                if length(curr_irrep_stats) > 1
                    fprintf('   --- Irrep Sub-Blocks for K-Point #%d ---\n', ik);
                    for idx = 1:length(curr_irrep_stats)
                        istats = curr_irrep_stats(idx);
                        if istats.dav_iter > 0 || istats.t_total > 0
                            fprintf('      [Irrep #%d] Time: %6.4fs | Iters: %4d | H*psi: %4d', ...
                                idx, istats.t_total, istats.dav_iter, istats.n_hpsi_ops);
                            if istats.n_restarts > 0
                                u_dims = unique(istats.restart_dims);
                                dim_str = sprintf('%d, ', u_dims);
                                fprintf(' | Restarts: %d (dim: [%s])\n', istats.n_restarts, dim_str(1:end-2));
                            else
                                fprintf('\n');
                            end
                        end
                    end
                end
                fprintf(' ------------------------------------------------------\n');
            end
        end
    end
    fprintf('\n');
end