function [Psi_obj, lambada, prof_data] = davidson_qe_spw_profile(mol, H, Psi_obj, tol, maxit, iter, h_diag, s_diag)
% precondition comprises Ekin,Vloc and Vnlc
    npol = 1 + mol.noncolin;
    
    % ================= 新增：动态识别波函数字段名 =================
    if isa(Psi_obj,'Wavefun')|| isfield(Psi_obj, 'psi')
        fname = 'psi';
    elseif isa(Psi_obj, 'SpwWavefun') || isfield(Psi_obj, 'gpsi')
        fname = 'gpsi';
    else
        error('Unsupported object type: Must have either ''psi'' or ''gpsi'' property.');
    end
    % ==========================================================

    % set parameters 
    cgoptions = struct;
    cgoptions.tol = tol;
    cgoptions.dav_iter = 0;
    cgoptions.maxiter = maxit;
    cgoptions.h_diag = h_diag;
    cgoptions.s_diag = s_diag;
    
    % 使用动态字段名获取未收敛波函数的初始数量
    cgoptions.notcnv = size(Psi_obj.(fname), 2); 
    cgoptions.npol = npol;
    
    % ==== [探查器: 每次调用全新初始化] ====
    cgoptions.prof_data = struct(...
        't_total', 0, ...      % 当次对角化外壳总耗时
        't_hpsi', 0, ...       % 当次 H*X 操作总时间
        't_proj', 0, ...       % 当次子空间投影与构建 Hc, Sc 总时间
        't_eig', 0, ...        % 当次稠密小矩阵对角化 diag_h 总时间
        't_upd', 0, ...        % 当次波函数更新、残差计算及预条件时间
        'n_hpsi_ops', 0, ...   % 当次等效 H*psi 调用的单带向量总数
        'n_restarts', 0,...
        'restart_dims', []);      % 当次触发最大子空间维度导致重启的次数
    % ======================================
    
    irep = Psi_obj.irep;
    if  ~isempty(H.eband_spw)
        % 乘以 2 是为了适配原有 davidson 求解器中的 Rydberg/Hartree 单位转换约定
        e = 2 * H.eband_spw{irep}; 
    else
        error('H.eband_spw does not exist or does not contain data for block irep = %d', irep);
    end

    if isempty(mol.efermi)
        cgoptions.btype = ones(1,cgoptions.notcnv);
    else
        cgoptions.btype = (e < 2*(mol.efermi + 5e-3));
    end
    
    conv = false;
    ntry = 0;
    while ~conv
        cgoptions.lrot = (iter == 1);
        % 调用已适配好的 diag_david，传入多态对象 Psi_obj
        [Psi_obj, e, cgoptions] = diag_david_spw_profile(H, Psi_obj, e, cgoptions);
        ntry = ntry + 1;
        if ntry > 5 || cgoptions.notcnv == 0
            conv = true;
        end
    end
    
    % output
    lambada = e / 2;
    
    % ==== [探查器: 提取探查数据并返回] ====
    prof_data = cgoptions.prof_data;
    prof_data.dav_iter = cgoptions.dav_iter; % 将迭代步数也存入结构体备用
    % ======================================
    
    fprintf('Total number of iteration is %d\n', cgoptions.dav_iter);
end