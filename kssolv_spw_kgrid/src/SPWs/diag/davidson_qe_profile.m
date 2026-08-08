function [X, lambada, prof_data] = davidson_qe_profile(mol, H, X0, tol, maxit, iter, h_diag, s_diag)
% precondition comprises Ekin,Vloc and Vnlc
    %diagt = tic;
    npol = 1 + mol.noncolin;
    
    % set parameters 
    cgoptions = struct;
    cgoptions.tol = tol;
    cgoptions.dav_iter=0;
    cgoptions.maxiter = maxit;
    cgoptions.h_diag = h_diag;
    cgoptions.s_diag = s_diag;
    cgoptions.notcnv = size(X0.psi,2);
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

    if ~isa(mol,'Crystal')&&mol.lsda
        e = 2*H.eband((1:mol.nbnd)+mol.nbnd*(X0.ispin-1));
    else
        e = 2*H.eband;
    end
    
    if isempty(mol.efermi)
        cgoptions.btype = ones(1,mol.nbnd);
    else
        cgoptions.btype = (e < 2*(mol.efermi + 5e-3));
    end
    
    conv = false;
    ntry = 0;
    while ~conv
        cgoptions.lrot = (iter == 1);
        [X0, e, cgoptions] = diag_david_profile(H,X0,e,cgoptions);
        ntry = ntry + 1;
        if ntry>5 || cgoptions.notcnv == 0
            conv = true;
        end
    end
    
    % output
    X = X0;
    lambada = e/2;
    
    % ==== [探查器: 提取探查数据并返回] ====
    prof_data = cgoptions.prof_data;
    prof_data.dav_iter = cgoptions.dav_iter; % 将迭代步数也存入结构体备用
    % ======================================
    
    fprintf('Total number of iteration is %d\n', cgoptions.dav_iter);
end