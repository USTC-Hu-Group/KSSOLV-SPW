function [Psi_obj, lambada] = davidson_qe_spw(mol, H, Psi_obj, tol, maxit, iter, h_diag, s_diag)
% precondition comprises Ekin,Vloc and Vnlc
    %diagt = tic;
    %fprintf("Start davidson diagonalization, convergence threshold is set to %20.13e\n"...
        %,tol);
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
    % cgoptions.btype = ones(1,cgoptions.notcnv);
    
    irep = Psi_obj.irep;
    if  ~isempty(H.eband_spw)
        % 乘以 2 是为了适配原有 davidson 求解器中的 Rydberg/Hartree 单位转换约定
        e = 2 * H.eband_spw{irep}; 
    else
        error('H.eband_spw does not exist or does not contain data for block irep = %d', irep);
    end

    % if ~isa(mol, 'Crystal') && mol.lsda
    %     % 这里假定无论是 X 类还是 Y 类，都保留了 ispin 属性来标识自旋
    %     e = 2 * H.eband((1:mol.nbnd) + mol.nbnd * (Psi_obj.ispin - 1));
    % else
    %     e = 2 * H.eband;
    % end

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
        [Psi_obj, e, cgoptions] = diag_david_spw(H, Psi_obj, e, cgoptions);
        ntry = ntry + 1;
        if ntry > 5 || cgoptions.notcnv == 0
            conv = true;
        end
    end
    
    % output (已通过函数签名第一项输出 [Psi_obj, lambada])
    lambada = e / 2;
    
    fprintf('Total number of iteration is %d\n', cgoptions.dav_iter);
    %fprintf('Total time for diag = %20.4e\n',toc(diagt));
end