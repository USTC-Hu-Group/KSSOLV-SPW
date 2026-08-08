function [X0, H] = genY0_random_spw(mol, H, smat_cell)
% The random generation method of wavefunctions
% Add k-grid and GPA logic by Liangyu Wang
% Adapted for nspin == 1 (spinless) ONLY

ishybrid = H.ishybrid;
H.ishybrid = 0;

if isa( mol, 'Crystal' )
    % nXcols = mol.nbnd * ones(mol.nkpts, 1);
    nkpts = mol.nkpts;
    X0 = cell(nkpts, 1);
    idxnzcell = cell(nkpts,1);
    H.eband = cell(nkpts, 1);
else
    % nXcols = mol.nbnd;
    nkpts = 1;
end

n1   = mol.n1;
n2   = mol.n2;
n3   = mol.n3;
tpiba = 2*pi/mol.alat;
info.first = true;

for ik = 1:nkpts  
    grid = Ggrid(ik, mol);
    idxnz = grid.idxnz;
    idxnzcell{ik} = idxnz;
    % 获取当前 k 点对应的 GPA 精度控制结构
    Gprec = smat_cell{ik};
    gorder = size(Gprec.dim,1);
    Qcell = cell(gorder, 1);
    
    nspws = cellfun(@(x) size(x, 2), Gprec.gmat(:,1));
    nbnd_spw = mol.spw_nband{ik};
    
    if isa(mol, 'Crystal')
        gg = sumel(mol.kpts(ik,:).^2);
    else
        gg = 0;
    end
    
    % 遍历 GPA 阶数
    for ig = 1:gorder 
        npw = nspws(ig);     % Dimension of subspace G_ik
        if npw < 1
            continue;
        end
        M = Gprec.gmat{ig};
        % gridgkk = diag(M' * (M .* grid.gkk));
        gridgkk = (M.^2)' * grid.gkk;
        [gkk, id] = sort(gridgkk);
        
        % map = zeros(numel(id), 1);
        [~, map] = sort(id);
   
        psif = zeros(npw, nbnd_spw(ig));  
        
        for j = 1:nbnd_spw(ig)
            for ig1 = 1:npw
                [rr, info] = randy(info);
                [arg, info] = randy(info);
                arg = 2 * pi * arg;
                psif(ig1, j) = complex(rr * cos(arg), rr * sin(arg)) / ...
                               (1 + gkk(ig1) / tpiba^2 + gg / tpiba^2);
            end
        end
        Qcell{ig} = psif(map, :);
    end
    
    % Rotate wavefunctions
    if isa( mol, 'Crystal' )
        X0_ik = SpwWavefunSet(Qcell, n1, n2, n3, idxnzcell{ik}, mol.wks(ik));
        X0_ik = genGpaInfo(X0_ik, Gprec);
        H_eband_ik = cell(gorder, 1);
        
        for ig = 1:gorder
            % 在 genY0_random 原代码中屏蔽了 eband 返回值，这里保留接口规范
            [X0_ik{ig}, H_eband_ik{ig}] = rotate_wfc_spw(mol, H{ik}, X0_ik{ig});
        end   
        
        X0{ik} = X0_ik;
        H.eband_spw{ik} = H_eband_ik;
        idx_eband = repelem(1:gorder,Gprec.dim);
        H.eband{ik} = vertcat(H_eband_ik{idx_eband});
    else   
        % 对于分子/单 k 点的常规情况
        X0_ik = SpwWavefun(Qcell{1}, n1, n2, n3, idxnz);
        [X0, H.eband] = rotate_wfc_Gpa(mol, H, X0_ik); 
    end
end

H.ishybrid = ishybrid;
end


function [X_out, ev] = rotate_wfc_spw(mol, H, X)
% find the best wavefunctions from the subspace which
% is spaned by the initial random wavefunctions
%  -----------------------------------------------
% only change X.psi to X.gpsi (and HX to HY)
%  GPA part by Liangyu Wang 2025
HX = H*X;
hc = X'*HX;
sc = X'*X;
eigsopts.isreal = false;
eigsopts.maxit  = 300;
eigsopts.tol    = 1e-16;
[V,D,flag] = eigs(hc,sc,ncols(X),'SR',eigsopts);
d = real(diag(D));
[sd,id]=sort(d);
ev = sd;
V = X.gpsi*V(:,id);

if (flag~=0)
   fprintf('Convergence not reached in eigs!, pause...\');
   pause;
end
X_out = SpwWavefun(V,mol.n1,mol.n2,mol.n3,H.idxnz,X.ispin);
end