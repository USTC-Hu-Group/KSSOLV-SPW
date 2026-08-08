%{
function HY = mtimesY(H,Y)
    % MTIMES  Overload multiplication operator for Ham class
    %    HX = H*X returns a wavefun corresponding to H*X.
    %
    %    See also Ham, Wavefun.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    
    %  -----------------------------------------------
    % mtimesY multiplication operator for Ham*GpaWavefun Class, return
    % GpaWavefun Class.
    % mtimesY gives the result of Ham, GpaWavefun multiplication in GPA subspace.
    % The coordinates transformation relation between two spaces is:
    %    vec = (a_1,...,a_n1)Y = (e_1,...,e_ng)MY = (e_1,...,e_ng)X
    % HY is the matrix/vetcor representation of H*Y on sub-space spanned by
    % GPA basis M.
    % It returns a vector the same size as M (dims of GPA subspace)
    %    M: Transformation matrix form original to GPA  i.e. G~ = G * M (size(M) = ng n1)
    %    H: Hamiltonian on original PW basis (size ng)
    %    Y: Wavefunctions on Symmetrized GPA PW basis (size n1 ~= ng/ngpa)
    % Edit by Liangyu Wang 2025/02/28
    
    ncol = ncols(Y);
    npw = numel(Y.idxnz);
    %npol = size(Y.psi,1)/npw;
    npol = 1;
    nspin = H.nspin;
    lspinorb = H.lspinorb;
    X = Gpa2WF(Y);
    M = Y.gmat;
    % X.psi = M*Y;
    % n1 = size(M,2);
    
    % Apply Laplacian
    KinX =repmat(H.gkin,npol,ncol).*X;
    
    % Apply total local potential
    if npol == 1
        Xr3d  = ifft3(X);
        if nspin == 1
            vXr3d = repmat(H.vtot(:),1,ncol).*Xr3d;
        elseif nspin == 2
            vtot = H.vtot{Y.ispin};
            vXr3d = repmat(vtot(:),1,ncol).*Xr3d;
        end
        VX    = fft3(vXr3d);
        VtotX = VX(Y.idxnz,:);
    elseif npol == 2
        Xup = Wavefun(Y.psi(1:npw,:),Y.n1,Y.n2,Y.n3,Y.idxnz);
        Xdw = Wavefun(Y.psi(npw+1:end,:),Y.n1,Y.n2,Y.n3,Y.idxnz);
        Xrup = ifft3(Xup);
        Xrdw = ifft3(Xdw);
        if iscell(H.vtot)
            sup = repmat(H.vtot{1}(:)+H.vtot{4}(:),1,ncol).*Xrup + ...
                repmat(H.vtot{2}(:)-1i*H.vtot{3}(:),1,ncol).*Xrdw;
            sdw = repmat(H.vtot{1}(:)-H.vtot{4}(:),1,ncol).*Xrdw + ...
                repmat(H.vtot{2}(:)+1i*H.vtot{3}(:),1,ncol).*Xrup;
        else
            sup = repmat(H.vtot(:),1,ncol).*Xrup;
            sdw = repmat(H.vtot(:),1,ncol).*Xrdw;
        end
        VXup = fft3(sup);
        VXdw = fft3(sdw);
        VtotX = zeros(npw*npol,ncol);
        VtotX(1:npw,:) = VXup.psi(Y.idxnz,:);
        VtotX(npw+1:end,:) = VXdw.psi(Y.idxnz,:);
    end
    % Apply nonlocal pseudopotential
    if npol == 1
        VnlX = H.vnlmat*(repmat(H.vnlsign,1,ncol).*(H.vnlmat'*X));
    elseif npol == 2
        becp_up = H.vnlmat'*Xup;
        becp_dw = H.vnlmat'*Xdw;
    
        if lspinorb
            psup = H.vnlsign{1}*becp_up + H.vnlsign{2}*becp_dw;
            psdw = H.vnlsign{3}*becp_up + H.vnlsign{4}*becp_dw;
        else
            % H.vnlsign{2}{3} = 0 without soc
            psup = H.vnlsign{1}*becp_up;
            psdw = H.vnlsign{4}*becp_dw;
        end
        VnlX = zeros(npw*npol,ncol);
        VnlX(1:npw,:) = H.vnlmat*psup;
        VnlX(npw+1:end,:) = H.vnlmat*psdw;
    end
    
    mHY = M'*(KinX + VtotX + VnlX);
    HY = Y;
    HY.gpsi = mHY; 
    % Apply Fock exchange operator for hybrid functional
    if H.ishybrid
        if ~isempty(Y.ik) && Y.ik~=H.ik
            error('k points of H and X are different! %d %d',Y.ik,H.ik);
        end
        Y.ik = H.ik;
        Y.wks= H.wks;
        Vexx = H.vexx;
        HY = HY + Vexx(Y);
    end
    
    end
    %}
    function HY = mtimesY(H, Y)
    % MTIMESY  Subspace multiplication operator for Ham * SpwWavefun
    % This refactored version maps Y to full space, computes standard H*X,
    % and projects the result back, matching the basic physics implementation.
    
    global t_spw_forward;
    global t_spw_backward;

    t_start = tic;
    X = Gpa2WF(Y); 
    t_spw_forward = t_spw_forward + toc(t_start);

    X.ik = Y.ik;
    X.wks = Y.wks;
    X.ispin = Y.ispin; 
    
    HX = H * X; 
    
    M = Y.gmat;
    HY = Y; 

    t_start = tic;
    HY.gpsi = M' * HX.psi; 
    t_spw_backward = t_spw_backward + toc(t_start);
    
end