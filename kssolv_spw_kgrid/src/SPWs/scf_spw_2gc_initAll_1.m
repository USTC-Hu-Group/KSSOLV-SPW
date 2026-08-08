function [mol,H,X,info] = scf_spw_2gc_initAll_1(mol,options,varargin)
% Self Consistent Field iteration (Modified for Two-Grid SPW method)
% (HSE/Hybrid functional logic has been removed for simplicity)
% The SPW-SCF code.
% This code use the Symmetrized Plane waves method to speedup the
% diagonalization part of SCF for Large System with Symmetry.
% A pilot SCF calculation will be performed first with reduced cutoff and
% low convergence limit to get a good prediction of band distribution
% {n_\alpha} for each Irrep \alpha.
% Then the SPW-SCF will start with the band distribution, solving the
% eigenvalue problem of each subspace of Irrep independently.

% Liangyu Wang 2026
% University of Science and Technology of China

if (nargin < 2)
    options = setksopt();
end

global verbose;
verbose = options.verbose;


if ~isfield(options, 'usespw')
    options.usespw = 0;
end

% ==========================================
%  (Two-Grid Method) Pilot Calculation to SPW
% ==========================================
if options.coarse
    fprintf('======================================\n');
    fprintf(' [Two-Grid SCF] Step 1: Coarse Grid \n');
    fprintf('======================================\n');
    
    % 1. Set crystal input for reduced cufoff SC
    mol_coarse = set(mol, 'ecut', mol.ecut_coarse);
    
    % 2. Setting of reduced cufoff SCF
    opt_coarse = options;
    opt_coarse.usespw = 0; 
    if isfield(options, 'scftol_coarse')
        opt_coarse.scftol = options.scftol_coarse;
    end
    opt_coarse.cgtol = options.cgtol_c;
    opt_coarse.maxcgiter = options.maxcgiter_c;
    opt_coarse.fixcgtol = 0;
    opt_coarse.test_proj = 1;
    opt_coarse.cbreak = 1;
    % 3. Run reduced cufoff regular SCF
    [mol_coarse, H_c, X_c, info_c, opt_coarse] = scf0(mol_coarse, opt_coarse,varargin{:});
    
    fprintf('======================================\n');
    fprintf(' [Two-Grid SCF] Step 2: Interpolation \n');
    fprintf('======================================\n');
    
    % 4. Reduced cutoff -> Target cutoff
    
    % X_fine   = interp_X_coarse2fine(X_c, mol_coarse, mol);

    % rho_fine = interp_rho_coarse2fine(H_c.rho, mol_coarse, mol);
    % rho_fine = getcharge(mol,X_fine,info_c.occ);
    fprintf('X and Rho are Not used!\n')
    if options.usespw
        S = [];

        isSym = cellfun(@(x) isa(x,'Sym'),varargin);
        
        if isfield(options,'symcharge') && options.symcharge ==1
            fprintf('Use Symmetrized Charge Densify.\n')
            idxSym = find(isSym,1);
            if ~isempty(idxSym)
                S = varargin{idxSym};
            else
                S = Sym(mol_coarse);
            end
        end

        if ~isempty(S.bilbao)
            infos = S.bilbao;
            smat_cell = SpwprecSet(mol_coarse,infos{:});
        else
            error('Need syminfo for spw method!')
        end
        nk = mol_coarse.nkpts;
        spw_nband = cell(nk,1);

        for ik = 1:mol_coarse.nkpts
            [Y0,~] = wfcX2Yspw_eband(smat_cell{ik},H_c{ik},X_c{ik});
            ycols = ncols(Y0).'; 
            fprintf('------------------Coarse Proj--------------------\n')
            fprintf('Coarse Ecut = %2.2f Ha\n',mol_coarse.ecut);
            fprintf('Total calculated bands for ik %02d = %03d\n',ik,sum(ycols));
            fprintf('Band Distribution :');  
            disp(ycols);
            spw_nband{ik} = ycols;
            clear Y
        end
        clear smat_cell
        mol.spw_nband = spw_nband;
    end

    fprintf('======================================\n');
    fprintf(' [Two-Grid SCF] Step 3: Fine Grid \n');
    fprintf('======================================\n');
    
    % 5. Do not heritage the Info from reduced SCF calculation.
    % options.rho0 = rho_fine;
    options.rho0 = [];
    options.X0   = [];
   
    
    % 6. Run the Target cutoff SCF 
    [mol, H, X, info, options] = scf0(mol, options,varargin{:});
    
    info.info_c = info_c;
% ==========================================
% Regular SCF
% ==========================================
else
    fprintf('Regular SCF for Pure DFT\n');
    [mol, H, X, info, options] = scf0(mol, options,varargin{:});
end
fprintf('======================================\n');
fprintf(' [Two-Grid SCF]Calculation Details: \n');
fprintf('======================================\n');
fprintf('Eigmethod: %s \n',options.eigmethod);
fprintf('maxcgiter: %02d \n',options.maxcgiter);
fprintf('cgtol: %.2e \n',options.cgtol);
fprintf('------------Fine Grid---------------\n');
fprintf(' Cutoff: %2.2f Ha\n',mol.ecut);
fprintf(' scftol: %.2e \n',options.scftol);
if options.coarse
    fprintf('-----------Coarse Grid--------------\n');
    fprintf(' Cutoff: %2.2f Ha\n',mol.ecut_coarse);
    fprintf(' scftol: %.2e \n',options.scftol_coarse);
    fprintf(' cgtol: %.2e \n',options.cgtol_c);
    fprintf(' maxcgiter: %02d \n',options.maxcgiter_c);
end

fprintf('======================================\n');
fprintf(' [Two-Grid SCF]Eigenvalue Problem Total Time Report: \n');
fprintf('======================================\n');

if options.usespw
    fprintf('Coarse Grid Wavefun Initial Time Use(s): %f \n',sum(info.info_c.t_direig));
    fprintf('Spws-Subspace Eigenvalue Time Use(s): %f \n',sum(info.t_spweig));
    fprintf('Spws-Subspace Projection time use(s): %f \n ',info.t_y2x);
else
    if options.coarse
        fprintf('Coarse Grid Wavefun Initial Time Use(s): %f \n',sum(info.info_c.t_direig));
    end
    fprintf('Regular Eigenvalue Time Use(s): %f\n',sum(info.t_direig));
end


end

function [mol,H,X,info,options] = scf0(mol,options,varargin)
% SCF Self Consistent Field iteration, both for semiconductor and metal.
%    [mol,H,X,info] = SCF(mol,options) adopts Self Consistent Field (SCF)
%    iteration to find the ground state minimum total energy and the
%    corresponding wave functions. mol is a Molecule object and options is
%    the options for running the SCF. Please read setksopt for detailed
%    information about options. SCF returns the molecule mol with/without
%    force, the Hamiltonian H, the wave functions X, and the information
%    for each iteration in info.

%This file is a merged version of the old scf, scf4m and scf4c.

global t_spw_forward;
global t_spw_backward;
t_spw_forward = 0;
t_spw_backward = 0;


fprintf('Total Number of Bands: %03d.\n',mol.nbnd);

if (nargin < 2)
    options = setksopt();
end

S = [];

isSym = cellfun(@(x) isa(x,'Sym'),varargin);

if isfield(options,'symcharge') && options.symcharge ==1
    fprintf('Use Symmetrized Charge Densify.\n')
    idxSym = find(isSym,1);
    if ~isempty(idxSym)
        S = varargin{idxSym};
    else
        S = Sym(mol);
    end
end

if (isfield(options,'usespw') && options.usespw ==1) || options.test_proj
    fprintf('Use Symmetrized Plane-Wave Method.\n')
    smat_cell = [];
    if ~isempty(S.bilbao)
        infos = S.bilbao;
        smat_cell = SpwprecSet(mol,infos{:});
    else
        error('Need syminfo for spw method!')
    end
end

scfstart  = tic;

% Initialize input variables
global verbose;
force      = options.force;
maxscfiter = options.maxscfiter;
scftol     = options.scftol;
what2mix   = options.what2mix;
mixtype    = options.mixtype;
mixdim     = options.mixdim;
betamix    = options.betamix;
brank      = options.brank;
X          = options.X0;
rho        = options.rho0;

iscryst    = isa(mol,'Crystal');
ishybrid   = options.ishybrid;
nspin      = mol.nspin;
domag      = mol.domag;
smear      = mol.smear;
Tbeta      = mol.temperature*8.6173324e-5/13.6;
%Tbeta      = 315774.67 / mol.temperature;
%options.cgtol = 1e-9;
%options.cgtol = 1e-2;
%The original value of cgtol is 1e-2, it seems not correct, and the final result is wrong.
%Anyway, 1e-9 is safe.

% Initialize Hamiltonian, Wavefun, and Preconditioners
if ~options.usespw 
    [mol,H,X,Hprec,nocc] = iterinit(mol,rho,X);
    nXcols     = ncols(X);
    ev         = zeros(sumel(nXcols),1);
else
    [mol,H,YS,HprecYS,nocc] = iterinit_spw(mol,smat_cell,rho,X);
    X = BlochWavefun();
    % H.eband = options.ev0;
end
% calculate Ewald and Ealphat
Eewald     = getEewald(mol);
Ealphat    = getEalphat(mol);

vion       = H.vion;
vext       = H.vext;
vtot       = H.vtot;
rho        = H.rho;

% Initialize output variables
Etotvec    = zeros(maxscfiter,1);
scferr     = zeros(maxscfiter,1);
dfmat      = [];
dvmat      = [];
cdfmat     = [];

if iscryst
	nkpts      = mol.nkpts;
	wks        = mol.wks;
    info.t_direig = [];
    if options.usespw
        nXcols = zeros(mol.nkpts,1);
        fprintf('======================================\n');
        fprintf(' [Two-Grid SCF] Step 4: Subspace Wavefun Projection \n');
        fprintf('======================================\n');
        info.t_yinit = 0;
        info.t_y2x = 0;
        info.t_spweig = []; 
        fprintf('------------Subspace Info------------\n');
        for ik = 1:nkpts
            YS{ik}.ik = ik; 
            ycols = ncols(YS{ik});
            fprintf('Total calculated bands = %03d\n',sum(ycols));
            nXcols(ik) = ycols'*YS{ik}.dim;
            fprintf('Band Distrubution ik %02d:\n',ik);
            disp(ycols');
            
        end
        ev = zeros(sumel(nXcols),1);
     end
end

if ishybrid
	Vexx       = options.vexx;
	H.vexx     = Vexx;
end

vhart = getVhart(mol,rho);
Ecoul_old = getEcoul(mol,rho,vhart);

scf_prof_history = cell(maxscfiter, nkpts);

fprintf('Beging SCF calculation for %s...\n',mol.name);
for iterscf = 1:maxscfiter
    
    fprintf('SCF iter %3d:\n', iterscf);
    options.iterscf = iterscf;
    
    rhoin  = rho;
    vtotin = vtot;
    
    first = (iterscf == 1);
    if first&&isfield(options,'ev0')
        H.eband = options.ev0;
    end
    if ~options.usespw
        if iscryst 
        	idx = 0;
            cbreak = 1;
            nbndtol = ceil(min(0.01*mol.nbnd,5));
            for ik = 1:nkpts
        	    idx = idx(end) + (1:nXcols(ik));
                timed = tic;
        	    [X{ik}, ev(idx), options, step_prof] = updateX(mol, H{ik}, X{ik}, Hprec{ik}, options);
                timed1 = toc(timed);
                scf_prof_history{iterscf,ik} = step_prof;
                info.t_direig = [info.t_direig,timed1];
                fprintf('Step %3d updateX time: %f\n',iterscf,timed1);
                H.eband{ik} = ev(idx);
                if options.test_proj
                    [Y0,~] = wfcX2Yspw_eband(smat_cell{ik},H{ik},X{ik});
                    ycols = ncols(Y0).'; 
                    totalband = ycols*smat_cell{ik}.dim;
                    fprintf('------------------Proj Test--------------------\n')
                    fprintf('Total bands for ik %02d = %03d\n',ik,totalband );
                    fprintf('Total calculated bands for ik %02d = %03d\n',ik,sum(ycols));
                    fprintf('Band Distribution :');  
                    disp(ycols);
                    clear Y
                end
                if options.cbreak
                    fineproj = ( totalband - mol.nbnd) < nbndtol;
                    cbreak = cbreak * fineproj;
                end      
            end

            if options.cbreak
                if cbreak
                    fprintf('Get fine band projection. Break Coarse Grid SCF.\n');
                    info.converge = true;
                    info.scfiter = iterscf;
                    info.prof_history = scf_prof_history(1:iterscf,:);
                    break;
                end
            end
 
            % Perform another Hamiltonian diagonalization for spin-down electrons in spin-unrestricted case
            if nspin == 2
                for ik = 1:nkpts
                    idx = idx(end) + (1:nXcols(ik+nkpts));
                    [X{ik+nkpts}, ev(idx),options] = updateX(mol, H{ik}, X{ik+nkpts}, Hprec{ik}, options);
                    H.eband{ik+nkpts} = ev(idx);
                end
            end
        else
            if nspin == 1||nspin == 4
        	    [X, ev, options] = updateX(mol, H, X, Hprec, options);
                H.eband = ev;
            elseif nspin == 2
                [X{1}, evup,options] = updateX(mol, H, X{1}, Hprec, options);
                [X{2}, evdw,options] = updateX(mol, H, X{2}, Hprec, options);
                ev = [evup; evdw];
                H.eband = ev;
            end
    
        end % crystal/mol for regular DFT
    else
        idx = 0;
        for ik = 1:nkpts
            idx = idx(end) + (1:nXcols(ik));
            timespw = tic;
            [YS{ik}, ev(idx),eband_spw, options,step_prof] = updateYspw(mol,H{ik},YS{ik},HprecYS{ik},options);
            ty =  toc(timespw);
            scf_prof_history{iterscf,ik} = step_prof;
            H.eband{ik} = ev(idx);
            H.eband_spw{ik} = eband_spw;
            fprintf('UpdateY Time for ik %02d: %.3f\n',ik,ty);
            info.t_spweig  = [info.t_spweig,ty];
            timey2x = tic; 
            X{ik} = wfcY2X_spw(YS{ik});
            % X{ik} = wfcY2X(YS{ik});
            info.t_y2x = info.t_y2x + toc(timey2x);
        end

    end % usespw?

 % scfiter
    [occ,mol.efermi] = getocc(ev,nocc,Tbeta,smear);
    info.occ = occ;

    if iscryst
    	idx = 0;
    	for ik = 1:nkpts*(mol.lsda+1)
    	    idx = idx(end) + (1:nXcols(ik));
            % the weight is multiplied after calculation of Entropy
            % ev(idx) = ev(idx)*wks(ik);
            X{ik}.occ = occ(idx);
    	end
    else
        if nspin == 1||nspin == 4
            X.occ = occ;
        elseif nspin == 2
            X{1}.occ = occ(1:mol.nbnd);
            X{2}.occ = occ(mol.nbnd+1:end);
        end
    end

    % Update density function rho
    rho = getcharge(mol,X,occ);

    if options.symcharge == 1
       rho = symcharge(mol,rho,S);
    end

    H.rho = rho;
    
    if strcmpi(what2mix,'rho')
        if nspin == 1
            %rhoerr = norm(rho(:)-rhoin(:))/norm(rhoin(:));
            rhoerr = 2*rhoerr_qe(mol,rho,rhoin);
        else
            rhoerr = 2*rhoerr_qe(mol,rho,rhoin);
        end
        scferr(iterscf) = rhoerr;
        fprintf('Rel Rho Err     = %20.3e\n',rhoerr);
        
        [rho,dfmat,dvmat,cdfmat] =... 
            potmixing(mol,rhoin,rho,iterscf,mixtype,...
            betamix, dfmat, dvmat, cdfmat, mixdim, brank);
    end
    
    if nspin == 1
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta*2;
    else
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta;
    end

    if iscryst
        Entropy = Entropy/nkpts;
    end
    % Kinetic energy and some additional energy terms
    if iscryst
        idx = 0;
        for ik = 1:nkpts*(mol.lsda+1)
            idx = idx(end) + (1:nXcols(ik));
            % the weight is multiplied after calculation of Entropy
            ev(idx) = ev(idx)*wks(ik);
        end
    end

    if nspin == 1
        Ekin = 2*sum(ev.*occ);
    else
        Ekin = sum(ev.*occ);
    end

    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    Ecor = getEcor(mol, rho, vtot, vion, vext);
    
    % Compute Hartree and exchange correlation energy and potential
    % using the new charge density; update the total potential
    [vhart,vxc,uxc2,rho,uxcsr]=getVhxc(mol,rho);
    
    % Update total potential
    vtot = getVtot(mol, vion, vext, vhart, vxc);
    if strcmpi(what2mix,'pot')
        if nspin == 1
            %vtoterr = norm(vtot(:)-vtotin(:))/norm(vtotin(:));
            vtoterr = rhoerr_qe(mol,vtot,vtotin);
        elseif nspin == 2||nspin == 4
            vtoterr = rhoerr_qe(mol,vtot,vtotin);
        end
        scferr(iterscf) = vtoterr;
        fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);
        
        [vtot,dfmat,dvmat,cdfmat] = ...
            potmixing(mol,vtotin,vtot,iterscf,mixtype,...
            betamix,dfmat,dvmat,cdfmat,mixdim,brank);
    end
    H.vtot = vtot;
    
    % Calculate the potential energy based on the new potential
    Ecoul = getEcoul(mol,rho,vhart);
    Exc   = getExc(mol,rho,uxc2,uxcsr);
    Etot  =  Entropy + Ekin + Eewald + Ealphat + Ecor + Ecoul + Exc;

    if ishybrid
    % Exchange energy is double counted in Ekin, need correction.
    	Exx  = getExx(X,Vexx,mol);
        if nspin == 2 || nspin == 4
            Exx = Exx/2;
        end
	Etot = Etot - Exx;
    end
    Etotvec(iterscf) = Etot;
    
    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    [cvg,resfro] = reportconverge(H,X,iterscf,maxscfiter, ...
        scferr(iterscf),scftol,verbose);
    if cvg
        info.converge = true;
        info.scfiter = iterscf;
        info.prof_history = scf_prof_history(1:iterscf,:);
        break;
    end
    deltaE = abs(Ecoul - Ecoul_old);
    Ecoul_old = Ecoul;
    if ~options.fixcgtol
        if 0
            options.cgtol = min(options.cgtol,0.01*deltaE/max(1.0,nocc));
            options.cgtol = max(options.cgtol,1e-9);
        else
            if 0
             new_cgtol = 0.01 * scferr(iterscf) / max(1.0, nocc);
             new_cgtol = min(new_cgtol, 1e-2);
             options.cgtol = max(new_cgtol, 1e-9);
            else
                options.cgtol = min(options.cgtol, 0.1*scferr(iterscf)/max(1,mol.nel));  
                options.cgtol = max(options.cgtol, 1e-10); 
            end
            
        end
    end
     fprintf('cgtol is set to %.2e\n',options.cgtol);
end

if (nspin == 4 && ~domag) || nspin == 1
    H.dv = vtot - vtotin;
else
    H.dv = rho_mix(mol,-1,vtot,vtotin); % V(out)- V(in) used to correct the forces
end

if iscryst
    X = assignoccs(X,occ);
else
    if nspin == 1||nspin == 4
        X.occ = occ;
    elseif nspin == 2
        X{1}.occ = occ(1:mol.nbnd);
        X{2}.occ = occ(mol.nbnd+1:end);
    end
end

if force
    mol.xyzforce = getFtot(mol,H,X,rho);
    info.Force = mol.xyzforce;
end

%{
if iscryst
    info.Eigvals = reshape(ev,mol.nbnd,mol.nkpts*(mol.lsda+1));
    options.ev0 = H.eband;
else
    info.Eigvals = ev;
    info.rho=rho;
    options.ev0 = ev;
end
%}
if iscryst
    info.Eigvals = cell(nkpts,1);
    startIdx = 1;
    for i = 1:nkpts
        len = nXcols(i);
        info.Eigvals{i} = sort(ev(startIdx : startIdx + len - 1));
        startIdx = startIdx + len;
    end
        options.ev0 = H.eband;
else
    info.Eigvals = ev;
    info.rho=rho;
    options.ev0 = ev;
end


if domag
    Mag = calmag(mol,rho);
end

info.Etotvec = Etotvec(1:iterscf);
info.SCFerrvec = scferr(1:iterscf);
info.Etot = Etot;
info.Eentropy = Entropy;
info.Efree = Etotvec(end) - Tbeta*Entropy;
info.Eoneelectron = Ekin+Ecor;
info.Ehart = Ecoul;
info.Exc = Exc;
info.Eewald = Eewald;
info.Efermi = mol.efermi;
if domag
    info.Mag = Mag;
end

timetot = toc(scfstart);
fprintf('Etot            = %20.13e\n', Etot);
fprintf('Entropy         = %20.13e\n', Entropy);
fprintf('Ekin            = %20.13e\n', Ekin);
fprintf('Eewald          = %20.13e\n', Eewald);
fprintf('Ealphat         = %20.13e\n', Ealphat);
fprintf('Ecor            = %20.13e\n', Ecor);
fprintf('Ehart           = %20.13e\n', Ecoul);
fprintf('Exc             = %20.13e\n', Exc);

if mol.abnd
    %ha2ev = 4.3597447222071e-18/1.602176634e-19;
    %nocc_max = find(occ>eps,1,'last');
    fprintf('Efermi         = %20.13e\n',mol.efermi);
 %   fprintf('HO energy       = %20.13e\n',ev(nocc_max));
%    fprintf('LO energy       = %20.13e\n',ev(nocc_max+1));
end

if force
    fprintf('----------------Forces for atoms----------------\n');
    force = mol.xyzforce;
    for it = 1:numel(mol.alist)
        fprintf('%5s : [%20.13e %20.13e %20.13e]\n',mol.atoms(mol.alist(it)).symbol,...
            force(it,1),force(it,2),force(it,3));
    end
end

if domag
    fprintf('----------------Magnetization----------------\n');
    if mol.lsda
        fprintf('Magtot          = %20.13e\n', Mag.totmag);
        fprintf('Magabs         = %20.13e\n', Mag.absmag);
    elseif mol.noncolin
        fprintf('Magabs          = %20.13e\n', Mag.absmag);
        fprintf('Magtot          = [%20.13e %20.13e %20.13e]\n', Mag.mx, Mag.my, Mag.mz);
    end
end

fprintf('--------------------------------------\n');
fprintf('Total time used = %20.3e\n', timetot);
fprintf('||HX-XD||_F     = %20.3e\n', resfro);

fprintf('Etot(Ry)        = %20.13e\n', Etot*2);
fprintf('Entropy(Ry)     = %20.13e\n', Entropy*2);
fprintf('Eoneelectron(Ry)= %20.13e\n', (Ekin+Ecor)*2);
fprintf('Ehart(Ry)       = %20.13e\n', Ecoul*2);
fprintf('Exc(Ry)         = %20.13e\n', Exc*2);
fprintf('Eewald(Ry)      = %20.13e\n', Eewald*2);

if iscryst % sort the H.eband and X.psi
    for ik = 1:mol.nkpts
        [H.eband{ik},sort_idx] = sort(H.eband{ik});
        X{ik}.psi = X{ik}.psi(:,sort_idx);
    end
else
    [H.eband,sort_idx] = sort(H.eband);
    X.psi = X.psi(:,sort_idx);
end

% t_transform_total = t_spw_forward + t_spw_backward;
fprintf('  │   ├─ Forward Transform (S*Y)     :  %10.4f s  [!!!]\n', t_spw_forward);
fprintf('  │   ├─ Backward Transform (S''*HX)   :  %10.4f s  [!!!]\n', t_spw_backward);
end
