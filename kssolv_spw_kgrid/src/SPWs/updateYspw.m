function [YS,ev,eband_spw,options,prof] = updateYspw(mol, H, YS, precYS, options)
%
% usage: [X,ev] = updateX(mol, H, X, prec, options);
%
% purupse: update the wavefunctions by computing the invariant
%          subspace associated with the lowest eigenvalues of H
%
% This is a version of GPA. GpaWF input, Gpa eigenpair output. 
% By now, only 'lobpcg' is used.
% The Transformation from Y to X is not included in this func.
global version;
%global time;
eigmethod = options.eigmethod;

% disp(eigmethod);
verbose   = options.verbose;
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
   verbose = 0;
else
   verbose = 1;
end
nYcols = ncols(YS);
nirep = size(YS.dim,1);
prof = cell(nirep,1);
eband_spw = cell(nirep,1);
switch lower(eigmethod) 
   case {'lobpcg'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      % nYcols = ncols(YS);
      ev = zeros(nYcols'*YS.dim,1);
      % ev_set = cell(Gprec.gorder,1);
      lvec_set = cell(nirep,1);
      rvec_set = cell(nirep,1);
      idx = 0;

      
     for ig = 1:nirep
        if nYcols(ig) > 0
           for j = 1:YS.dim(ig)
               if isempty(mol.efermi)
                   btype = ones(1,nYcols(ig));
               else
                   % btype = ones(1,nYcols(ig));
                   btype = (H.eband_spw{ig} < (mol.efermi + 5e-3));
               end
              if j < 2
                     idx = idx(end) + (1:nYcols(ig));
                     prec = precYS{ig};
                     if ~options.diag_profile
	                 [YS{ig}, ev(idx), lvec_set{ig}, rvec_set{ig}] = ...
                     lobpcg_cpu(H, YS{ig}, prec, cgtol, maxcgiter, verbose);
                     else
                      [YS{ig}, ev(idx), lvec_set{ig}, rvec_set{ig},prof_nirep] = ...
                     lobpcg_cpu_profile(H, YS{ig}, prec, cgtol, maxcgiter, verbose);
                      prof{ig} = prof_nirep;
                     end
                     eband_spw{ig} = ev(idx); 
              else
                     ev(idx(end) + (1:nYcols(ig))) = ev(idx);
                     idx = idx(end) + (1:nYcols(ig));
              end
           end
        end
     end
   
    case {'cg'}
      cgtol = options.cgtol;
      maxcgiter = options.maxcgiter;
      % nYcols = ncols(YS);
      ev = zeros(nYcols'*YS.dim,1);
      % ev_set = cell(Gprec.gorder,1);
      lvec_set = cell(nirep,1);
      rvec_set = cell(nirep,1);
      idx = 0;

      
     for ig = 1:nirep
        if nYcols(ig) > 0
           for j = 1:YS.dim(ig)
               if isempty(mol.efermi)
                   btype = ones(1,nYcols(ig));
               else
                   % btype = ones(1,nYcols(ig));
                   btype = (H.eband_spw{ig} < (mol.efermi + 5e-3));
               end
              if j < 2
                     idx = idx(end) + (1:nYcols(ig));
                     prec = precYS{ig};
                     if ~options.diag_profile
	                 [YS{ig}, ev(idx), lvec_set{ig}, rvec_set{ig}] = ...
                     cgdiagg_qe_matlab(H,YS{ig}, prec, cgtol, maxcgiter, btype);
                     else
                      [YS{ig}, ev(idx), lvec_set{ig}, rvec_set{ig},prof_nirep] = ...
                     cgdiagg_qe_matlab_profile(H,YS{ig}, prec, cgtol, maxcgiter, btype); 
                        prof{ig} = prof_nirep; 
                     end
                     eband_spw{ig} = ev(idx); 
              else
                     ev(idx(end) + (1:nYcols(ig))) = ev(idx);
                     idx = idx(end) + (1:nYcols(ig));
              end
           end
        end
     end

    
   case {'davidson_qe'}
       cgtol = options.cgtol;
       maxcgiter = options.maxcgiter;
       
       gorder = size(YS.gmat,1);
       if ~options.spwdavidprec || ~isfield(options,'h_diag')
           options.h_diag = cell(mol.nkpts,gorder);
           options.s_diag = cell(mol.nkpts,gorder);
       end
       
       ik = YS.ik;
       if isempty(options.h_diag{ik,1}) 

           npol = 1 + mol.noncolin;
           % vion_g = fft3(H.vion)/mol.n1/mol.n2/mol.n3;
           % vion_g1 = vion_g(1);
           vion_g1 = mean(H.vion(:));
           h_diag0 = repmat(H.gkin,npol,1) + vion_g1;
           [h_diag0, s_diag0] = addnlc(mol,h_diag0,H.vnlsign,H.vnlmat);
           
           for ig = 1:gorder
               Pmat = YS.gmat{ig,1};
               options.h_diag{ik,ig} = sum(abs(Pmat).^2 .* h_diag0).';
               options.s_diag{ik,ig} = sum(abs(Pmat).^2 .* s_diag0).';
           end
           options.spwdavidprec = 1;
       end
       
     ev = zeros(nYcols'*YS.dim,1); 
     idx = 0;
     for ig = 1:nirep
        if nYcols(ig) > 0
           for j = 1:YS.dim(ig)
              if j < 2
                     idx = idx(end) + (1:nYcols(ig));
                     if ~options.diag_profile
	                 [YS{ig}, ev(idx)] = ...
            davidson_qe_spw(mol, H,  YS{ig}, cgtol, maxcgiter,options.iterscf,options.h_diag{ik,ig},options.s_diag{ik,ig});
                     else
                         [YS{ig}, ev(idx),prof_nirep] = ...
            davidson_qe_spw_profile(mol, H,  YS{ig}, cgtol, maxcgiter,options.iterscf,options.h_diag{ik,ig},options.s_diag{ik,ig});
                     prof{ig} = prof_nirep;
                     end
                     eband_spw{ig} = ev(idx); 
              else
                     ev(idx(end) + (1:nYcols(ig))) = ev(idx);
                     idx = idx(end) + (1:nYcols(ig));
              end
           end
        end
     end
    
    
    
    
    case {'eigs'}
      eigstol = options.eigstol;
      maxeigsiter = options.maxeigsiter;
      idx = 0;
      ev = zeros(sumel(nYcols),1);
    for ig = 1:nirep
      if YS.gconj{ig} > 0
             idx = idx(end) + (1:nYcols(ig));
             %prec = precYS{ig};
      [YS{ig}, ev(idx)] = diagbyeigs(mol, H, nYcols(ig), eigstol, maxeigsiter);
      else
             ev(idx(end) + (1:nYcols(YS.gconjref{ig}))) = ev(idx);
             idx = idx(end) + (1:nYcols(YS.gconjref{ig}));
      end 
    end
   otherwise
      disp('Unknown method for diagonalizing H! Use eigs');
      [X, ev] = diagbyeigs(mol, H, ncol, eigstol, maxeigsiter);
end
%
if 0%( verbose ==1 )
   HX = H*X;
   G = X'*HX;
   R = HX-X*G;
   ev = sort(real(eig(G)));
   for j = 1:ncol
      resnrm(j) = norm(R(:,j));
      %fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', ...
      %        j, ev(j), resnrm(j));
   end
end
