function [X0, H] = genX0_random(mol,H)
% The random generation method of wavefunctions
%  Add k-grid by Liangyu Wang, USTC, 260323
ishybrid = H.ishybrid;
H.ishybrid = 0;
if isa( mol, 'Crystal' )
    nXcols = mol.nbnd*ones(mol.nkpts,1);
    nkpts = mol.nkpts;
else
    nXcols = mol.nbnd;
    nkpts = 1;
end
n1   = mol.n1;
n2   = mol.n2;
n3   = mol.n3;
nspin = mol.nspin;

tpiba = 2*pi/mol.alat;
info.first = true;
npol = 1+mol.noncolin;

if nspin == 2
    Qcell = cell(nkpts*2,1);
else
    Qcell = cell(nkpts,1);
end


nglist = zeros(mol.nkpts,1);
idxnzcell = cell(mol.nkpts,1);

for ik = 1:nkpts  

  grid = Ggrid(ik,mol);
  idxnz = grid.idxnz;
  idxnzcell{ik} = idxnz;

  [gkk,id] = sort(grid.gkk);
  % map = zeros(numel(id),1);
  % for i = 1:numel(id)
  %     for j = 1:numel(id)
  %         if id(j) == i
  %             map(i) = j;
  %         end
  %     end
  % end
  [~, map] = sort(id);
  
  npw = size(idxnz,1);
  nglist(ik) = npw;

  nbnd = nXcols(ik);
  psif = zeros(npw,npol,nbnd);  
  for j = 1:nbnd
      if isa(mol,'Crystal')
          gg = sumel(mol.kpts(ik,:).^2);
      else
          gg = 0;
      end
      for ipol = 1:npol
          for ig = 1:npw
              [rr,info] = randy(info);
              [arg,info] = randy(info);
              arg = 2*pi*arg;
              psif(ig,ipol,j) = complex(rr*cos(arg),rr*sin(arg))/(1+gkk(ig)/...
                  tpiba^2+gg/tpiba^2);
          end
      end
  end

  if ~mol.noncolin
      Qcell{ik} = squeeze(psif(map,1,:));
  else
      psif = psif(map,:,:);
      Qcell{ik} = zeros(npw*npol,nbnd);
      Qcell{ik}(1:npw,1:nbnd) = psif(:,1,:);
      Qcell{ik}(npw+1:end,1:nbnd) = psif(:,2,:);
  end

  if nspin == 2
     psif_spin2 = zeros(npw, nbnd);
     for j = 1:nXcols(ik)  
         for ig = 1:npw
             % 保持原代码逻辑：两次 randy(info) 获取 rr
             [~, info] = randy(info); 
             [rr, info] = randy(info);
             [arg, info] = randy(info);
             arg = 2*pi * arg;
             psif_spin2(ig, j) = complex(rr*cos(arg), rr*sin(arg)) / ...
                 (1 + gkk(ig)/tpiba^2 + gg/tpiba^2);
         end
     end
     Qcell{ik + nkpts} = psif_spin2(map, :);
  end
end

% Rotate wavefunctions
if nspin == 2
    if isa( mol, 'Crystal' )
        X0 = BlochWavefun(Qcell,n1,n2,n3,idxnzcell,mol.wks,mol.nspin);
        H.eband = cell(nkpts*2,1);
        for ik = 1:nkpts
            fprintf('step1 %d \n',ik)
            [X0{ik},H.eband{ik}] = rotate_wfc(mol, H{ik}, X0{ik});
            fprintf('step2 %d \n',ik)
            [X0{ik+nkpts}, H.eband{ik+nkpts}] = rotate_wfc(mol, H{ik}, X0{ik+nkpts});
        end
    else
        X0 = cell(2,1);
        X0{1} = Wavefun(Qcell{1},n1,n2,n3,idxnzcell{1},1);
        X0{2} = Wavefun(Qcell{2},n1,n2,n3,idxnzcell{1},2);
        H.eband = zeros(nXcols*2,1);
        [X0{1}, H.eband(1:nXcols)] = rotate_wfc(mol, H, X0{1});
        % yinxl 2025.12.31
        [X0{2}, H.eband(nXcols+1:end)] = rotate_wfc(mol, H, X0{2});
        % X0{2} = X0{1};
        % H.eband(nXcols+1:end) = H.eband(1:nXcols);
        % end yinxl 2025.12.31
    end
else
    if isa( mol, 'Crystal' )
        X0 = BlochWavefun(Qcell,n1,n2,n3,idxnzcell,mol.wks,mol.nspin);
        H.eband = cell(nkpts,1);
        for ik = 1:nkpts
            [X0{ik}, H.eband{ik}] = rotate_wfc(mol, H{ik}, X0{ik});
        end   
    else   
        X0 = Wavefun(Qcell{1},n1,n2,n3,idxnzcell{1});
        [X0, H.eband] = rotate_wfc(mol, H, X0); 
    end
end
H.ishybrid = ishybrid;
end











