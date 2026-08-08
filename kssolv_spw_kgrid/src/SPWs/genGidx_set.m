function rkidx_set = genGidx_set(cry,Xgc,ik)
    
    ngen = size(Xgc,1);
    rkidx_set = cell(ngen,1);
    
    % Transform matrix representation of sym_ops to Reciprocal lattice vector basis.
    Rec = transpose(2*pi*inv(cry.supercell));
    kpt = cry.kpts(ik,:)/Rec;

    real2rec = cry.supercell*cry.supercell';
    Xgc_r = cellfun(@(x) round(real2rec*x(1:3,1:3)*inv(real2rec)),Xgc(:,1),'UniformOutput', false);
    
    % generate plane-wave basis and its coordinates under Rec lattice vector.
    grid = Ggrid(ik,cry);
    kcoef = cat(2,grid.gkx,grid.gky,grid.gkz);  
    klabel =  round(kcoef / Rec);
    fprintf("Start generating <Rotation matrix> representation of %d generators.\n",ngen);
    for i = 1:ngen
        rotmat = Xgc_r{i,1}';
        G0 = round(kpt*rotmat - kpt);
        klabelc = klabel*rotmat + G0;
        [found, rkidx] = ismember(klabelc,klabel,'rows');
        rkidx_set{i} = int32(rkidx);
        if ~all(found)
            error(' %d can not find conter part at No.%d generator', sum(~found),i);
        end
    end
end