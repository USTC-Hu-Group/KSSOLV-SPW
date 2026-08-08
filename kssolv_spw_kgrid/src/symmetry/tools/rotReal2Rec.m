function rot_rec = rotReal2Rec(cry,rotmat)
if ~isa(cry,'Crystal')
   error('Input cry should be Crystal class.');  
end

real2rec = cry.supercell*cry.supercell';
real2recinv = inv(real2rec);
if isnumeric(rotmat)
    rot_lat = rotmat;
    rot_rec = round(real2rec*rot_lat*real2recinv);
elseif iscell(rotmat)
    rot_rec = cellfun(@(x) round(real2rec*x*real2recinv),rotmat(:,1), ...
        'UniformOutput',false);
elseif isa(rotmat,'Sym')
    rot_rec = cellfun(@(x) round(real2rec*x*real2recinv), ...
        rotmat.symops(:,1), ...
        'UniformOutput',false);
end

end