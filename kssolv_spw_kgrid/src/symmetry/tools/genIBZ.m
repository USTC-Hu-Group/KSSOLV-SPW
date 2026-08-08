function cry = genIBZ(cry,S)
% genIBZ calculates the Irreducible k-points of the selected k-points
% according to the symmetry of the system.
% By now, this code is valid for nspin = 1 and 2.
% When nspin = 1, we use the full rotation symmetry.
% When nspin = 2, we consider the rotation and time-reversal symmetry (TRS)
% of the initial magnetic momentum. Only those rotation without TRS will be
% used here.
% Notice that, the choice of k-points is different but equivlent from IBZ.

% Liangyu Wang, Yuanfan Xiong 2026/01/12
   if ~isa(cry,'Crystal') && isa(S,'Sym')
       error('IBZ Input should be Crystal & Sym class.\n');
   end
   rec = 2*pi*inv(cry.supercell)';
   kpts = cry.kpts/rec;
   rot_rec = rotReal2Rec(cry,S);
   n_sym = size(rot_rec,1);
   n_kpts = cry.nkpts;
   tol = 1e-5;


   visited = false(n_kpts, 1);
   ibz_kpts = [];
   ibz_counts = [];
   if sum(abs(cry.initmag.amag)) > 0
       use_trs = 0;
   else
       use_trs = 1;
   end
   fprintf('Start generating IBZ. nkpts=%3d, nsym=%2d, TRS=%2d\n', ...
       n_kpts,n_sym,use_trs);
   ik_count = 0;
   for i = 1:n_kpts
      if visited(i)
          continue;
      end
      current_k = kpts(i,:);
      ibz_kpts = [ibz_kpts;current_k];
      star_indices = [];
      star_indices_trs = [];
      for s = 1:n_sym
          RS = rot_rec{s};
          k_rot = current_k*RS.';
          idx = find_equivalent_k(k_rot,kpts,tol);
          star_indices = [star_indices;idx];

          if use_trs
              idx_trs = find_equivalent_k(k_rot,-kpts,tol);
              star_indices = [star_indices;idx_trs];
          end
      end
      unique_indices = unique(star_indices);
      count = length(unique_indices);
      ibz_counts = [ibz_counts;count];
      visited(unique_indices) = true;
      ik_count = ik_count + 1;
      ibz_ref = build_kref(cry,star_indices,use_trs);
      cry.ibz_ref{ik_count,1} = ibz_ref; 
   end
   if sum(ibz_counts) == n_kpts
       fprintf('IBZ check passed! nkpts_ibz=%3d \n',length(ibz_counts));
   else
       error('IBZ counts number != nkpts.\n');
   end
   ibz_weights = ibz_counts / n_kpts;
    
   % if cry.nspin == 1
   %     ibz_weights = ibz_weights;
   % end
   cry.kpts = ibz_kpts * rec;
   cry.nkpts = size(ibz_kpts,1);
   cry.wks = ibz_weights;
   if cry.nspin == 2
       cry.wks = repmat(cry.wks,2,1);
   end
end 

function indices = find_equivalent_k(target_k, all_kpts, tol)
    delta = all_kpts - target_k;
    delta = delta - round(delta);
    dist_sq = sum(delta.^2, 2);
    indices = find(dist_sq < tol);
end

function data_store = build_kref(cry,star_index,usetrs)
kpts = cry.kpts/(2*pi*inv(cry.supercell)');
if usetrs
    star_index = reshape(star_index,2,[])';
    star_index_trs = star_index(:,2);
    star_index = star_index(:,1);
end
unique_idx = unique(star_index);
nkpt= length(unique_idx);
data_store = cell(nkpt,3);
for i = 1:nkpt
    data_store{i,1} = kpts(unique_idx(i),:);
    data_store{i,2} = find(star_index == unique_idx(i));
end
if usetrs
    unique_idx = unique(star_index_trs);
   for i = 1:nkpt
       data_store{i,3} = find(star_index_trs == unique_idx(i));
   end
end

end