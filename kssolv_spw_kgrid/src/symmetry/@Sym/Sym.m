classdef Sym
    % Sym KSSOLV class for symmetry.
    properties (SetAccess = public)
       sgidx
       symops
       symops_trs
       time_reversal
       issupercell
       ismagnetic
       bilbao
    end
    
    methods
        function S = Sym(cry)
            if ~isa(cry,'Crystal')
                error('Input should be "Crystal" class for "Sym" class.\n');
            else
                [rotmat,sgidx] = ks_pymatgen(cry);
                S.sgidx = sgidx;
                S.symops = rotmat(:,1:2);
            end
            sym_e_idx = cell2mat(rotmat(:,3)) == 3;
            trans_mat = cell2mat(rotmat(sym_e_idx,2));
            issupercell = sum(max(abs(trans_mat))) > 0;
            if issupercell
                S.issupercell = 1; 
                if 0
                    rot_data = cell2mat(cellfun(@(x) x(:), rotmat(:,1), ...
                        'UniformOutput', false));
                    rot_data = reshape(rot_data, 9, []);     
                    [~, rot_ref_idx] = rref(rot_data);
                    num_mats = size(rot_data, 2);
                    S.rotref = false(num_mats,1);
                    S.rotref(rot_ref_idx) = true;  
                else
                    trans = sum(abs(cell2mat(rotmat(:,2))),2);
                    idx_ztrans = trans < 1e-5;
                    S.symops = rotmat(idx_ztrans,1:2);
                end
            else
                S.issupercell = 0;
            end
            

            %-------mag-------%
            if cry.nspin == 2
                S.ismagnetic = 1;
            elseif cry.nspin == 4
                S.ismagnetic = 4;
            else
                S.ismagnetic = 0;
            end
            
            if sum(abs(cry.initmag.amag)) < 1e-5
                S.time_reversal = 1;
            else
                S.time_reversal = 0;
            end
            nsym = size(S.symops,1);
            
            xyzdlist = mod(cry.xyzlist/cry.supercell,1);
            xyzdlist = round(xyzdlist,8);
            xyzdlist = mod(xyzdlist,1);
            amag0 = cry.initmag.amag;
            allowed_sym = false(nsym, 1);
            allowed_sym_trs = false(nsym,1);
            
            for i = 1:nsym
                rotmat = S.symops{i,1};
                tau = S.symops{i,2};
                rxyzdlist = mod(xyzdlist * rotmat.' + tau,1); 
                rxyzdlist = round(rxyzdlist,8);
                rxyzdlist = mod(rxyzdlist,1);
                
                delta = mod(reshape(rxyzdlist, [], 1, 3) - reshape(xyzdlist, 1, [], 3) + 0.5, 1) - 0.5;
                
                [min_dist, idxr] = min(max(abs(delta), [], 3), [], 2);
                
                idxr(min_dist > 1e-4) = 0;
            
                if any(idxr == 0)
                    continue; % 
                end
                
                % Equal to 
                %[~,idxr] = ismembertol(rxyzdlist,xyzdlist,1e-4,'ByRows',true);

                % 磁性对称性判断 (常规对称性 和 附加时间反演对称性)
                if sum(abs(amag0(idxr) - amag0)) < 1e-5
                    allowed_sym(i) = true;
                elseif sum(amag0(idxr) + amag0) < 1e-5
                    allowed_sym_trs(i) = true;
                end                
            end

           % 根据筛选结果更新对称操作
           S.symops_trs = S.symops(allowed_sym_trs,:);
           S.symops = S.symops(allowed_sym,:);
        end

        function S = Sbilbaoinfo(S,idx,varargin)
           R = [];
           t = [];

           n_extra = length(varargin);

           if n_extra == 0
               R = eye(3);
               t = zeros(1,3);
           
           elseif n_extra == 1
               input1 = varargin{1};
               sz = size(input1);

               if isequal(sz,[4,4])
                   R = input1(1:3,1:3);
                   t = input1(1:3,4).';
               elseif isequal(sz,[3,3])
                   R = input1;
                   t = zeros(1,3);
               else
                   error('input error!')
               end
           elseif n_extra == 2
               R = varargin{1};
               t = varargin{2};
               t = reshape(t,1,3);
           end

           info = cell(3,1);
           info{1} = idx;
           info{2} = R;
           info{3} = t;

           S.bilbao = info;              
        end
    end
end