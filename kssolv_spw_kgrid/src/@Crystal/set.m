%{
function cry = set(cry,varargin)
% CRYSTAL/SET Set function for crystal class
%    cry = SET(cry,str1,field1,str2,field2,...) returns a crystal class of
%    the given fields with respect to the name strings.
%
%    See also Molecule, Crystal.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Crystal.set');
end

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,cry.alist] = unique([value.anum]);
        cry.atoms = value(IA);
        for ia = 1:numel(value)
            cry.initmag.amag = [value.amag];
            cry.initmag.theta = [value.theta];
            cry.initmag.phi = [value.phi];
        end
        continue;
    end
    
    if strcmpi(attr_name,'autokpts')
        nkx = value(1);
        nky = value(2);
        nkz = value(3);
        if numel(value) ==3
            skx = 0;
            sky = 0;
            skz = 0;
        else
            skx = value(4);
            sky = value(5);
            skz = value(6);
        end
		[I,J,K] = ndgrid((0:nkx-1)-((0:nkx-1) >= nkx/2)*nkx, ...
        	(0:nky-1)-((0:nky-1) >= nky/2)*nky, ...
        	(0:nkz-1)-((0:nkz-1) >= nkz/2)*nkz);
        pregkx = (I(:)+skx)/nkx;
        pregky = (J(:)+sky)/nky;
        pregkz = (K(:)+skz)/nkz;
        cry.nkxyz = [nkx,nky,nkz];
        cry.kpts = [pregkx pregky pregkz];
        continue;
    end
    
    cry.(attr_name) = value;
end

end
%}
function cry = set(cry,varargin)
% CRYSTAL/SET Set function for crystal class
%    cry = SET(cry,str1,field1,str2,field2,...) returns a crystal class of
%    the given fields with respect to the name strings.
%
%    See also Molecule, Crystal.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Crystal.set');
end

% [新增]：记录是否修改了网格参数
need_finalize = false;

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,cry.alist] = unique([value.anum]);
        cry.atoms = value(IA);
        for ia = 1:numel(value)
            cry.initmag.amag = [value.amag];
            cry.initmag.theta = [value.theta];
            cry.initmag.phi = [value.phi];
        end
        continue;
    end
    
    if strcmpi(attr_name,'autokpts')
        nkx = value(1);
        nky = value(2);
        nkz = value(3);
        if numel(value) ==3
            skx = 0;
            sky = 0;
            skz = 0;
        else
            skx = value(4);
            sky = value(5);
            skz = value(6);
        end
		[I,J,K] = ndgrid((0:nkx-1)-((0:nkx-1) >= nkx/2)*nkx, ...
        	(0:nky-1)-((0:nky-1) >= nky/2)*nky, ...
        	(0:nkz-1)-((0:nkz-1) >= nkz/2)*nkz);
        pregkx = (I(:)+skx)/nkx;
        pregky = (J(:)+sky)/nky;
        pregkz = (K(:)+skz)/nkz;
        cry.nkxyz = [nkx,nky,nkz];
        cry.kpts = [pregkx pregky pregkz];
        continue;
    end
    
    % 常规属性赋值
    cry.(attr_name) = value;
    
    % [新增逻辑]：监控 ecut 或 supercell 是否被修改，并清空依赖项
    if strcmpi(attr_name, 'ecut') || strcmpi(attr_name, 'supercell')
        cry.ecut2 = [];
        cry.n1 = [];
        cry.n2 = [];
        cry.n3 = [];
        cry.vext = [];
        need_finalize = true;
    end
end

% [新增逻辑]：如果核心参数被修改，重新调用 finalize
if need_finalize && ~isempty(cry.ecut) && ~isempty(cry.supercell)
    % 【关键修复】：备份 k 点信息，防止 finalize 内部的坐标变换污染数据
    temp_kpts  = cry.kpts;
    temp_nkpts = cry.nkpts;
    temp_wks   = cry.wks;
    
    % 将 kpts 置为 [0 0 0] 骗过 finalize，这样乘以矩阵后依然是 0
    cry.kpts = [0 0 0]; 
    
    % 执行底层网格的重新推导
    cry = finalize(cry);
    
    % 【关键修复】：原样恢复 k 点信息
    cry.kpts  = temp_kpts;
    cry.nkpts = temp_nkpts;
    cry.wks   = temp_wks;
end

end