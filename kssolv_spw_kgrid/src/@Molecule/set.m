%{
function mol = set(mol,varargin)
% MOLECULE/SET Set function for molecule class
%    mol = SET(mol,str1,field1,str2,field2,...) returns a molecule class of
%    the given fields with respect to the name strings.
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Molecule.set');
end

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,mol.alist] = unique([value.anum]);
        mol.atoms = value(IA);
        % Store atomic magnetic moments in mol
        mol.initmag.amag = [value.amag];
        mol.initmag.theta = [value.theta];
        mol.initmag.phi = [value.phi];
        continue;
    end
    
    mol.(attr_name) = value;
end

end
%}
function mol = set(mol,varargin)
% MOLECULE/SET Set function for molecule class
%    mol = SET(mol,str1,field1,str2,field2,...) returns a molecule class of
%    the given fields with respect to the name strings.
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Molecule.set');
end

% 新增标志位：记录是否修改了会影响网格尺寸的核心参数
need_finalize = false; 

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,mol.alist] = unique([value.anum]);
        mol.atoms = value(IA);
        % Store atomic magnetic moments in mol
        mol.initmag.amag = [value.amag];
        mol.initmag.theta = [value.theta];
        mol.initmag.phi = [value.phi];
        continue;
    end
    
    % 执行标准赋值
    mol.(attr_name) = value;
    
    % [新增逻辑]：如果修改了 ecut，则清空相关的衍生变量以触发重新计算
    if strcmpi(attr_name, 'ecut')
        mol.ecut2 = [];
        mol.n1 = [];
        mol.n2 = [];
        mol.n3 = [];
        mol.vext = []; % vext 的矩阵大小依赖于 n1,n2,n3，也必须清空重置
        need_finalize = true;
    end
    
    % [新增逻辑]：如果修改了 supercell，同样会影响空间离散点数量 n1,n2,n3
    if strcmpi(attr_name, 'supercell')
        mol.n1 = [];
        mol.n2 = [];
        mol.n3 = [];
        mol.vext = [];
        need_finalize = true;
    end
end

% [新增逻辑]：如果修改了联动参数，且必要属性齐备，则重新调用 finalize 自动推导缺省值
if need_finalize && ~isempty(mol.ecut) && ~isempty(mol.supercell)
    mol = finalize(mol);
end

end