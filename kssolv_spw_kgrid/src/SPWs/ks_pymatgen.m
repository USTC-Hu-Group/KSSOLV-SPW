function [rotation_matrices,spacegroup_number ] = ks_pymatgen(cry)
% ks_pymatgen uses Python package Pymatgen to generate the symmetry
% operators (sym_ops) of input structure in KSSOLV format (Crystal Class). 
% It displays the Space Group ID, and outputs the sym_ops in matrix form.
%-------------------------------------------------------------------------%
% Mind: python & its package pymatgen are needed.
%-------------------------------------------------------------------------%
% 'rotation_matrices' contains the matrix form of Rotation part of sym_ops
% (3*3 matrix in the basis of lattice vectors),  the factor translation part of sym_ops.
% (1*3 vector in the basis of lattice vectors), trace, determinant. (cell)
%
% 'space_group_number' means as it reads.
%
% Edit by Liangyu Wang (20251009)

% 0. generate atomlist, C, coefs form cry
natom = sum(cry.natoms);
atomlist = zeros(1,natom,'Atom');
for i = 1:natom
    atomlist(i) = cry.atoms(cry.alist(i),1);
end
C = cry.supercell;
coefs = cry.xyzlist/C;

% 1. define atomlist(py)
species = arrayfun(@(atom) atom.symbol, atomlist, 'UniformOutput', false); 
species = py.list(species);

% 2. convert lattice parameters (C) to lattice_matrix(py)
lattice_matrix = py.list({ ...
    py.list({C(1,1), C(1,2), C(1,3)}), ...
    py.list({C(2,1), C(2,2), C(2,3)}), ...
    py.list({C(3,1), C(3,2), C(3,3)}) ...
});

% 3. convert atom coordinates (coefs) to coords(py)
coords = py.list();  
for i = 1:size(coefs, 1)
    coord = py.list({coefs(i, 1), coefs(i, 2), coefs(i, 3)});
    coords.append(coord);
end

% 4. convert to pymatgen format & add tolerance
pymatgen = py.importlib.import_module('pymatgen.core');
Lattice = pymatgen.Lattice;
Structure = pymatgen.Structure;

% create lattice
lattice = Lattice(lattice_matrix);

% create structure
structure = Structure(lattice, species, coords);

% import SpacegroupAnalyzer
sga_module = py.importlib.import_module('pymatgen.symmetry.analyzer');
SpacegroupAnalyzer = sga_module.SpacegroupAnalyzer;

% set tolerance (symprec)
% tolerance = 1e-3; 

%  SpacegroupAnalyzer 
sga = SpacegroupAnalyzer(structure);
%sga = SpacegroupAnalyzer(structure, tolerance);

% generate space group information
spacegroup = sga.get_space_group_symbol();
spacegroup_number = int32(sga.get_space_group_number());

% generate symmetry operators
symm_ops = sga.get_symmetry_operations();

% Output space group info
disp("------------------Symmetry Group-----------------")
fprintf('Space Group Notation: %s\n', char(spacegroup));
fprintf('Space Group ID: %d\n', spacegroup_number);

% save symmetry operators & save as matrix in MATLAB format 
rotation_matrices = cell(length(symm_ops), 3);
translation_vectors = cell(length(symm_ops), 1);

    for i = 1:length(symm_ops)
        op = symm_ops{i};
        
        % generate rotation & translation matrix
        rotation_matrix = op.rotation_matrix;
        translation_vector = op.translation_vector;
        
        % convert to MATLAB format 
        rot_mat = double(py.array.array('d', py.numpy.nditer(rotation_matrix)));
        rotation_matrices{i,1} = reshape(rot_mat,3,3)';
        rotation_matrices{i,3} = trace(rotation_matrices{i,1});
        rotation_matrices{i,4} = det(rotation_matrices{i,1});
       rotation_matrices{i,2} = double(py.array.array('d', py.numpy.nditer(translation_vector)));
        % rotation_matrices{i,6} = translation_vectors{i};
    end
%reshaped_rotation_mat = cellfun(@(x) reshape(x, 3, 3),  rotation_matrices, 'UniformOutput', false);
% zeroIdx = cellfun(@(x) isequal(x, [0,0,0]), translation_vectors);
% indices = find(zeroIdx);
end