import json
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# name of structure file
structure_file = 'pystructure.json'  

# read structure file
try:
    with open(structure_file, 'r') as f:
        structure_data = json.load(f)
except FileNotFoundError:
    print(f"Error: The file {structure_file} does not exist.")
    exit(1)
except json.JSONDecodeError:
    print(f"Error: The file {structure_file} is not a valid JSON file.")
    exit(1)

# create structure
structure = Structure.from_dict(structure_data)

# sym_analyze
analyzer = SpacegroupAnalyzer(structure)

# space group info
spacegroup = analyzer.get_space_group_symbol()
spacegroup_number = analyzer.get_space_group_number()
operations = analyzer.get_symmetry_operations()

# print
print(f"Spacegroup: {spacegroup}, Number: {spacegroup_number}")

# creat matlab variables
matlab_script = f"""
disp("------------------Symmetry Group-----------------")
fprintf('Space Group Symbol: %s\\n', '{spacegroup}');
fprintf('Space Group Number: %d\\n', {spacegroup_number});

% Extract symmetry operation data and save as MATLAB variables
rotation_matrices = cell({len(operations)}, 3);
translation_vectors = cell({len(operations)}, 1);
"""

# generate sym_ops
for i, op in enumerate(operations):
    rotation_matrix = np.array(op.rotation_matrix).tolist()  # 转换为列表
    translation_vector = np.array(op.translation_vector).tolist()  # 转换为列表
    
    # rotation?
    transposed_rotation_matrix = rotation_matrix

    # add to matlab variables
    matlab_script += f"""
rotation_matrices{{{i+1}, 1}} = [{', '.join(map(str, transposed_rotation_matrix[0]))}; {', '.join(map(str, transposed_rotation_matrix[1]))}; {', '.join(map(str, transposed_rotation_matrix[2]))}];  % 3x3 矩阵（转置）
rotation_matrices{{{i+1}, 2}} = trace(rotation_matrices{{{i+1}, 1}});
rotation_matrices{{{i+1}, 3}} = det(rotation_matrices{{{i+1}, 1}});
translation_vectors{{{i+1}}} = [{', '.join(map(str, translation_vector))}];
"""

# pure rotation
matlab_script += """
% Filter operations with zero translation vectors
zeroIdx = cellfun(@(x) isequal(x, [0,0,0]), translation_vectors);
indices = find(zeroIdx);
%disp('Indices of operations with zero translation vectors:');
%disp(indices);
"""

# 
with open('symmetry_operations.m', 'w') as f:
    f.write(matlab_script)

print("MATLAB script has been saved as symmetry_operations.m")
