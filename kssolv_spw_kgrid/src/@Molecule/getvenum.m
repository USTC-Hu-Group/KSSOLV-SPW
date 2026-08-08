function venums = getvenum(mol)
atoms = mol.atoms;
ntypes = length(atoms);
venums = zeros(1,ntypes);
for it = 1:ntypes
    pp = PpData(atoms(it),mol.lspinorb);
    venums(it) = pp.venum;
end

end