% Si64 setup
% construct atoms

a2 = Atom('Si');
natom =2;
atomlist = zeros(1,natom,'Atom');
for i = 1:natom
    atomlist(i) = a2;
end

% primitive cell
coefs = [     0.000000000         0.000000000         0.000000000
     0.250000000         0.250000000         0.250000000
];
atob = 1.8897161646320724;
C = [        3.8492999077         0.0000000000         0.0000000000
        1.9246499538         3.3335915068         0.0000000000
        1.9246499538         1.1111971689         3.1429402136]*atob;

% coordinates of atoms
xyzlist = coefs*C;