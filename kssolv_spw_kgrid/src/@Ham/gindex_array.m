function glabel0 = gindex_array(glabel,cry)
n1 = cry.n1;
n2 = cry.n2;
n3 = cry.n3;
nn = [n1,n2,n3];
glabel0 = zeros(size(glabel));
for i = 1:size(glabel,1)
    for j = 1:size(glabel,2)
        glabel0(i,j) = gindex(glabel(i,j),nn(j));
    end
end
