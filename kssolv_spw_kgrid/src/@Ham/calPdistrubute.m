n = size(eigv, 1);
maxn = 0;
for i = 1:n
    maxn = max(length(find(eigv(:, i) ~= 0)), maxn);
end
maxm = 0;
for i = 1:n
    maxm = max(length(find(eigv(i, :) ~= 0)), maxm);
end
maxm
maxn