prepare
filename = [datapath 'poisson_square_implicit_sigma_scan.h5'];
info = h5info(filename);
groups = info.Groups;
sbn = length(groups);
swn = length(groups(1).Groups);

sbs = zeros(sbn, 1);
sws = zeros(swn, 1);

pos = h5read(filename, sprintf('%s/pos', groups(1).Groups(1).Name));
N = length(pos);
x = pos(1, :);
y = pos(2, :);
anal = poisson_square_analytical(x, y);

cutoff = zeros(sbn, swn);
err = zeros(sbn, swn);
for i = 1:sbn
    for j = 1:swn
        name = groups(i).Groups(j).Name;
        fprintf('%s\n', name)
        sbs(i) = h5readatt(filename, name, 'sigmaB');
        sws(j) = h5readatt(filename, name, 'sigmaW');
        cutoff(i, j) = h5readatt(filename, name, 'cutoff');

        sol = h5read(filename, sprintf('%s/sol', name))';
        err(i, j) = max(abs(sol - anal));
    end
end

save([plotdatapath 'poisson_sigma_dependence_gau.mat'], 'err', 'cutoff', 'sbs', 'sws')