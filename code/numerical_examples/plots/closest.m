function [idxs, val, othervals] = closest(pos, dim, val)
% pos = positions
% dim = search over dimnesion dim
% val = search for this value
[d, ~] = size(pos);
assert(1 <= dim && dim <= d, 'Error dimension %d out of range %d.', dim, d);
coord = pos(dim, :);
errs = abs(coord - val);
idxs = find(errs == min(errs));
val = pos(dim, idxs(1));
idxs = find(pos(dim, :) == val);
othervals = pos(3-dim, idxs);
[othervals, I] = sort(othervals);
idxs = idxs(I);
end

