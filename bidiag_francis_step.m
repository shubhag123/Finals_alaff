function [B_next, bulge] = bidiag_francis_step(B)


[m, n] = size(B);
[k, l] = find(B(1:min(m-1,n),2:min(m,n)-1), 1);

if isempty(k)
    B_next = B;
    bulge = 0;
    return
end

% Compute the bulge
c = B(k, k) / norm(B(k:k+1, l));
s = -B(k+1, k) / norm(B(k:k+1, l));
G = [c s; -s c];

% Apply the rotation to the right of the bulge
B(k:k+1, l:n) = G' * B(k:k+1, l:n);

% Apply the rotation to the left of the bulge
B(1:k+1, k:k+1) = B(1:k+1, k:k+1) * G;

% Find the new bulge
if l < n
    B(k, l+1) = 0;
    bulge = l+1;
else
    % Otherwise, the new bulge is in the same column as the old one
    bulge = l;
end

% Compute the Givens rotation matrix that zeros out the new bulge
c = B(k, bulge) / norm(B(k:k+1, bulge));
s = -B(k+1, bulge) / norm(B(k:k+1, bulge));
H = [c -s; s c];

% Apply the Givens rotation to the left of the new bulge
if bulge < n
    B(k:k+1, bulge+1:n) = H' * B(k:k+1, bulge+1:n);
end

% Apply the Givens rotation to the right of the new bulge
B(1:k+1, bulge:bulge+1) = B(1:k+1, bulge:bulge+1) * H';

% Return the updated matrix and the new bulge position
B_next = B;
end