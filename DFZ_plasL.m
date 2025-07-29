function l_plas = DFZ_plasL(tauP_in, K_local, l_DFZ, l_max)
% calculate the length of the plastic zone based on the DFZ condition

arguments (Input)
    tauP_in % Lattice friction stress, used in residue calculation
    K_local % Stress intensity factor, material property
    l_DFZ % Initial length of the plastic zone
    l_max % Maximum allowable length of the plastic zone
end

arguments (Output)
    l_plas
end

max_iter = 10000; % maximum number of iterations
iter_count = 0; % iteration counter

tol = 1e-6; % convergence tolerance
lp_min = l_DFZ; % start from l_DFZ
lp_max = l_max; % upper bound for l_p

L_max = lp_max + l_DFZ;
k_max = sqrt(lp_max / L_max);
m_max = k_max^2;
[FK_max, ~] = ellipke(m_max);
Residue_max = tauP_in * pi - K_local /  sqrt(2 * pi * L_max) * FK_max;

% If Residue_max < 0, expand the range of lp_max
while Residue_max < 0
    lp_max = lp_max * 2; % double the upper bound
    L_max = lp_max + l_DFZ;
    k_max = sqrt(lp_max / L_max);
    m_max = k_max^2;
    [FK_max, ~] = ellipke(m_max);
    Residue_max = tauP_in * pi - K_local / sqrt(2 * pi * L_max) * FK_max;
    iter_count = iter_count + 1;
    if iter_count > max_iter
        disp('Maximum iterations reached without convergence.');
        break
    end
end

while (lp_max - lp_min) / lp_max > tol
    lp_guess = (lp_min + lp_max) / 2; % midpoint
    L_guess = lp_guess + l_DFZ; % update L
    k = sqrt(lp_guess / L_guess);
    m = k^2;
    [FK, ~] = ellipke(m); % calculate FK
    Residue = tauP_in * pi - K_local / sqrt(2 * pi * L_guess) * FK;

    if Residue > 0
        lp_max = lp_guess; % narrow the range to the lower half
    else
        lp_min = lp_guess; % narrow the range to the upper half
    end
    iter_count = iter_count + 1;
    if iter_count > max_iter
        fprintf('Maximum number of iterations reached.\n');
        break;
    end
end

l_plas = lp_guess;

end