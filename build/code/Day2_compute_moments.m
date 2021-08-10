function m = Day2_compute_moments(prof, inv)
    m1 = mean(prof, 'all');
    m2 = mean(inv, 'all');
    m = [m1 m2];
end