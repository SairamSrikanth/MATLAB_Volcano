function [G_idx, rho, mu_idx, rc, M_idx] = randomVolcano(n)
    rng("shuffle");
    G_idx=round(9+rand(1)*3, 2);
    rho=round(rand(1)*1000+2000,-2);
    mu_idx=round(rand(1)*3+2, 2);
    rc=round(rand(1)*20+10);
    for i = 1:n
        M_idx(i)=round(rand(1)*2+4.7, 2);
    end
end
