function writeData(filename, data, G_idx, rho, mu_idx, rc, M_idx, sigma, maxtilt)
    cell1 = {"v", 0.25; "Patm", 100000; "g", 9.81; "r", 500; "G", "10^"+G_idx; "rho", rho; "mu", "10^"+mu_idx; "rc", rc};
    writecell(cell1, filename, "WriteMode", "append", "Delimiter", ",");
    writecell({"",""}, filename, "WriteMode", "append", "Delimiter", " ");
    cell2 = {"M", "10^"+M_idx; "sigma", sigma};
    writecell(cell2, filename, "WriteMode", "append", "Delimiter", ",");
    writecell({"",""}, filename, "WriteMode", "append", "Delimiter", " ");
    cell3 = {"tilt_erupt", maxtilt+"nrad"};
    writecell(cell3, filename, "WriteMode", "append", "Delimiter", ","); 
    writecell({"",""}, filename, "WriteMode", "append", "Delimiter", " ");
    writematrix(data, filename, "WriteMode", "append");
end
