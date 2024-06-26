function output = read_Scaffold_in(filename)

    fid = importdata(filename);
    Scaffold.Scaffold_name = fid.textdata(1);
    Scaffold.inner_radius = fid.data(1);
    Scaffold.thickness = fid.data(21);
    Scaffold.cp1_1 = fid.data(2);
    Scaffold.cp1_2 = fid.data(12);
    Scaffold.cp2_1 = fid.data(22);
    Scaffold.cp2_2 = fid.data(13);
    Scaffold.eta1 = fid.data(4);
    Scaffold.eta2 = fid.data(14);
    Scaffold.g1 = fid.data(5);
    Scaffold.g2 = fid.data(25);
    Scaffold.rho_hat1 = fid.data(6);
    Scaffold.rho_hat2 = fid.data(16);
    Scaffold.vfrac1 = fid.data(7);
    Scaffold.vfrac2 = fid.data(17);
    Scaffold.k1 = fid.data(8);
    Scaffold.zeta1 = fid.data(18);
    Scaffold.k2 = fid.data(28);
    Scaffold.zeta2 = fid.data(9);
    Scaffold.fd1 = fid.data(10);
    Scaffold.fd2 = fid.data(20);
   
    output = Scaffold;
end