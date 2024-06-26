function write_Scaffold_in(filename)

    fid = fopen('Scaffold_in_', 'w');
    a = filename.Scaffold_name;
    b = filename.inner_radius;
    c = filename.thickness;
    d = filename.cp1_1;
    e = filename.cp1_2;
    f = filename.cp2_1;
    g = filename.cp2_2;
    h = filename.eta1;
    i = filename.eta2;
    j = filename.g1;
    k = filename.g2;
    l = filename.rho_hat1;
    m = filename.rho_hat2;
    n = filename.vfrac1;
    o = filename.vfrac2;
    p = filename.k1;
    q = filename.zeta1;
    r = filename.k2;
    s = filename.zeta2;
    t = filename.fd1;
    u = filename.fd2;
   
    fprintf(fid, '%s \n', a{:});
    fprintf(fid, '%f', b);
    fprintf(fid, '		%f \n', c);
    fprintf(fid, '%f', d);
    fprintf(fid, '	%f', e);
    fprintf(fid, '	%f', f);
    fprintf(fid, '		%f \n', g);
    fprintf(fid, '%f', h);
    fprintf(fid, '	%f \n', i);
    fprintf(fid, '%f', j);
    fprintf(fid, '		%f \n', k);
    fprintf(fid, '%f', l);
    fprintf(fid, '	%f \n', m);
    fprintf(fid, '%f', n);
    fprintf(fid, '	%f \n', o);
    fprintf(fid, '%f', p);
    fprintf(fid, '	%f', q);
    fprintf(fid, '	%f', r);
    fprintf(fid, '		%f \n', s);
    fprintf(fid, '%f', t);
    fprintf(fid, '	%f', u);

    fclose(fid);

end