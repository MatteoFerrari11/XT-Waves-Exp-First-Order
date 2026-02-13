function dphi = dphi_FEM(x,p)

switch p
    case(1)
        dphi = [-1,-1; 1,0; 0 1];
    case(2)
        dphi_dx1 = [
            4*x(1) + 4*x(2) - 3;
            4*x(1) - 1;
            0;
            4 - 8*x(1) - 4*x(2);
            4*x(2);
            -4*x(2)
            ];

        dphi_dx2 = [
            4*x(1) + 4*x(2) - 3;
            0;
            4*x(2) - 1;
            -4*x(1);
            4*x(1);
            4 - 4*x(1) - 8*x(2)
            ];

        dphi = [dphi_dx1, dphi_dx2];

        bdphi = dphi;
        dphi(4,:) = bdphi(5,:);
        dphi(6,:) = bdphi(4,:);
        dphi(5,:) = bdphi(6,:);

    case(3)
        
        L1 = 1 - x(1) - x(2);
        L2 = x(1);
        L3 = x(2);

        dL_v = @(L) 0.5 * (27*L^2 - 18*L + 2);

        dL_edge_i = @(Li, Lj) 4.5 * Lj * (6*Li - 1);
        dL_edge_j = @(Li, Lj) 4.5 * Li * (3*Li - 1);

        dphi_dx = zeros(10, 1);
        dphi_dy = zeros(10, 1);

        val = dL_v(L1);
        dphi_dx(1) = -val; dphi_dy(1) = -val;

        val = dL_v(L2);
        dphi_dx(2) = val;  dphi_dy(2) = 0;

        val = dL_v(L3);
        dphi_dx(3) = 0;    dphi_dy(3) = val;

        d_L1 = dL_edge_i(L1, L2);
        d_L2 = dL_edge_j(L1, L2);
        dphi_dx(4) = d_L2 - d_L1;
        dphi_dy(4) = -d_L1;

        d_L2 = dL_edge_i(L2, L1);
        d_L1 = dL_edge_j(L2, L1);
        dphi_dx(5) = d_L2 - d_L1;
        dphi_dy(5) = -d_L1;

        d_L2 = dL_edge_i(L2, L3);
        d_L3 = dL_edge_j(L2, L3);
        dphi_dx(6) = d_L2;
        dphi_dy(6) = d_L3;

        d_L3 = dL_edge_i(L3, L2);
        d_L2 = dL_edge_j(L3, L2);
        dphi_dx(7) = d_L2;
        dphi_dy(7) = d_L3;

        d_L3 = dL_edge_i(L3, L1);
        d_L1 = dL_edge_j(L3, L1);
        dphi_dx(8) = -d_L1;
        dphi_dy(8) = d_L3 - d_L1;

        d_L1 = dL_edge_i(L1, L3);
        d_L3 = dL_edge_j(L1, L3);
        dphi_dx(9) = -d_L1;
        dphi_dy(9) = d_L3 - d_L1;

        term = 27;
        d_L1 = term * L2 * L3;
        d_L2 = term * L1 * L3;
        d_L3 = term * L1 * L2;

        dphi_dx(10) = d_L2 - d_L1;
        dphi_dy(10) = d_L3 - d_L1;

        dphi = [dphi_dx, dphi_dy];

        bdphi = dphi;
        dphi(4:5,:) = bdphi(6:7,:);
        dphi(8:9,:) = bdphi(4:5,:);
        dphi(6:7,:) = bdphi(8:9,:);

end