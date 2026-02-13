clc
clear
close all
format long

p1 = 4;
p2 = 4;
reg = 'max';
T = 4;
a = 0;
b = 1;
Nq = 16;
N1 = 2;
N2 = 3;

flag = 'no_exp';

inf_sup = zeros(p2-p1+1,N2-N1+1);

c = @(x) 1;
cont_p = 0;

for p_x = p1:p2

    cont_p = cont_p+1;
    p_t = p_x;

    if strcmp(reg,'max')
        r_x = p_x-1;
        r_t = p_t-1;
    else
        r_x = 1;
        r_t = 1;
    end

    h_x = zeros(1,N2-N1+1);

    cont = 0;

    for ii = N1:N2

        cont = cont + 1;
        N_x = 2^(ii+1);
        N_t = 2^(ii+1)*T;

        knots_t = [zeros(r_t+1,1); repelem(linspace(0,T,N_t+1),p_t-r_t)'; ones(r_t+1,1)*T]';
        knots_x = [a*ones(r_x+1,1); repelem(linspace(a,b,N_x+1),p_x-r_x)'; b*ones(r_x+1,1)]';

        t = ((0:N_t)/N_t)*T;
        x = a + ((0:N_x)/N_x)*(b-a);
        h_t = t(2)
        h_x(cont) = x(2)-x(1)

        B_t_e = mat_splines_exp(N_t,p_t,r_t,1,1,knots_t,T);
        C_t_e = mat_splines_exp(N_t,p_t,r_t,0,1,knots_t,T);
        M_t_e = mat_splines_exp(N_t,p_t,r_t,0,0,knots_t,T);

        if strcmp(flag,'no_exp')
            B_t = mat_splines(N_t,p_t,r_t,1,1,knots_t,0,T);
            C_t = mat_splines(N_t,p_t,r_t,0,1,knots_t,0,T);
            M_t = mat_splines(N_t,p_t,r_t,0,0,knots_t,0,T);
        else
            B_t = mat_splines_exp(N_t,p_t,r_t,1,1,knots_t,T);
            C_t = mat_splines_exp(N_t,p_t,r_t,0,1,knots_t,T);
            M_t = mat_splines_exp(N_t,p_t,r_t,0,0,knots_t,T);
        end
        M_x = mat_splines(N_x,p_x,r_x,0,0,knots_x,a,b);
        B_x = mat_splines_with_c(N_x,p_x,r_x,1,1,knots_x,a,b,c);

        B_t = B_t(2:end,2:end);
        M_t = M_t(2:end,2:end);
        C_t = C_t(2:end,2:end);
        B_t_e = B_t_e(2:end,2:end);
        M_t_e = M_t_e(2:end,2:end);
        C_t_e = C_t_e(2:end,2:end);
        M_x = M_x(2:end-1,2:end-1);
        B_x = B_x(2:end-1,2:end-1);


        R_Bx = chol(sparse(B_x));
        M_R = sparse(M_x) / R_Bx;
        Nx = M_R * M_R';

        sBt = sparse(B_t); sMt = sparse(M_t); sCt = sparse(C_t);
        sBt_e = sparse(B_t_e); sMt_e = sparse(M_t_e); sCt_e = sparse(C_t_e);
        sBx = sparse(B_x); sMx = sparse(M_x); sNx = sparse(Nx);

        H1_11 = kron(sBt, sMx) + kron(sMt, sBx);
        H1_22 = kron(sMt, sMx) + kron(sBt, sNx);

        H2_11 = kron(sBt, sMx);
        H2_22 = kron(sBt, sNx);

        H1 = sparse(blkdiag(H1_11, H1_22));
        H2 = sparse(blkdiag(H2_11, H2_22));

        S11 = kron(sCt_e, sBx);
        S12 = kron(sBt_e, sMx);
        S21 = S12;
        S22 = -kron(sCt_e, sMx);

        S = sparse([S11, S12; S21, S22]);

        LHS = S' * (H2 \ S);

        opts.issym = 1;
        lambda_min = eigs(LHS, H1, 1, 'smallestabs', opts);

        inf_sup(cont_p,cont) = sqrt(real(lambda_min))

    end
end