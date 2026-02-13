% Accompanying code for the Paper:
% "Inf-sup stable space--time discretization of the wave equation based on a first-order-in-time variational formulation"
%
% Authors: M. Ferrari and I. Perugia and E. Zampa
%
% This code computes an approximate solution of the wave problem:
%
%     d^2_t U(x,t) - div (c(x) nabla U(x,t)) = F(x,t),  x,t in Q_T = Omega x (0,T)
%     U(x,0) = U0(x),  dt U(x,0) = V0(x),           x   in Omega
%     U(x,t) = 0                                    x,t   in partial Omega x (0,T)
%
%  Omega = (Lx1,Lx2) x (Ly1,Ly2)
%
% in a first order formulation
%
%     d_t V(x,t) - div (c(x) nabla U(x,t)) = F(x,t),    x,t in Q_T
%     d_t U(x,t) - V(x,t) = 0,                      x,t in Q_T
%     U(x,0) = U0(x),  V(x,0) = V0(x),              x in (0,L)
%     U(t,0) = V(t,0) = 0,  V(t,L) = U(t,L) = 0,    t in (0,T)
%
% using the variational formulation:
%
% find U in H^1(0,T;H_0^1(Omega)),  V in H^1(0,T;L^2(Omega)) with U(x,0)=V(x,0)=0
% such that
%
%     (d_t V, lambda)_{L_e^2(Q_T)} + (c nabla U, nabla lambda)_{L_e^2(Q_T)} = (F,lambda)_{L^2_e(Q_T)} - (c nabla U0, nabla lambda)_{L^2_e(Q_T)}
%     (d_t U, chi)_{L_e^2(Q_T)} - (V, chi)_{L_e^2(Q_T)} = (V0,chi)_{L^2_e(Q_T)
%
% for all lambda in L^2(0,T;H_0^1(Omega)) and chi in L^2(Q_T)
%
% Here, (U,V)_{L^2_e(Q_T)} = (U,V e^{-./T})_{L^2(Q_T)} 
%
% The code employs arbitrary regularity splines in time and piecewise
% continuous polynomials on a simplicial mesh in space

% and calculates errors in the norm L^2 of U, V
% Note: this code is far from being optimized

clc
clear
close all
format long

p_x = 3;
dim_p_x = (p_x+1)*(p_x+2)/2;
dim_p_x_3 = (p_x-1)*(p_x-2)/2;
L_x1 = 0;
L_x2 = 1;
L_y1 = 0;
L_y2 = 1;

p_t = 3;
r_t = 2;
T = 1;

Nq_t = 8;
Nq_x = 7;

N1 = 1;
N2 = 2;

err_U = zeros(1,N2-N1+1);
err_V = zeros(1,N2-N1+1);

U = @(x,y,t) sin(pi*x)*sin(pi*y)*sin(t*x*y)^2;
F = @(x,y,t) 2*pi^2*sin(pi*x)*sin(pi*y) + 2*t^2*x^2*sin(pi*x)*sin(pi*y) + 2*t^2*y^2*sin(pi*x)*sin(pi*y) - 2*x^2*y^2*sin(pi*x)*sin(pi*y) - 2*pi^2*cos(t*x*y)^2*sin(pi*x)*sin(pi*y) - 4*t^2*x^2*cos(t*x*y)^2*sin(pi*x)*sin(pi*y) - 4*t^2*y^2*cos(t*x*y)^2*sin(pi*x)*sin(pi*y) + 4*x^2*y^2*cos(t*x*y)^2*sin(pi*x)*sin(pi*y) - 4*t*x*pi*cos(t*x*y)*sin(t*x*y)*cos(pi*y)*sin(pi*x) - 4*t*y*pi*cos(t*x*y)*sin(t*x*y)*cos(pi*x)*sin(pi*y);
dtU =  @(x,y,t) x*y*sin(pi*x)*sin(pi*y)*sin(2*t*x*y);

cont = 0;

for ii = N1:N2

    tic 

    cont = cont + 1;
    N_x = 2^(ii+1);
    N_t = 2^(ii+1);

    knots_t = [zeros(r_t+1,1); repelem(linspace(0,T,N_t+1),p_t-r_t)'; ones(r_t+1,1)*T]';
    mesh_x = generateSquareMesh(L_x1,L_x2,L_y1,L_y2,N_x,N_x,[0 0]);
    mesh_x = addMeshStructures(mesh_x);
    siz_x = mesh_x.nNode + (p_x-1)*mesh_x.nFace + dim_p_x_3*mesh_x.nElem;
    
    ele2dof = zeros(mesh_x.nElem,dim_p_x);
    bndFlagDof = zeros(siz_x,1);
    for i_elem = 1 : mesh_x.nElem
        ele2dof(i_elem,1:3) = mesh_x.tri(i_elem,:);
        for i = 1 : 3
            iDof = mesh_x.tri(i_elem,i);
            if mesh_x.bndFlagNode(iDof) == 1
                bndFlagDof(iDof) = 1;
            end
        end
        for i_loc_face = 1 : 3
            i_face = mesh_x.triFace(i_elem,i_loc_face);
            for i = 1 : p_x - 1
                if mesh_x.element2faceSign(i_elem,i_loc_face) == 1
                    iDof = mesh_x.nNode + (p_x-1)*(i_face-1) + i;
                else
                    iDof = mesh_x.nNode + (p_x-1)*(i_face-1) + p_x-i;
                end
                bndFlagDof(iDof) = mesh_x.bndFlagFace(i_face);
                ele2dof(i_elem,3+(p_x-1)*(i_loc_face-1)+i) = iDof;
            end
        end
        for i = 1 : dim_p_x_3
            iDof = mesh_x.nNode + (p_x-1)*mesh_x.nFace + dim_p_x_3*(i_elem-1)+i;
            ele2dof(i_elem,3*p_x+i) = iDof;
        end
    end

    t = ((0:N_t)/N_t)*T;

    h_t = t(2);
    h_x = sqrt((L_x1-L_x2)^2+(L_y1-L_y2)^2)/N_x; 

    B_t = mat_splines_exp(N_t,p_t,r_t,1,1,knots_t,T);
    C_t = mat_splines_exp(N_t,p_t,r_t,0,1,knots_t,T);

    [M_x,B_x] = mat_FEM(mesh_x,p_x,Nq_x,siz_x,ele2dof);

    siz_t = N_t*(p_t-r_t)+r_t+1;
    F_vec_2 = zeros(siz_x*siz_t,1);

    cont_v = [];
    for elem_t = 2 : siz_t
        cont_v = [cont_v elem_t];
    end
    B_t = B_t(cont_v,cont_v);
    C_t = C_t(cont_v,cont_v);

    cont_v = [];
    for elem_x = 1 : siz_x
        if bndFlagDof(elem_x) ~= 1
            cont_v = [cont_v elem_x];
        end
    end
    M_x = M_x(cont_v,cont_v);
    B_x = B_x(cont_v,cont_v);

    time_matrices = toc

    [xs_x,ys_x,ws_x] = quadratureSimplex(Nq_x);
    F_vec_1 = zeros(siz_x*siz_t,1);

    tic
    for elem_x = 1 : mesh_x.nElem
        B = squeeze(mesh_x.Tmat(elem_x,:,:));
        x0 = mesh_x.coord(elem_x,:,1)';
        for i_q_x = 1 : Nq_x
            x_hat = [xs_x(i_q_x),ys_x(i_q_x)]';
            x_real = x0 + B*x_hat;
            phi_x = phi_FEM(x_hat,p_x);
            for elem_t = 1 : siz_t
                for k_t = max(1,floor(elem_t/(p_t-r_t))-(p_t+elem_t)):min(elem_t,N_t)
                    [xs_t,ws_t] = lgwt(Nq_t,t(k_t),t(k_t+1));
                    for i_q_t = 1 : Nq_t
                        phi_t = sp_and_der(p_t,knots_t,elem_t-1,xs_t(i_q_t),1)*exp(-xs_t(i_q_t)/T);
                        for i = 1 : dim_p_x
                            iDof = ele2dof(elem_x,i);
                            F_vec_1((elem_t-1)*siz_x+iDof) = F_vec_1((elem_t-1)*siz_x+iDof)...
                                + phi_t*ws_t(i_q_t)*F(x_real(1),x_real(2),xs_t(i_q_t))*phi_x(i)*ws_x(i_q_x)*mesh_x.vol(elem_x);
                        end
                    end
                end
            end
        end
    end

    cont_v = [];
    for elem_t = 1 : siz_t
        for elem_x = 1 : siz_x
            if ~(bndFlagDof(elem_x) == 1 || elem_t == 1)
                cont_v = [cont_v (elem_t-1)*siz_x+elem_x];
            end
        end
    end

    F_vec_1 = F_vec_1(cont_v);
    F_vec_2 = F_vec_2(cont_v);

    RHS = [F_vec_1 ; F_vec_2];

    time_F = toc


    tic
    [U_app_coeff,V_app_coeff] = solver(sparse(B_t),sparse(C_t),sparse(M_x),sparse(B_x),F_vec_1);
    U_app_coeff = real(U_app_coeff);
    V_app_coeff = real(V_app_coeff);
    time_solver = toc

    % tic
    % S = [kron(sparse(C_t),sparse(B_x)) kron(sparse(B_t),sparse(M_x)); ...
    %       kron(sparse(B_t),sparse(M_x)) -kron(sparse(C_t),sparse(M_x))];
    % SOL = S \ RHS;
    % U_app_coeff_2 = SOL(1:length(SOL)/2);
    % V_app_coeff_2 = SOL(length(SOL)/2+1:end);
    % time_system_without_solver = toc

    
    U_app_coeff_w_b = zeros(siz_x*siz_t,1);
    V_app_coeff_w_b = zeros(siz_x*siz_t,1);
    cont_v = 0;
    for elem_t = 1 : siz_t
        for elem_x = 1 : siz_x
            if elem_t == 1
                U_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = 0;
                V_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = 0;
            elseif bndFlagDof(elem_x) == 1
                U_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = 0;
                V_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = 0;
            else
                cont_v = cont_v+1;
                U_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = U_app_coeff(cont_v);
                V_app_coeff_w_b((elem_t-1)*siz_x+elem_x) = V_app_coeff(cont_v);
            end
        end
    end

    U_app_coeff = U_app_coeff_w_b;
    V_app_coeff = V_app_coeff_w_b;

    tic
    [err_U(cont),err_V(cont)] = err_L2_2D(U_app_coeff,V_app_coeff,mesh_x,siz_t,siz_x,p_x,Nq_x,p_t,Nq_t,N_t,t,knots_t,U,dtU,ele2dof)
    time_error = toc

end

if N1 < N2 && N_x == N_t
    plot_error({}, err_U,T./2.^(N1+1:N2+1),1);
    plot_error({'L2 error U','L2 error V'}, err_V, T./2.^(N1+1:N2+1),1);
    title(['Errors in 2D for splines with px=' num2str(p_x) ' rx=0 and pt=' num2str(p_t) ' rt=' num2str(r_t)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end