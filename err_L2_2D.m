function [L2_U,L2_V] = err_L2(U_app_coeff,V_app_coeff,mesh_x,siz_t,siz_x,p_x,Nq_x,p_t,Nq_t,N_t,t,knots_t,U_ex,V_ex,ele2dof)

L2_U = 0;
L2_V = 0;

dim_p_x = (p_x+1)*(p_x+2)/2;


for elem_x = 1 : mesh_x.nElem
    B = squeeze(mesh_x.Tmat(elem_x,:,:));
    x0 = mesh_x.coord(elem_x,:,1)';
    [xs_x,ys_x,ws_x] = quadratureSimplex(Nq_x);

    for k_t = 1 : N_t
        [xs_t,ws_t] = lgwt(Nq_t,t(k_t),t(k_t+1));

        for i_q_x = 1 : Nq_x
            x_hat = [xs_x(i_q_x),ys_x(i_q_x)]';
            x_real = x0 + B*x_hat;
            phi_x = phi_FEM(x_hat,p_x);

            for i_q_t = 1:Nq_t
                t_val = xs_t(i_q_t);

                U_h_val = 0;
                V_h_val = 0;

                for j_t = 1 : siz_t 

                    phi_t_val = sp_and_der(p_t,knots_t,j_t-1,t_val,0);

                    if abs(phi_t_val) < 1e-14
                        continue; 
                    end

                    for i_x = 1 : dim_p_x
                        iDof = ele2dof(elem_x,i_x);
                        idx_global = (j_t-1)*siz_x + iDof;

                        U_h_val = U_h_val + U_app_coeff(idx_global) * phi_x(i_x) * phi_t_val;
                        V_h_val = V_h_val + V_app_coeff(idx_global) * phi_x(i_x) * phi_t_val;
                    end
                end

                diff_U = U_h_val - U_ex(x_real(1),x_real(2),t_val);
                diff_V = V_h_val - V_ex(x_real(1),x_real(2),t_val);

                w_tot = ws_x(i_q_x) * ws_t(i_q_t) * mesh_x.vol(elem_x);
                L2_U = L2_U + (diff_U)^2 * w_tot;
                L2_V = L2_V + (diff_V)^2 * w_tot;
            end
        end
    end
end

L2_U = sqrt(L2_U);
L2_V = sqrt(L2_V);