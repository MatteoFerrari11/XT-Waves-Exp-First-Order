function [L2_U,L2_V] = err_L2_1D(U_app_coeff,V_app_coeff,N_x,N_t,siz_x,siz_t,p_x,p_t,Nq,x,t,knots_x,knots_t,U_ex,V_ex,U0,V0)

L2_U = 0;
L2_V = 0;
norm_U = 0;
norm_V = 0;

for k_x = 1 : N_x
   
   [xs_x,ws_x] = lgwt(Nq,x(k_x),x(k_x+1));

    for k_t = 1 : N_t

        [xs_t,ws_t] = lgwt(Nq,t(k_t),t(k_t+1));

        for i_q_x = 1 : Nq

            x_val = xs_x(i_q_x);

            for i_q_t = 1:Nq

                t_val = xs_t(i_q_t);

                U_h_val = 0;
                V_h_val = 0;

                for j_t = 1 : siz_t 

                    phi_t_val = sp_and_der(p_t,knots_t,j_t-1,t_val,0);
                    dphi_t_val = sp_and_der(p_t,knots_t,j_t-1,t_val,1);

                    if abs(phi_t_val) < 1e-14
                        continue; 
                    end

                    for j_x = 1 : siz_x

                        phi_x_val = sp_and_der(p_x,knots_x,j_x-1,x_val,0);
                        dphi_x_val = sp_and_der(p_x,knots_x,j_x-1,x_val,1);

                        id_global = (j_t-1)*siz_x + j_x;

                        U_h_val = U_h_val + U_app_coeff(id_global) * phi_x_val * phi_t_val;
                        V_h_val = V_h_val + V_app_coeff(id_global) * phi_x_val * phi_t_val;
                        
                    end
                end

                U_ex_loc = U_ex(t_val,x_val);
                V_ex_loc = V_ex(t_val,x_val);

                diff_U = U_h_val - U_ex_loc + U0(x_val);
                diff_V = V_h_val - V_ex_loc + V0(x_val);

                w_tot = ws_x(i_q_x) * ws_t(i_q_t);
                L2_U = L2_U + (diff_U)^2 * w_tot;
                norm_U = norm_U + U_ex_loc^2 * w_tot;
                L2_V = L2_V + (diff_V)^2 * w_tot;
                norm_V = norm_V + V_ex_loc^2 * w_tot;
            end
        end
    end
end

L2_U = sqrt(L2_U/norm_U);
L2_V = sqrt(L2_V/norm_V);