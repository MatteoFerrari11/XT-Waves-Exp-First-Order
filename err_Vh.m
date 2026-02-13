function err_norm_Vh = err_Vh(U_app_coeff,V_app_coeff,dtV_N_coeff,V_N_app_coeff,N_x,N_t,siz_x,siz_t,p_x,p_t,Nq,x,t,knots_x,knots_t,V_ex,dx_U_ex,V0,dx_U0,c)

L2_dt_U = 0;
L2_dx_U = 0;
L2_dt_VN = 0;
L2_V = 0;

norm_dt_U = 0;
norm_dx_U = 0;
norm_dt_VN = 0;
norm_V = 0;

for k_x = 1 : N_x
   
   [xs_x,ws_x] = lgwt(Nq,x(k_x),x(k_x+1));

    for k_t = 1 : N_t

        [xs_t,ws_t] = lgwt(Nq,t(k_t),t(k_t+1));

        for i_q_x = 1 : Nq

            x_val = xs_x(i_q_x);

            for i_q_t = 1:Nq

                t_val = xs_t(i_q_t);

                dx_U_h_val = 0;
                dt_U_h_val = 0;
                V_h_val = 0;
                dt_VN_h_val = 0;
                dt_VN_ex_loc = 0;

                for j_t = 1 : siz_t 

                    phi_t_val = sp_and_der(p_t,knots_t,j_t-1,t_val,0);
                    dphi_t_val = sp_and_der(p_t,knots_t,j_t-1,t_val,1);

                    for j_x = 1 : siz_x

                        phi_x_val = sp_and_der(p_x,knots_x,j_x-1,x_val,0);
                        dphi_x_val = sp_and_der(p_x,knots_x,j_x-1,x_val,1);

                        id_global = (j_t-1)*siz_x + j_x;

                        V_h_val = V_h_val + V_app_coeff(id_global) * phi_x_val * phi_t_val;
                        dx_U_h_val = dx_U_h_val + U_app_coeff(id_global) * dphi_x_val * phi_t_val;
                        dt_U_h_val = dt_U_h_val + U_app_coeff(id_global) * phi_x_val * dphi_t_val;
                        dt_VN_h_val = dt_VN_h_val + V_N_app_coeff(id_global) * dphi_x_val * dphi_t_val;
                        
                    end
                end

                %special for our problem
                aux_t_val = (25*pi^2*cos((5*pi*t_val)/2))/8;

                for j_x = 1 : siz_x

                    dphi_x_val = sp_and_der(p_x,knots_x,j_x-1,x_val,1);
                    dt_VN_ex_loc = dt_VN_ex_loc + dtV_N_coeff(j_x) * dphi_x_val * aux_t_val;

                end

                dx_U_ex_loc = dx_U_ex(t_val,x_val);
                V_ex_loc = V_ex(t_val,x_val);
                dt_U_ex_loc = V_ex(t_val,x_val);

                diff_V = V_h_val + V0(x_val) - V_ex_loc;
                diff_dt_VN = (dt_VN_h_val - dt_VN_ex_loc)*c(x_val);
                diff_dx_U = (dx_U_h_val + dx_U0(x_val) - dx_U_ex_loc)*c(x_val);
                diff_dt_U = dt_U_h_val + V0(x_val) - dt_U_ex_loc;

                w_tot = ws_x(i_q_x) * ws_t(i_q_t);

                L2_dt_U = L2_dt_U + (diff_dt_U)^2 * w_tot;
                L2_dx_U = L2_dx_U + (diff_dx_U)^2 * w_tot;
                L2_dt_VN = L2_dt_VN + (diff_dt_VN)^2 * w_tot;
                L2_V = L2_V + (diff_V)^2 * w_tot;

                norm_V = norm_V + V_ex_loc^2 * w_tot;
                norm_dt_U = norm_dt_U + dt_U_ex_loc^2 * w_tot;
                norm_dx_U = norm_dx_U + dx_U_ex_loc^2 * w_tot;
                norm_dt_VN = norm_dt_VN + dt_VN_ex_loc^2 * w_tot;
            end
        end
    end
end

err_norm_Vh = sqrt(L2_dt_U + L2_dx_U + L2_dt_VN + L2_V)/sqrt(norm_V + norm_dt_U + norm_dx_U + norm_dt_VN);