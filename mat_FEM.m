function [M,S] = mat_FEM(MESH,p,Nq,siz_x,ele2dof)

S = sparse(siz_x,siz_x);
M = sparse(siz_x,siz_x);

dim_P = (p+1)*(p+2)/2;

[x,y,w] = quadratureSimplex(Nq);

M_ref = zeros(dim_P,dim_P);

for i_q = 1 : Nq
    phi_x = phi_FEM([x(i_q), y(i_q)],p);
    for i = 1:dim_P
        for j = 1:dim_P
            M_ref(i,j) = M_ref(i,j) + w(i_q)*phi_x(i)*phi_x(j);
        end
    end
end

for i_el = 1 : MESH.nElem
    B = squeeze(MESH.Tmat(i_el,:,:));
    for i_q = 1 : Nq
        x_hat = [x(i_q),y(i_q)]';
        dphi = dphi_FEM(x_hat,p);
        for i = 1 : dim_P
            iDof = ele2dof(i_el,i);
            for j = 1 : dim_P
                jDof = ele2dof(i_el,j);
                S(iDof,jDof) = S(iDof,jDof) + w(i_q)*MESH.vol(i_el)*(inv(B)'*dphi(i,:)')'*inv(B)'*dphi(j,:)';
            end
        end
    end

    for i = 1 : dim_P
        iDof = ele2dof(i_el,i);
        for j = 1 : dim_P
            jDof = ele2dof(i_el,j);
            M(iDof,jDof) = M(iDof,jDof) + M_ref(i,j)*MESH.vol(i_el);
        end
    end
    
end
