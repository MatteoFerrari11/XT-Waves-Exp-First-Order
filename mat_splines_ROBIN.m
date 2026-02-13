function K = mat_splines_ROBIN(Nt,p,r,k1,k2,knots,T1,T2)

siz = Nt*(p-r)+r+1;
K = zeros(siz,siz);

for jj = 1 : siz
    for ii = 1 : siz
        K(jj,ii) = sp_and_der(p,knots,ii-1,T2-10*eps,k1)*sp_and_der(p,knots,jj-1,T2-10*eps,k2);
    end
end
