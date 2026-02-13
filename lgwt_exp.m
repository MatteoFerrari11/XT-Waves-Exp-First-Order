function [x,w] = lgwt_exp(N,T,a,b)

%Computes the nodes and weights for a Gaussian quadrature 
%weighted by the exponential function exp(-t/T) over the interval [a,b]
%
%input:
%N ---> N+1 quadrature points (exact for polynomials of degree 2N*1) 
%
%T ---> final time
%
%output:
%x --> N+1 quadrature nodes
%w --> N+1 quadrature weights

syms x
weight = exp(-x/T);
p0 = 1;
delta = zeros(N+1,1);
gamma = zeros(N,1);
delta(1) = int(x*weight*p0*p0,a,b)/int(weight*p0*p0,a,b);
p1 = x-delta(1);
for j = 1:N
   delta(j+1) = double(int(x*weight*p1*p1,a,b)/int(weight*p1*p1,a,b));
   gamma(j) = double(int(p1*p1*weight,a,b)/int(p0*p0*weight,a,b));
   p2 = (x-delta(j+1))*p1 - gamma(j)*p0;
   p0 = p1;
   p1 = p2;
end

Tri = diag(delta) + diag(sqrt(gamma),1) + diag(sqrt(gamma),-1);
[X,D] = eig(Tri);
x = eig(D);
F = T*(1-exp(-b/T));
w = X(1,:)'.^2*F;
