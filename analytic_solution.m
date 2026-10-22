function result=analytic_solution(x,basis_der_x)
% 本题精确解 u(x)=x*cos(x) 及其导数
%   basis_der_x=0: u(x)   = x cos x
%   basis_der_x=1: u'(x)  = cos x - x sin x
%   basis_der_x=2: u''(x) = -2 sin x - x cos x

if basis_der_x==0

    result=x.*cos(x);

elseif basis_der_x==1

    result=cos(x)-x.*sin(x);

elseif basis_der_x==2

    result=-2*sin(x)-x.*cos(x);

else

    error('analytic_solution: 不支持的导数阶数');

end

end
