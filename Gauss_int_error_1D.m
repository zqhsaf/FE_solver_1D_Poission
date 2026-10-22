function int_value=Gauss_int_error_1D(analytic_fun,Gauss_weights,Gauss_nodes,vertices,uh_local_vec,basis_type,basis_der_x)
% Gauss_int_error_1D 在一个单元上用高斯积分计算
%   (u^(basis_der_x) - u_h^(basis_der_x))^2 的积分
%
%   analytic_fun : 精确解函数，调用形式 analytic_fun(x, basis_der_x)
%   basis_der_x==0: L2 误差被积函数；basis_der_x==1: H1 半范数误差被积函数

Gpn=length(Gauss_weights);

int_value=0;

for k=1:Gpn

    u_exact_der=feval(analytic_fun,Gauss_nodes(k),basis_der_x);

    u_h_der=local_FE_function_1D(Gauss_nodes(k),uh_local_vec,vertices,basis_type,basis_der_x);

    int_value=int_value+Gauss_weights(k)*(u_exact_der-u_h_der)^2;

end

end
