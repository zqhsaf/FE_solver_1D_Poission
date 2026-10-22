function Hs_error=compute_Hs_error(analytic_fun,basis_der_x,number_of_elements,P,T,Tb,solution,basis_type)
% compute_Hs_error 计算 H^s 误差范数（逐单元高斯积分后求和再开平方）
%
%   analytic_fun       : 精确解函数，调用形式 analytic_fun(x, basis_der_x)
%   basis_der_x==0     : 返回 L2 范数误差 ||u-u_h||_0
%   basis_der_x==1     : 返回 H1 半范数误差 |u-u_h|_1
%   number_of_elements : 单元数
%   P,T                : 网格端点坐标矩阵、单元端点拓扑矩阵
%   Tb                 : 全局基函数编号矩阵
%   solution           : 有限元节点解
%   basis_type         : 101 一次元 / 102 二次元

Hs_error=0;

Gpn=5;

for n=1:number_of_elements

    vertices=P(:,T(:,n));

    uh_local_vec=solution(Tb(:,n));

    [Gauss_weights,Gauss_nodes]=generate_Gauss_local(vertices,Gpn);

    int_value=Gauss_int_error_1D(analytic_fun,Gauss_weights,Gauss_nodes,vertices,uh_local_vec,basis_type,basis_der_x);

    Hs_error=Hs_error+int_value;

end

Hs_error=sqrt(Hs_error);

end
