function result=local_FE_function_1D(x,uh_local_vec,vertices,basis_type,basis_der_x)
% local_FE_function_1D 在一个单元内计算有限元函数（或其 basis_der_x 阶导数）
%   x             : 求值点
%   uh_local_vec  : 该单元上的局部自由度（节点值）
%   vertices      : 单元两个端点坐标
%   basis_type    : 101 一次元 / 102 二次元
%   basis_der_x   : 基函数求导阶数（0 为函数本身，1 为一阶导数）

result=0;

number_of_local_basis_fun=length(uh_local_vec);

for k=1:number_of_local_basis_fun

    result=result+uh_local_vec(k)*FE_basis_local_fun_1D(x,vertices,basis_type,k,basis_der_x);

end

end
