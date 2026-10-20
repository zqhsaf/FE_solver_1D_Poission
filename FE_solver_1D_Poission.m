function solution=FE_solver_1D_Poission(left,right,basis_type_trial,basis_type_test,number_of_element,value1,value2,exact_fun)
%basis_type_trial==101:1D linear
%exact_fun 为可选参数：传入精确解函数时自动计算并打印节点最大误差

N=number_of_element;

[P,T]=generate_PbTb(left,right,N,basis_type_trial);

if basis_type_trial==101
    Pb_trial=P;
    Tb_trial=T;
    number_of_local_basis_fun_trial=2;
elseif basis_type_trial==102
        
end

if basis_type_test==101
    Pb_test=P;
    Tb_test=T;
    number_of_local_basis_fun_test=2;
elseif basis_type_test==102
        
end

left_type=1;
right_type=1;

boundarynodes=generate_boundarynodes(Pb_test,left_type,right_type);

matrix_size=[N+1,N+1];
number_of_elements=size(T,2);
basis_der_x_trial=1;
basis_der_x_test=1;

A=assemble_matrix_1D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test);

vector_size=[N+1,1];

b=assemble_vector_1D('function_f',vector_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,0);

dirichlet_tag=1;

g=generate_boundary_value(Pb_test,boundarynodes,value1,value2);

[A,b]=treat_Dirichlet_boundary(A,b,boundarynodes, Pb_test, g, dirichlet_tag);

solution=A\b;

if nargin>=8 && ~isempty(exact_fun)
    max_error=compute_max_FE_nodes(exact_fun,Pb_test,solution);
    fprintf('节点最大误差 max|u_h-u| = %.3e (N=%d)\n',max_error,N);
end