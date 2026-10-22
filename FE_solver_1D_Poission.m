function solution=FE_solver_1D_Poission(left,right,basis_type_trial,basis_type_test,number_of_element,value1,value2,exact_fun,left_type,right_type,q_value1,q_value2)
%basis_type_trial==101:1D linear
%basis_type_trial==102:1D quadratic
%exact_fun 为可选参数：传入精确解函数时自动计算并打印节点最大误差
%left_type/right_type 为可选参数：边界类型 1=Dirichlet，2=Neumann，3=Robin，默认两端均为 Dirichlet
%value1/value2 为左/右端边界值（Dirichlet 给函数值；Neumann 给导数值 r；Robin 给 p）
%q_value1/q_value2 为左/右端 Robin 系数 q（仅在边界类型为 3 时使用，默认 0）

N=number_of_element;

[P,T]=generate_PT(left,right,N);

[Pb,Tb]=generate_PbTb(left,right,N,basis_type_trial);

if basis_type_trial==101
    Pb_trial=Pb;
    Tb_trial=Tb;
    number_of_local_basis_fun_trial=2;
elseif basis_type_trial==102
    Pb_trial=Pb;
    Tb_trial=Tb;
    number_of_local_basis_fun_trial=3;
        
end

if basis_type_test==101
    Pb_test=Pb;
    Tb_test=Tb;
    number_of_local_basis_fun_test=2;
elseif basis_type_test==102
    Pb_test=Pb;
    Tb_test=Tb;
    number_of_local_basis_fun_test=3;
        
end

if nargin<9 || isempty(left_type)
    left_type=1;
end

if nargin<10 || isempty(right_type)
    right_type=1;
end

if nargin<11 || isempty(q_value1)
    q_value1=0;
end

if nargin<12 || isempty(q_value2)
    q_value2=0;
end

boundarynodes=generate_boundarynodes(Pb_test,left_type,right_type);

matrix_size=[length(Pb_test),length(Pb_trial)];
number_of_elements=size(T,2);
basis_der_x_trial=1;
basis_der_x_test=1;

A=assemble_matrix_1D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test);

vector_size=[length(Pb_test),1];

b=assemble_vector_1D('function_f',vector_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,0);

dirichlet_tag=1;

g=generate_boundary_value(Pb_test,boundarynodes,value1,value2);

gq=sparse(1,length(Pb_test));
gq(boundarynodes(2,1))=q_value1;
gq(boundarynodes(2,2))=q_value2;

[A,b]=treat_Robin_boundary_1D('function_c',A,b,boundarynodes,Pb_test,g,gq);

b=treat_Neumann_boundary_1D('function_c',b,boundarynodes,Pb_test,g);

[A,b]=treat_Dirichlet_boundary(A,b,boundarynodes, Pb_test, g, dirichlet_tag);

solution=A\b;

