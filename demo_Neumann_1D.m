%% Example 3: 混合边界条件 u(0)=0 (Dirichlet), u'(1)=cos(1)-sin(1) (Neumann)
%  -d/dx(e^x du/dx) = -e^x[cos(x)-2sin(x)-xcos(x)-xsin(x)], 0<=x<=1
%  精确解 u(x)=x cos(x)
%
%  Neumann 边界处理(弱形式):
%    ∫ c u'v' dx = ∫ f v dx + r_b c(b) v(b)
%  只需在载荷向量 x=b 对应节点分量上加 r_b*c(b)
clear; clc;

rb=cos(1)-sin(1);          % 右端 Neumann 边界值 u'(1)
N=4;

%% 1. 一次元(101)求解
u1=FE_solver_1D_Poission(0,1,101,101,N,0,rb,'function_u_exact',1,2);

Pb1=0:1/N:1;
fprintf('\n===== 一次元 N=%d 节点结果 =====\n',N);
disp(table(Pb1',u1,Pb1'.*cos(Pb1'),u1-Pb1'.*cos(Pb1'), ...
    'VariableNames',{'x','u_h','u_exact','error'}));

%% 2. 二次元(102)求解
u2=FE_solver_1D_Poission(0,1,102,102,N,0,rb,'function_u_exact',1,2);

Pb2=0:1/(2*N):1;
fprintf('===== 二次元 N=%d 节点结果 =====\n',N);
disp(table(Pb2',u2,Pb2'.*cos(Pb2'),u2-Pb2'.*cos(Pb2'), ...
    'VariableNames',{'x','u_h','u_exact','error'}));

%% 3. 直接验证 Neumann 条件：比较 u_h'(1) 与 r_b
[P,T]=generate_PT(0,1,N);
[~,Tb1]=generate_PbTb(0,1,N,101);
[~,Tb2]=generate_PbTb(0,1,N,102);
vertices=P(:,T(:,N));
du1_b=local_FE_function_1D(1,u1(Tb1(:,N)),vertices,101,1);
du2_b=local_FE_function_1D(1,u2(Tb2(:,N)),vertices,102,1);
fprintf('Neumann 条件检验: u_h''(1) 一次元=%.8f, 二次元=%.8f, 精确值=%.8f\n', ...
    du1_b,du2_b,rb);

%% 4. 收敛阶验证
fprintf('\n===== 一次元收敛阶 (理论: 节点/L2 为 O(h^2), H1 为 O(h)) =====\n');
print_order(101,rb);
fprintf('\n===== 二次元收敛阶 (理论: L2 为 O(h^3), H1 为 O(h^2)) =====\n');
print_order(102,rb);

function print_order(basis_type,rb)

fprintf('%4s %12s %8s %12s %8s %12s %8s\n', ...
    'N','max节点误差','阶','L2误差','阶','H1误差','阶');

em_prev=0; eL_prev=0; eH_prev=0;

for N=[4 8 16 32]

    [P,T]=generate_PT(0,1,N);
    [Pb,Tb]=generate_PbTb(0,1,N,basis_type);

    u=FE_solver_1D_Poission(0,1,basis_type,basis_type,N,0,rb,[],1,2);

    em=compute_max_FE_nodes('function_u_exact',Pb,u);
    eL=compute_Hs_error('analytic_solution',0,N,P,T,Tb,u,basis_type);
    eH=compute_Hs_error('analytic_solution',1,N,P,T,Tb,u,basis_type);

    if em_prev==0
        fprintf('%4d %12.3e %8s %12.3e %8s %12.3e %8s\n',N,em,'-',eL,'-',eH,'-');
    else
        fprintf('%4d %12.3e %8.2f %12.3e %8.2f %12.3e %8.2f\n', ...
            N,em,log2(em_prev/em),eL,log2(eL_prev/eL),eH,log2(eH_prev/eH));
    end

    em_prev=em; eL_prev=eL; eH_prev=eH;
end

end
