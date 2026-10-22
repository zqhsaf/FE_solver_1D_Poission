%% 一维二次元(basis_type=102)求解二阶椭圆方程
%  -d/dx(e^x du/dx) = -e^x[cos(x)-2sin(x)-xcos(x)-xsin(x)], 0<=x<=1
%  u(0)=0, u(1)=cos(1), 精确解 u(x)=x cos(x)
clear; clc;

%% 1. N=4 个二次元单元求解
N=4;
u=FE_solver_1D_Poission(0,1,102,102,N,0,cos(1),'function_u_exact');

Pb=0:1/(2*N):1;
u_exact=Pb.*cos(Pb);

fprintf('\n===== 二次元 N=%d 各节点结果 =====\n',N);
disp(table(Pb',u,u_exact',u-u_exact', ...
    'VariableNames',{'x','u_h','u_exact','error'}));

%% 2. 收敛阶验证（N 依次加倍；二次元理论阶：节点最大/L2 为 O(h^3)，H1 为 O(h^2)）
fprintf('===== 二次元收敛阶 =====\n');
fprintf('%4s %12s %8s %12s %8s %12s %8s\n', ...
    'N','max节点误差','阶','L2误差','阶','H1误差','阶');

e_max_prev=0; e_L2_prev=0; e_H1_prev=0;

for N=[4 8 16 32]

    [P,T]=generate_PT(0,1,N);
    [Pb,Tb]=generate_PbTb(0,1,N,101);

    u=FE_solver_1D_Poission(0,1,101,101,N,0,cos(1));

    e_max=compute_max_FE_nodes('function_u_exact',Pb,u);
    e_L2 =compute_Hs_error('analytic_solution',0,N,P,T,Tb,u,101);
    e_H1 =compute_Hs_error('analytic_solution',1,N,P,T,Tb,u,101);

    if e_max_prev==0
        fprintf('%4d %12.3e %8s %12.3e %8s %12.3e %8s\n', ...
            N,e_max,'-',e_L2,'-',e_H1,'-');
    else
        fprintf('%4d %12.3e %8.2f %12.3e %8.2f %12.3e %8.2f\n', ...
            N,e_max,log2(e_max_prev/e_max), ...
            e_L2,log2(e_L2_prev/e_L2), ...
            e_H1,log2(e_H1_prev/e_H1));
    end

    e_max_prev=e_max; e_L2_prev=e_L2; e_H1_prev=e_H1;
end
