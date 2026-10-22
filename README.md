# 一维二阶椭圆方程有限元通用求解器（MATLAB）

使用一维线性、二次有限元求解含 Dirichlet / Neumann / Robin 边界条件的二阶椭圆方程：

$$
-\frac{d}{dx}\left(c(x)\frac{du}{dx}\right) = f(x),\qquad a\le x\le b
$$

## 算例

$$
-\frac{d}{dx}\left(e^x\frac{du}{dx}\right)
=-e^x\big[\cos x-2\sin x-x\cos x-x\sin x\big],\quad 0\le x\le 1
$$

- 精确解：$u(x)=x\cos x$
- 边界条件：$u(0)=0$，右端可为 Dirichlet $u(1)=\cos 1$、Neumann $u'(1)=\cos 1-\sin 1$ 或 Robin $u'(1)+q_b u(1)=p_b$

## 功能特性

- 一次元（`basis_type=101`）与二次元（`basis_type=102`）
- 三类边界条件：Dirichlet（类型 1）、Neumann（类型 2）、Robin（类型 3）
- 5 点高斯积分装配刚度矩阵与载荷向量
- 节点最大误差、L2 误差、H1 半范数误差计算
- 收敛阶实测：一次元 $O(h^2)$（H1 为 $O(h)$），二次元 $O(h^3)$（H1 为 $O(h^2)$，节点处呈 $O(h^4)$ 超收敛）

## 快速开始

```matlab
% 两端 Dirichlet，一次元
u = FE_solver_1D_Poission(0,1,101,101,N,0,cos(1),'function_u_exact');

% 两端 Dirichlet，二次元
u = FE_solver_1D_Poission(0,1,102,102,N,0,cos(1),'function_u_exact');

% 左 Dirichlet + 右 Neumann：u'(1)=cos(1)-sin(1)
rb = cos(1)-sin(1);
u = FE_solver_1D_Poission(0,1,102,102,N,0,rb,'function_u_exact',1,2);

% 左 Dirichlet + 右 Robin：u'(1)+q_b*u(1)=p_b
qb = 1;
pb = 2*cos(1)-sin(1);
u = FE_solver_1D_Poission(0,1,102,102,N,0,pb,'function_u_exact',1,3,0,qb);
```

演示脚本（MATLAB 中直接输入文件名运行）：

| 脚本 | 内容 |
|---|---|
| `demo_quadratic_1D.m` | 二次元求解与收敛阶验证 |
| `demo_Neumann_1D.m` | Neumann 混合边界算例 |
| `demo_Robin_1D.m` | Robin 混合边界算例 |

## 调用参数

```matlab
FE_solver_1D_Poission(left,right,basis_type_trial,basis_type_test,N, ...
                      value1,value2,exact_fun, ...
                      left_type,right_type,q_value1,q_value2)
```

- `value1/value2`：左/右端边界值（Dirichlet 给函数值，Neumann 给 $r$，Robin 给 $p$）
- `exact_fun`：可选，精确解函数名，传入后自动打印节点最大误差
- `left_type/right_type`：可选，1=Dirichlet，2=Neumann，3=Robin，默认 1
- `q_value1/q_value2`：可选，Robin 系数 $q$，默认 0

## 文件结构

| 文件 | 功能 |
|---|---|
| `FE_solver_1D_Poission.m` | 主求解程序 |
| `generate_PT.m` / `generate_PbTb.m` | 网格端点矩阵 P、T 与节点/基函数矩阵 Pb、Tb |
| `generate_boundarynodes.m` | 边界节点信息 |
| `generate_boundary_value.m` | 边界值向量 |
| `assemble_matrix_1D.m` / `assemble_vector_1D.m` | 刚度矩阵、载荷向量装配 |
| `Gauss_quad_1D_trial_test_matrix.m` / `..._vector.m` | 单元矩阵/向量高斯积分 |
| `generate_Gauss_local.m` | 区间高斯点与权重 |
| `FE_basis_local_fun_1D.m` | 单元局部基函数及其导数 |
| `treat_Dirichlet_boundary.m` | Dirichlet 边界处理 |
| `treat_Neumann_boundary_1D.m` | Neumann 边界处理 |
| `treat_Robin_boundary_1D.m` | Robin 边界处理 |
| `compute_max_FE_nodes.m` | 节点最大误差 |
| `compute_Hs_error.m` / `Gauss_int_error_1D.m` / `local_FE_function_1D.m` | L2/H1 误差 |
| `function_c.m` | 系数函数 $c(x)=e^x$ |
| `function_f.m` | 右端项 $f(x)$ |
| `function_u_exact.m` / `analytic_solution.m` | 精确解及其导数 |
