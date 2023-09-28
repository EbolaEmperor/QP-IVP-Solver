# 四精度浮点数IVP求解器

编译：

```bash
make solve
```

运行：

```bash
./solve 文件名
```

## 说明

由于jsoncpp并不原生支持四精度浮点数，因此您的json文件会被默认以double的精度读取，若您需要高精度的初值，请将您的数值打上双引号，这样我们会读取字符串并转化为四精度浮点数。

## 支持的RK方法列表

| 名称 | 计算类型 | 步长类型 | 支持的阶数 |
| :-: | :-: | :-: | :-: |
| Classical RK | 显式 | 固定 | 4 |
| ESDIRK | 隐式 | 固定 | 4(3) | 4 |
| Gauss-Legendre | 隐式 | 固定 | 2s |
|Radau-IIA | 隐式 | 固定 | 2s-1 |
|Fehlberg | 隐式 | 自适应 | 4(5) |
|Dormand-Prince | 显式 | 自适应 | 5(4) |
|Dormand-Prince 8(7) | 显式 | 自适应 | 8(7) |
|Adaptive ESDIRK | 隐式 | 自适应 | 4 |
|Adaptive Gauss-Legendre | 隐式 | 自适应 | 2s |
|Adaptive Radau-IIA | 隐式 | 自适应 | 2s-1 |
|Embedded ESDIRK | 隐式 | 自适应 | 4(3), 5(4), 6(5) |

带括号的阶数，括号内数字表示其 embedded 的方法的阶数。

## 参考文献（部分）

1. Kennedy, Christopher & Carpenter, Mark. (2019). Diagonally Implicit Runge–Kutta Methods for Stiff ODEs. Applied Numerical Mathematics. 146. 10.1016/j.apnum.2019.07.008. 

2. P.J. Prince and J.R. Dormand. “High order embedded Runge-Kutta formulae”. In: Journal of Computational and
Applied Mathematics 7.1 (1981), pp. 67–75. ISSN: 0377-0427. DOI: https://doi.org/10.1016/0771-050X(81)
90010-3.