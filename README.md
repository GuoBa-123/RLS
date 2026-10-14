# RLS

**面向嵌入式的递归最小二乘（Recursive Least Squares）在线参数估计库。**

header-only，纯 C++17，无第三方依赖。矩阵维度由模板参数在编译期确定，运行期零动态内存分配，适合在 MCU 上做系统辨识、自适应滤波和时变参数跟踪。

## 特性

- **header-only** —— 只有 `rls.h` 与 `matrix.h`，`#include` 即用，无需编译链接。
- **编译期定型** —— 尺寸是模板参数，维度不匹配在编译期就会报错，运行期没有额外开销。
- **零堆分配** —— 内部仅使用 `std::array`，没有 `new`/`malloc`，可用于禁止动态内存的工程。
- **单精度** —— 全程 `float`，适配 MCU 的单精度 FPU。
- **MIMO** —— `M` 为输入维度，`N` 为输出维度，`N = 1` 即常见的 SISO 辨识。
- **遗忘因子** —— 支持跟踪时变参数。
- **数值保护** —— 分母保护、协方差矩阵对称化、`trace(P)` 上限，抑制单精度下的发散。

## 文件结构

```text
.
├── rls.h            # RLS 估计器
├── matrix.h         # 轻量定长矩阵库 
├── CMakeLists.txt   # header-only 目标 
└── README.md
```

## 快速开始

```cpp
#include "rls.h"

int main() {
    // M = 3 输入维度，N = 1 输出维度
    rls<3, 1> est(0.99f, 1000.0f);   // (遗忘因子 lambda, P0 系数 alpha)

    Mat<3, 1> x{1.0f, 2.0f, 3.0f};   // 当前输入
    Mat<1, 1> y{4.0f};               // 当前测量值

    est.update(x, y);

    const rls<3, 1>::rls_t *s = est.getRlsMessage();
    float theta0 = s->theta(0, 0);   // 估计出的参数
}
```

CMake 工程中通过 `add_subdirectory` 引入：

```cmake
add_subdirectory(RLS)
target_link_libraries(your_target PRIVATE rls::rls)
```

不使用 CMake 时，直接把仓库目录加入头文件搜索路径即可：

```sh
g++ -std=c++17 -I<path>/RLS your_source.cc
```

## API

### 构造

```cpp
rls<M, N> est;                  // theta = 0, P = alpha * I
rls<M,N> est(lambda)            //只输入遗忘因子，alpha默认0.99
rls<M, N> est(lambda, alpha);   // 指定遗忘因子与 P0 系数
```

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `lambda` | `0.99f` | 遗忘因子。越大记忆越长、越稳；越小跟踪越快、越易受噪声影响 |
| `alpha` | `1000.0f` | `P` 初值系数，`P0 = alpha * I`。越大初始收敛越快 |

运行期可用 `setLambda(float)` 调整遗忘因子。

### 迭代

```cpp
void update(const Mat<M,1> &input, const Mat<N,1> &output);  // 推荐
void calculate();                                            // 先填 state.x / state.y，再计算
```

- `update(x, y)` 写入 `state` 后立即执行一次迭代，是最常用的形式。
- `calculate()` 只执行一次迭代，输入输出需事先放在 `state.x` / `state.y` 中。适合已经持有 `rls_t`、例如从消息队列直接取数的场景。

### 读取结果

```cpp
const rls_t *getRlsMessage() const;
```

返回内部状态结构体的指针：

```cpp
typedef struct {
    Mat<N,1> y;      // 当前测量值
    Mat<M,1> x;      // 当前输入
    Mat<M,N> theta;  // 参数向量（估计结果）
} rls_t;
```

### 复位

```cpp
void reset();   // theta 清零，P 恢复为 alpha * I，增益与误差清零
```

## 算法

标准 RLS 迭代。`x` 为输入列向量，`y` 为输出列向量，`θ` 为 `M×N` 参数矩阵，`λ` 为遗忘因子：

```text
Px    = P * x
denom = λ + xᵀ * Px
e     = y - θᵀ * x        // 先验误差
k     = Px / denom        // 增益
θ     = θ + k * eᵀ
P     = (P - k * Pxᵀ) / λ
```

实现上有三点与教科书写法不同，数学上等价：

1. **复用 `Px`** —— `k * xᵀ * P` 利用 `P` 的对称性等于 `k * (Px)ᵀ`，把每步乘法量从 O(M³) 降到 O(M²)。这是嵌入式场景下的主要优化。
2. **对称化** —— 每步结束后执行 `P = (P + Pᵀ) / 2`。单精度下反复迭代会让 `P` 缓慢失去对称性，进而影响正定性，这一步用于抑制累积误差。
3. **windup 保护** —— `λ < 1` 且输入激励不足时 `P` 会按 `1/λ` 指数增长，`float` 迅速溢出，因此限制 `trace(P)` 不超过 `MAX_P_TRACE`，超限时整体按比例缩放。

分母 `λ + xᵀPx` 小于 `1e-6` 时跳过本次迭代，避免除零。

## 注意事项

- **遗忘因子与激励** —— 被辨识对象长时间缺乏充分激励时，`λ < 1` 会让协方差膨胀。这类场景建议把 `λ` 取到 `0.99` 以上，或改用带方向遗忘 / 变遗忘因子的变体。
- **`alpha` 的选取** —— `P0` 过大使初始几步对噪声过度敏感，过小则收敛慢。常用区间为 `100 ~ 1e4`，按信噪比调整。
- **`MAX_P_TRACE`** —— 按实际参数数量与 `float` 动态范围调整。
- **单精度限制** —— 本库以 `float` 为准。若被辨识系统的参数尺度相差悬殊（例如同时存在 `1e-3` 与 `1e3` 量级），建议先做归一化。
- **矩阵求逆** —— `matrix.h` 提供 `inv()`（1×1 / 2×2 / 3×3 直接展开，其余阶数用带主元的高斯-约当消元），RLS 主循环并不需要。

## 环境要求

- C++17
- 仅依赖标准库头文件

## 另：
现在用了#define Mat，如果你**有变量命名为Mat或用了别的矩阵库会导致冲突**
（不过有理由相信你只有这一个矩阵库）