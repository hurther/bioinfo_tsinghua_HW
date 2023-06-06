## HMMbasic Tutorial
### #1 使用说明
该项目依赖numpy和scipy，使用前请确保已经安装这两个库。

一共实现了两个模型，分别是HMMv1和HMMv2，其中HMMv1是基于概率传递的模型，HMMv2是基于概率的对数值传递的模型。当观测序列过长时，HMMv1会出现概率下溢的情况，而HMMv2则不会。详细内容见测试HMMHMM

HMMv1和HMMv2的接口完全相同，只是内部实现不同，因此使用时可以根据需要选择使用哪个模型。**推荐使用HMMv2**。

#### #1.1 Initialization

使用时，首先需要导入HMMv1或者HMMv2，然后创建一个HMM对象。用户传入初始状态概率矩阵、状态转移概率矩阵、观测概率矩阵（或调用`init_hmm`方法进行随机初始化），可选择传入迭代次数和阈值。

```python
import numpy as np
from HMMv2 import HMMv2 
# 注意导入路径问题
# 可以将HMMv2.py和utils.py放在当前工作目录下，或者将HMMv2.py所在目录加入环境变量

# 加入环境变量方式为
# import sys
# sys.path.append("HMMv2所在目录")

states = ["box 1", "box 2", "box3"] # 状态集合
n_states = len(states)

observations = ["red", "white"] # 观测集合
n_obs = len(observations)

ip = np.array([0.2, 0.4, 0.4]) # 初始状态概率矩阵
A = np.array([[0.5, 0.2, 0.3], 
                [0.3, 0.5, 0.2],
                [0.2, 0.3, 0.5]]) # 状态转移概率矩阵
B = np.array([[0.5, 0.5],
                [0.4, 0.6],
                [0.7, 0.3]]) # 观测概率矩阵

hmm = HMMv2(ip, A, B, n_iter=100, tol=1e-3) # 创建HMM对象
hmm_random = HMMv2().init_hmm(n_state, n_obs) # 随机初始化HMM对象
```

#### #1.2 Scoring

对于输入序列，调用HMM对象的`score`方法，即可得到该序列在给定HMM模型参数下出现的概率。

```python
seq = ["red", "white", "red"]
trans_seq = [observations.index(i) for i in seq] # 将序列转换为数字序列
hmm.score(trans_seq) # -2.0385
```

#### #1.3 Decoding

对于输入序列，调用HMM对象的`decode`方法，即可得到该序列在给定HMM模型参数下的最优状态序列。

```python
box = hmm.decode(trans_seq)
print("The ball picked:", ", ".join(map(lambda x: observations[x], seen)))
print("The hidden box", ", ".join(map(lambda x: states[x], box)))
# The ball picked: red, white, red
# The hidden box: box 3, box 3, box 3
```

#### #1.4 Learning

对于输入序列，调用HMM对象的`fit`方法，即可根据输入序列学习调整初始HMM模型的参数。

```python
# 调整前
print("Initial state probability matrix:")
print(hmm_random.ip)
print("Transition probability matrix:")
print(hmm_random.a)
print("Emission probability matrix:")
print(hmm_random.b)
# 学习参数
X = np.array([[0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                [0, 0, 0, 0, 0, 1, 1, 0, 1, 1],
                [1, 1, 1, 1, 1, 0, 0, 1, 0, 0]])
hmm_random.fit(X)
# 调整后
print("Initial state probability matrix:")
print(hmm_random.ip)
print("Transition probability matrix:")
print(hmm_random.a)
print("Emission probability matrix:")
print(hmm_random.b)
```

### #2 设计说明

#### #2.1 接口说明

```python
def __init__(self, init_probability=None, transition_probability=None, emission_probability=None, n_iter=100, tol=1e-2)
```
构造函数，创建HMM对象。用户可以传入初始状态概率矩阵、状态转移概率矩阵、观测概率矩阵，可以选择传入迭代次数和阈值。其中，初始状态概率矩阵、状态转移概率矩阵、观测概率矩阵的形状分别为(n_states, )、(n_states, n_states)、(n_states, n_obs)。

```python
def init_hmm(self, nStates=None, nObservations=None)
```
随机初始化HMM对象。用户需传入状态数和观测数，随机生成初始状态概率矩阵、状态转移概率矩阵、观测概率矩阵。

```python
def forward(self, obs=None)
```
前向算法，计算观测序列的概率。用户需传入观测序列（一维数组），返回观测序列的概率。注意HMMv1中传递及输出的是概率值，HMMv2中传递和输出的是对数概率值。

```python
def backward(self, obs=None)
```
后向算法，计算观测序列的概率。用户需传入观测序列（一维数组），返回观测序列的概率。注意HMMv1中传递及输出的是概率值，HMMv2中传递和输出的是对数概率值。

```python
def score(self, obs=None)
```
计算观测序列的概率。用户需传入观测序列（一维数组），返回观测序列的概率。注意HMMv1中传递及输出的是概率值，HMMv2中传递和输出的是对数概率值。

```python
def decode(self, obs=None)
```
Viterbi算法，计算观测序列的最优状态序列。用户需传入观测序列（一维数组），返回观测序列的最优状态序列。

```python
def fit(self, obs=None)
```
Baum-Welch算法，根据观测序列学习调整初始HMM模型的参数。用户需传入观测序列矩阵(n_samples, n_observe)，返回学习后的HMM对象。

### #3 测试说明

一共有三个测试文件`test_HMMv1.py`, `test_HMMv2.py`, `test_overflow.py`，分别对应HMMv1、HMMv2、概率下溢的测试。

此处使用hmmlearn库中的hmm.MultinomialHMM模型作为参考，对HMMv1和HMMv2进行测试。测试结果如下：

```bash
cd test
# test_HMMv1.py
python test_HMMv1.py
# output:
# ip_random: [0.18205878 0.46212909 0.35581214]
# a_random: [[0.65738127 0.17132261 0.17129612]
#  [0.03807826 0.56784481 0.39407693]
#  [0.41686469 0.01211874 0.57101657]]
# b_random: [[0.79676223 0.20323777]
#  [0.4978376  0.5021624 ]
#  [0.36699967 0.63300033]]
# -------------------- hmmlearn --------------------
# The ball picked: red, white, red
# The hidden box: box3, box3, box3
# -2.038545309915233
# 0.13021800000000003
# [0.34662761 0.63698553 0.01638686]
# [[3.53968296e-01 2.69385782e-01 3.76645921e-01]
#  [1.60218661e-02 4.59487938e-01 5.24490196e-01]
#  [1.09174404e-02 4.17056326e-04 9.88665503e-01]]
# [[0.65619559 0.34380441]
#  [0.70929824 0.29070176]
#  [0.31287601 0.68712399]]
# -2.5414169427405087
# -------------------- HMMv1 --------------------
# The ball picked: red, white, red
# The hidden box box3, box3, box3
# -2.038545309915233
# 0.130218
# Converged after 42 iterations.
# [1.69215864e-22 9.99999998e-01 1.95242973e-09]
# [[6.92954516e-01 1.18437910e-01 1.88607574e-01]
#  [6.98391358e-07 9.69604537e-01 3.03947645e-02]
#  [1.26105306e-01 1.46597776e-02 8.59234917e-01]]
# [[1.15436060e-05 9.99988456e-01]
#  [2.54851226e-01 7.45148774e-01]
#  [1.35390161e-01 8.64609839e-01]]
# -3.339911027271047

# test_HMMv2.py
python test_HMMv2.py
# output:
# ip_random: [0.18205878 0.46212909 0.35581214]
# a_random: [[0.65738127 0.17132261 0.17129612]
#  [0.03807826 0.56784481 0.39407693]
#  [0.41686469 0.01211874 0.57101657]]
# b_random: [[0.79676223 0.20323777]
#  [0.4978376  0.5021624 ]
#  [0.36699967 0.63300033]]
# -------------------- hmmlearn --------------------
# The ball picked: red, white, red
# The hidden box: box3, box3, box3
# -2.038545309915233
# 0.13021800000000003
# [0.34662761 0.63698553 0.01638686]
# [[3.53968296e-01 2.69385782e-01 3.76645921e-01]
#  [1.60218661e-02 4.59487938e-01 5.24490196e-01]
#  [1.09174404e-02 4.17056326e-04 9.88665503e-01]]
# [[0.65619559 0.34380441]
#  [0.70929824 0.29070176]
#  [0.31287601 0.68712399]]
# -2.5414169427405087
# -------------------- HMMv2 --------------------
# The ball picked: red, white, red
# The hidden box box3, box3, box3
# -2.038545309915233
# 0.13021800000000003
# Converged at iteration 18.
# [0.35223989 0.63360325 0.01415686]
# [[3.48738519e-01 2.83276569e-01 3.67984912e-01]
#  [1.50640658e-02 4.54298578e-01 5.30637357e-01]
#  [7.84604780e-03 3.05296923e-04 9.91848655e-01]]
# [[0.661235   0.338765  ]
#  [0.70852795 0.29147205]
#  [0.31289763 0.68710237]]
# -2.5390284209583673

# test_overflow.py
python test_overflow.py
# output:
# --------------------hmmlearn--------------------
# -1963.708757299727
# --------------------HMMv1--------------------
# 0.0
# --------------------HMMv2--------------------
# -1963.708757299727
```

从测试结果可以看出，HMMv1相比于hmmlearn的结果有一定的差距，而HMMv2的结果与hmmlearn的结果基本一致甚至效果更好。

### #4 改进方向

1. 目前HMMv2采用python实现，效率较低，核心算法可以考虑使用C++/Cython实现，提高效率。
2. `fit`目前仅支持输入相同长度的观测序列，可以考虑支持输入不同长度的观测序列。
3. `score`目前仅支持输入单个观测序列，可以考虑支持输入多个观测序列。
4. 类设计较为简单，可以考虑增加更多的接口和将任务拆分，提高灵活性。
5. 可以考虑实现更多的HMM模型，如高斯HMM、混合HMM等。