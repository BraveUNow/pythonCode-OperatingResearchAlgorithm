

# Benders Decomposition

## 1.算法推导

​	对于一个线性规划问题：
$$
\text{Minimize}\quad 
c^Tx+f^Ty\\
\text{s.t. }\quad Ax+By=b\\
x\ge0,y\in Y\subseteq \mathbb{R}^q
$$
​	当y有一个比较复杂的约束（如为整数），而x为的约束较为简单（线性约束）时，模型求解computational cost会难以忍受，如果令y固定，求解x，可以缩小问题的规模，因此考虑将原模型分解成2个较小的子模型，当固定y时，相当于求一个只有变量x的子模型，令q(y)为：
$$
q(y)=\quad
\text{Minimize}\quad c^Tx\\
\text{s.t.}\quad Ax=b-By\\
x\ge0
$$
​	则原问题可以表示为：
$$
\text{Minimize}\quad 
f^Ty+q(y)\\
\text{s.t. }\quad y\in Y
$$
​	考虑到y取不同的值，q(y)的可行域发生变化，可以将子模型q(y)转化成对偶形式：
$$
\text{Maximize}\quad \alpha^T (b-By)\\
\text{s.t. }\quad A^T\alpha\le c\\
\alpha\text{ is unrestricted}
$$
​	将该模型作为子问题(subproblem), 注意到不论y的值怎么发生变化，子问题的约束不变，也就是说子问题的可行域不随输入y的变化而变化，y的不同只会影响子问题的目标值。子问题的可行域是一个**多面体(polyhedron)**，该多面体有$I(\alpha_p^1,\alpha_p^2,\cdots,\alpha_p^I)$个**极点(extreme points)**,$J(\alpha_r^1,\alpha_r^2,\cdots,\alpha_r^J)$个**极方向(extreme rays)**，对于多面体内任何一点都可以由极点和极方向的组合表示$\forall \alpha=\{\sum\limits_i^I\lambda_i\alpha_p^i+\sum\limits_j^J\mu_j\alpha_r^j\mid\sum\limits_i^I\lambda_i=1,\lambda\ge0,\mu\ge0\}$。

 - 若可行域非空且无界，最优解一定可以由可行域的某个极方向表示。

 - 若可行域非空且有界，最优解一定在可行域的某个极点上。

 - 若可行域为空集，其对偶问题无界，则原问题也无界，模型目标值可以取到无穷小，失去求解意义。

  因此，原问题可以转化成：
$$
\text{Minimize}\quad 
f^Ty+max(\alpha^T (b-By))\\
\text{s.t. }\quad A^T\alpha\le c\\
y\in Y
$$
​	等价于:
$$
\text{(MP)=Minimize}\quad 
f^Ty+z\\
\text{s.t. }\quad (\alpha_r^j)^T(b-By)\le0\quad j=1,2,\cdots,J\\
(\alpha_p^i)^T(b-By)\le z\quad i =1,2,\cdots,I\\ \\ \\
\text{(SP)=Maximize}\quad \alpha^T (b-By)\\
\text{s.t. }\quad A^T\alpha\le c\\
\alpha\text{ is unrestricted}
$$
​	将该模型作为主问题(master problem)，由于子问题的可行域可以有很多个极点和极方向，全部列举出来并不现实（也并不需要），因此我们可以先将主问题的约束全部松弛掉，形成限制主问题(restrict master problem)，并在每次迭代过程中逐步添加约束，即在满足最优性的条件下追求原问题的可行性，当限制主问题的z值与子问题的目标值相等时，主问题的所有约束都得到了满足

- 对于给定的y值，任意极点乘以(b-By)都小于子问题的目标值=z, $(\alpha_p^i)^T(b-By)\le z$
- 对于给定的y值，没有极方向与y代表的向量小于90°，即（$\alpha_r^j)^T(b-By)\le0$)，若小于90°，子问题的目标值可以无限大（$\mid\alpha_r^j\mid$可以取任意值),子问题无界，其对偶问题无可行解，原问题也无可行解

​	首先，求解RMP，得到对应的$(y^*,z^*)$。对于给定$y*$，若SP得到无界解，说明对于该无界解对应的可行解是一个极方向，且此时对于该极方向j，与向量(b-By)的夹角小于90°，$(\alpha_r^j)^T(b-By)\ge0$，因此将这个极方向添加进RMP（add Benders feasibility cuts:$\alpha_r^T(b-By)\le 0$, $\alpha_r$ is the extreme ray corresponding to current solution)；若SP得到最优解，且最优解的目标值大于$z^*$，则说明存在某个极点$\alpha_p^i$,代入目标函数得到的值大于$z^*$，则违背了MP的约束$(\alpha_p^i)^T(b-By)\le z$,因此将这个极点添加进RMP（add Benders optimality cuts:$\alpha_p^T(b-By)\le z$, $\alpha_p$ is the extreme point corresponding to current solution),这样求解RMP得到的(b-By)与该极点向量相乘一定小于求解得到的$z^*$；当最优解的目标值等于$z^*$时，说明对于所有极点，与向量(b-By)的内积都小于$z^*$，符合MP的约束，也就是说RMP得到的最优解是MP的可行解，且是MP的最优解。（松弛问题的可行域大于原问题的可行域，其最优解的目标值是原问题的下界，如果松弛问题的最优解可行，一定是原问题的最优解），这样MP也就得到了最优的解。

<img src="feasible domain.png" style="zoom:30%;" />

​	设求解RMP得到y的最优解在多面体极点B上，求解SP，得到最优点在极点E上，目标值为z，则对于$\alpha$可行域任意极点$\alpha_p^i$，$(\alpha_p^i)^T(b-By)\le z$，此时添加新的极点约束进入RMP，相当于缩小y的可行域，且缩小后的y可行域一定包含点B（因为B是符合所有极点的约束的），则$y_{new}$的最优一定也在B，SP的最优一定还在F，对结果没有影响，可以不用添加该cut。

<img src="feasible domain1.png" style="zoom:30%;" />

##  2.伪代码
\begin{align}
\text{\textbf {Classic Benders Decomposition}}\\
\rule{20.3cm}{0.05em}\\
\text{Indentify the primal MP and SP}\\
\text{repeat:}\\
\quad \text{Solve MP, get z}\\
\quad \text{change SP's objective funciton}\\
\quad \text{Solve SP}\\
\quad \text{if SP.STATUS=OPTIMAL:}\\
\quad\quad \text{if SP.obj=z:}\\
\quad\quad\quad \text{stop}\\
\quad\quad \text{else:}\\
\quad\quad\quad \text{add Benders optimality cuts:}\alpha_p^T(b-By)\le z\\
\quad \text{else:}\\
\quad\quad\quad \text{add Benders feasibility cuts:}\alpha_r^T(b-By)\le 0\\
\end{align}


## 3.算法实例
本实例来自《Wiley Encyclopedia of Operations Research and
Management Science》的一个章节，原PDF见[Benders Decompositon](http://hacivat.ie.boun.edu.tr/~taskin/pdf/taskin_benders.pdf)
​	对MIP问题:.
$$
min\quad x_1+x_2+x_3+x_4+x_5+7y_1+7y_2+7y_3+7y_4+7y_5\\
\text{s.t. }\quad 
x_1+x_4+x_5=8\\
x_2+x_5=3\\
x_3+x_4=5\\
x_1\le 8y_1\\
x_2\le 3y_2\\
x_3\le 5y_3\\
x_4\le 5y_4\\
x_5\le 3y_5\\
x_1,x_2,x_3,x_4,x_5\ge0,y_1,y_2,y_3,y_4,y_5\in \{0,1\}
$$

**step1.** 识别MP和SP, 记$\overline{(b-By)}=(8,3,5,8y_1,3y_2,5y_3,5y_4,3y_5)$:
$$
\text{(RMP)=Minimize}\quad 
7y_1+7y_2+7y_3+7y_4+7y_5+z\\
\text{s.t. }\quad z\ge0,y_1,y_2,y_3,y_4,y_5\in \{0,1\}\\ \\ \\
\text{(SP)=Maximize}\quad 8\alpha_1+3\alpha_2+5\alpha_3+8y_1\alpha_4+3y_2\alpha_5+5y_3\alpha_6+5y_4\alpha_7+3y_5\alpha_8\\
\text{s.t. }\quad \alpha_1+\alpha_4\le1\\
\alpha_2+\alpha_5\le1\\
\alpha_3+\alpha_6\le1\\
\alpha_1+\alpha_3+\alpha_7\le1\\
\alpha_1+\alpha_2+\alpha_8\le1\\
\alpha_1,\alpha_2,\alpha_3\text{ is unrestricted, }\alpha_4,\alpha_5,\alpha_6,\alpha_7,\alpha_8\le0
$$

**step2.** 求解RMP, 得最优解为$y^*=(0,0,0,0,0),z^*=0$, 将$y^*$代入SP, 发现SP有无穷解, 且我们得到该解对应的极射线$(1,0,0,-1,0,0,-1,-1)$, 之后我们可以添加feasibility cut: $8-8y_1-5y_4-3y_5\le0$,等价于$8y_1+5y_4+3y_5\ge8$

**step3.** 更新RMP,得最优解为$y^*=(1,0,0,0,0),z^*=0$, SP仍然无界, 对应的极射线为$(0,1,1,0,-1,-1,-1,-1)$, 继续添加feasibility cut: $3y_2+5y_3+5y_4+3y_5\ge8$

**step4.** 更新RMP,得最优解为$y^*=(0,0,0,1,1),z^*=0$, SP有最优解$8>z=0$, 对应的极点为$(1,0,0,0,0,0,0,0)$, 添加optimality cut: $8\le z$

**step5.** 更新RMP, 得最优解为 $y^*=(0,0,0,1,1),z^*=8$, 当前的$y^*$与上一次迭代相同, 因此SP的最优值一定为8, 等于当前的$z^*$, 算法结束. 得到SP的对偶变量, 即原问题x的取值. 且最后的RMP为:
$$
\text{(RMP)=Minimize}\quad 
7y_1+7y_2+7y_3+7y_4+7y_5+z\\
\text{s.t. }\quad 
8y_1+5y_4+3y_5\ge8\\
3y_2+5y_3+5y_4+3y_5\ge8\\
z\ge8\\
z\ge0,y_1,y_2,y_3,y_4,y_5\in \{0,1\}
$$
​	最终结果为$x=(0,0,0,5,3),y=(0,0,0,1,1)$, objective=22

Gurobi代码见[code of the instance](benders.py)



