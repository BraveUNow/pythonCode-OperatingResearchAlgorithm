# pythonCode-OR_algorithm
本文档包含本人在学习运筹学过程中学习的一些算法的python代码实现，主要是线性规划&amp;网络流方面的一些经典算法，部分算法使用Gurobi solver进行求解。由于学习顺序较乱且部分算法在学习过程中没有编写代码，文档内容会较乱，之后会不断进行补充。

目前已完成的算法有:
- [单纯形法](simplex_method/simplex.py)
- [Dijkstra算法](Shortest_Path/classic_Dijkstra.py)
- [标号法求解带时间窗的最短路问题](Shortest_Path/shortestPath_labeling_TW.py)
- [列生成求解木料切割问题](Column_Generation/column_generation.py)
- Benders分解方法:[算法](Benders_Decomposition/benders.py),[个人笔记](Benders_Decomposition/Benders_Decomposition.md)