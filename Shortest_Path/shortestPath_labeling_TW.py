import numpy as np

class shortPathTW:
    def __init__(self) -> None:
        pass
    
    #寻找未处理标记
    def untreated(self, labelQ, labelP):
        labelT = {}
        for i in labelQ:
            labelT[i] = []
            for j in labelQ[i]:
                if j not in labelP[i]:   #若j不是已处理标记，加入未处理标记
                    labelT[i].append(j)
        return labelT

    #寻找最佳标记
    def findBestT(self, labelT):
        temp = {}
        #列举所有T标号，key为标号，value为对应的node（方便查找）
        for i in labelT:
            for j in labelT[i]:
                temp[j] = i
        tempKey = list(temp.keys())
        tempKey = sorted(tempKey, key=lambda x: (x[0], x[1]))#先根据时间，再根据成本排序
        bestT = tempKey[0]#最佳标号
        bestNode = temp[bestT]
        return bestT, bestNode

    #检查新的标号是否被dominated
    def EFF(self, newState, labelQ, node):
        for i in labelQ[node]:
            if newState[0] >= i[0] and newState[1] >= i[1]:
                return True
        return False
    
    #节点，各节点时间窗，节点间行驶时间，节点间行驶
    def solve(self, node, timeConstraint, timeConsume, cost):
    #初始化永久标记Q，已处理标记P
        labelQ = {}
        labelP = {}
        for i in node:
            labelQ[i] = []
            labelP[i] = []
        labelQ[1].append((0, 0))#初始点的标号为(0, 0)，时间，费用
        shortestPath = {(0, 0): [1]}#最短路径
        while labelP != labelQ:
            labelT = self.untreated(labelQ, labelP)#未处理标号集合
            bestT, bestNode = self.findBestT(labelT)#找对最小的标号进行路径拓展
            for arc in timeConsume:
                if bestNode == arc[0]:   #找到最佳标号对应Node的后继节点
                    #若该节点到后继节点没有违背时间窗约束的话
                    if bestT[0] + timeConsume[arc] <= timeConstraint[arc[1]][1]:
                        newState = (max(timeConstraint[arc[1]][0], bestT[0] + timeConsume[arc]), bestT[1] + cost[arc])
                        #如果新的标号没有被dominated
                        if self.EFF(newState, labelQ, arc[1]) == False:
                            labelQ[arc[1]].append(newState)#加入永久标号Q
                            oldPath = shortestPath[bestT].copy()
                            #该标号对应的最短路径
                            shortestPath[newState] = []
                            shortestPath[newState].extend(oldPath)
                            shortestPath[newState].append(arc[1])
            labelP[bestNode].append(bestT)#加入已处理标号
        self.labelQ = labelQ#所有永久标号
        self.minCost = {}#各点的最小成本
        self.shortestPath = {}#各点最短路
        for i in self.labelQ:
            minCost = min(labelQ[i], key=lambda x:x[1])#最小成本标号
            self.minCost[i] = minCost[1]
            self.shortestPath[i] = shortestPath[minCost]
        

if __name__ == '__main__':
    node = [1, 2, 3, 4]
    timeConstraint = {2: [3, 10], 3: [4, 6], 4: [4, 10]}#时间窗
    timeConsume = {(1, 2): 3, (1, 4): 5, (2, 3): 2, (4, 2): 1}
    cost = {(1, 2): 2, (1, 4): -7, (2, 3): 5, (4, 2): 1}
    shortPathTW = shortPathTW()
    shortPathTW.solve(node, timeConstraint, timeConsume, cost)
    #输出结果
    for i in shortPathTW.labelQ:
        print('到达{}的最小成本为{}'.format(i, shortPathTW.minCost[i]), end=',')
        print('对应路径为{}'.format(shortPathTW.shortestPath[i]))
