# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 20:38:00 2021

@author: hyx
"""

import numpy as np


class shortPath:
    #节点个数， 距离
    def __init__(self,nodeNum,road):
        self.nodeNum = nodeNum#节点数量
        dis = np.mat(np.zeros([nodeNum,nodeNum]))
        for i in road:
            dis[i[0],i[1]] = i[2]
        self.road = road
        self.disMat = dis#距离矩阵

    def shortestPath(self):#输入距离矩阵
        n = np.size(self.disMat,axis=0)
        self.disMat[self.disMat==0] = np.Infinity#无相连路径的给一个足够大的值
        labelT = []#临时标号
        labelP = {0:0}#永久标号
        route = {0:[0]}
        for i in range(1,n):
            labelT.append(self.disMat[0,i])#第一次迭代时，每个点的T标号为出发点0到其他点的距离
            route[i] = [0]
        labelT = dict(list(zip(np.arange(1,n),labelT)))
           
        for i in range(n-1):
            shortNode = sorted(labelT.items(),key= lambda x:x[1])[0][0]#根据value排序
            labelP[shortNode] = sorted(labelT.items(),key= lambda x:x[1])[0][1]#找到最小标号的节点，并入永久标号P
            '''a = np.array(list(zip(labelT.keys(),labelT.values())))#转为array
            shortNode = a[a[:,1].argmin(),0]#最小T标号对应的node
            shortNode = int(shortNode)
            labelP[shortNode] = a[:,1].min()#添加该node为P标号'''
            del labelT[shortNode]#T标号中删除该node       
            route[shortNode].append(shortNode)
            #对于剩下T标号的所有Node，重新计算最短距离
            for j in labelT.keys():
                if labelP[shortNode]+self.disMat[shortNode,j] <= labelT[j]:
                    labelT[j] = labelP[shortNode]+self.disMat[shortNode,j]
                    route[j] = route[shortNode].copy()#到达j之前的路径改为已被放入P标号的最优路径
            self.labelP = labelP
            self.route = route
        return labelP,route

if __name__ == '__main__':
    road = [[0,1,3],[1,2,3],[0,3,4],[1,4,2],[1,5,3],[4,5,3],[3,6,3],
            [2,8,5],[5,6,1],[6,5,1],[5,8,2.5],[6,8,2],[6,7,2],[7,8,4]]
    shortPath = shortPath(9,road)
    disMatrix = shortPath.disMat
    shortPath.shortestPath()
    for i in range(len(shortPath.route)):
        print('到第%d个点的最短路径为:点%d'%(i+1,shortPath.route[i][0]+1),end='')
        for j in range(1,len(shortPath.route[i])):
            print('→点%d'%(shortPath.route[i][j]+1),end='')
        print('\n')
        
    
    
    
            
