# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 19:06:40 2021

@author: hyx
"""

from gurobipy import *

type = [4,5,7]#需求样式
timber = [9,14,16]#木料长度
cost = [5,9,10]#木料成本
demand = [30,20,40]#各个样式需求量

mp = Model()
sp = []
for i in range(3):
    sp.append(Model()) #3个子问题

zp = mp.addVars(len(type),obj=cost[0],vtype=GRB.CONTINUOUS,name='z')

cp = mp.addConstrs(quicksum(zp[j]*(timber[0]//type[j]) for j in range(len(type)) if i==j) \
                              >=demand[i] for i in range(len(demand)))
#初始全部选第一种木料，且每种方案只切一种类型的木材          
mp.setParam(GRB.Param.LogToConsole,0)                 
mp.optimize()

dualSolution = mp.getAttr(GRB.Attr.Pi,mp.getConstrs())#得到对偶解

#定义子问题
s1 = sp[0].addVars(len(type),obj=dualSolution,vtype=GRB.INTEGER)
#子问题的目标系数为RMP的影子价格，求最大，令reducted cost最小
sp[0].addConstr(quicksum(s1[i]*type[i] for i in range(len(type))) <= timber[0])#木材长度小于木料
sp[0].ModelSense=GRB.MAXIMIZE#最大化
sp[0].setParam(GRB.Param.LogToConsole,0)#不显示log
sp[0].optimize()


s2 = sp[1].addVars(len(type),obj=dualSolution,vtype=GRB.INTEGER)
sp[1].addConstr(quicksum(s2[i]*type[i] for i in range(len(type))) <= timber[1])
sp[1].ModelSense=GRB.MAXIMIZE
sp[1].setParam(GRB.Param.LogToConsole,0)    
sp[1].optimize()

s3 = sp[2].addVars(len(type),obj=dualSolution,vtype=GRB.INTEGER)
sp[2].addConstr(quicksum(s3[i]*type[i] for i in range(len(type))) <= timber[2])
sp[2].ModelSense=GRB.MAXIMIZE
sp[2].setParam(GRB.Param.LogToConsole,0)    
sp[2].optimize()

rc = [cost[0]-sp[0].objVal,cost[1]-sp[1].objVal,cost[2]-sp[2].objVal]#得到reducted cost
#reducted cost<0, 则引入此变量
while min(rc) < 0:
    k = rc.index(min(rc))#得到最小rc的子问题索引
    columnCoeff = sp[k].getAttr('X',sp[k].getVars())
    column = Column(columnCoeff,mp.getConstrs())
    #MP问题该引入变量的系数列向量为其SP的基解
    mp.addVar(obj=cost[k],vtype=GRB.CONTINUOUS,column=column)
    mp.optimize()
    
    for i in range(len(type)):
        #RMP的的影子价格改变，更改所有子问题的目标系数构成新的子问题
        s1[i].obj = cp[i].Pi
        s2[i].obj = cp[i].Pi
        s3[i].obj = cp[i].Pi
    sp[0].optimize()
    sp[1].optimize()
    sp[2].optimize()
    rc = [cost[0]-sp[0].objVal,cost[1]-sp[1].objVal,cost[2]-sp[2].objVal]#再次得到recuted cost

for v in mp.getVars():
    v.vtype = GRB.INTEGER #对原问题添加回整数约束
mp.optimize()
#输出结果
if mp.status == GRB.OPTIMAL:
    print('objive value =',mp.objVal)
    for v in mp.getVars():
            print('%s:%g'%(v.varName,v.x))
    
    
    
    
    
    
    
    
