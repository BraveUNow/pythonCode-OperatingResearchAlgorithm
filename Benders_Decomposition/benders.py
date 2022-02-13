from gurobipy import *
import numpy as np

mp = Model()#主问题
sp = Model()#子问题

yNum = list(np.arange(1,6))
alphaNum = list(np.arange(1,9))

y = mp.addVars(yNum,obj=7,vtype=GRB.BINARY,name='y')
z = mp.addVar(obj=1,vtype=GRB.CONTINUOUS,name='z')

alpha = sp.addVars(alphaNum,vtype=GRB.CONTINUOUS,name='alpha')
for i in alphaNum:
    if i <= 3:
        alpha[i].lb = -GRB.INFINITY#等于约束的对偶变量无约束
    else:
        alpha[i].lb = -GRB.INFINITY#最小化问题小于约束对偶变量小于0
        alpha[i].ub = 0

c1 = sp.addConstr(alpha[1]+alpha[4]<=1)
c2 = sp.addConstr(alpha[2]+alpha[5]<=1)
c3 = sp.addConstr(alpha[3]+alpha[6]<=1)
c4 = sp.addConstr(alpha[1]+alpha[3]+alpha[7]<=1)
c5 = sp.addConstr(alpha[1]+alpha[2]+alpha[8]<=1)
sp.ModelSense = GRB.MAXIMIZE

sp.setParam(GRB.Param.InfUnbdInfo,1)#获取无界解信息

while 1:
    mp.optimize()
    sp.setObjective(8*alpha[1]+3*alpha[2]+5*alpha[3]+8*y[1].x*alpha[4]+ \
                    3*y[2].x*alpha[5]+5*y[3].x*alpha[6]+5*y[4].x*alpha[7]+3*y[5].x*alpha[8])
        #得到(b-By),设定subproblem的目标函数
    sp.optimize()
    if sp.status == GRB.OPTIMAL:
        if sp.objVal == z.x:
            break
        else:
            mp.addConstr(8*alpha[1].x+3*alpha[2].x+5*alpha[3].x+8*y[1]*alpha[4].x+ \
                         3*y[2]*alpha[5].x+5*y[3]*alpha[6].x+5*y[4]*alpha[7].x+3*y[5]*alpha[8].x<=z)
                #add optimality cut
    elif sp.status == GRB.UNBOUNDED:
        ray = sp.UnbdRay#obtain extreme ray
        mp.addConstr(8*ray[0]+3*ray[1] + 5*ray[2] + \
                     8*ray[3]*y[1] + 3*ray[4]*y[2] + 5*ray[5]*y[3] + \
                     5*ray[6]*y[4] + 3*ray[7]*y[5]<= 0)
        #add feasibility cut
    else:
        print (SP_Dual.status)

#输出变量值
x = sp.getAttr(GRB.Attr.Pi,sp.getConstrs())#x为sp对偶变量
for i in range(5):
    print('x[%d]:%g'%(i+1,x[i]))
for v in mp.getVars():
    print(v.varName,':',v.x)