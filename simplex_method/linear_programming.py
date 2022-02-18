from matplotlib import rc_file
import numpy as np
import copy

class simplex:
    def __init__(self, mat1, c, mat2=[], Inte=[]):
        self.bestObjective = -np.Infinity
        self.a = mat1[0]
        self.b = mat1[1]
        self.c = c
        s1, s2 = np.shape(self.a)#不等式约束个数，变量个数
        self.constrNum = s1
        self.unEqConstr = s1#不等式约束数量
        self.eqConstr = 0
        self.variableNum = s2
        if mat2:    #若存在等式约束
            s3 = np.size(mat2[1])#等式约束个数
            self.eqConstr = s3
            self.constrNum += s3#总约束为不等式约束+等式约束
            self.a = np.vstack((self.a, mat2[0]))#扩展矩阵
            self.b = np.hstack((self.b, mat2[1]))
        if Inte:
            self.Int = Inte#定义整数变量
        self.cOrigin = self.c.copy()
            
    def normalize(self):  #将矩阵标准化
        relaxedMat = np.identity(self.constrNum, dtype='float')
        self.a = np.hstack((self.a, relaxedMat))#将单位矩阵从右侧合并到系数矩阵a
        self.c = np.append(self.c, np.zeros((1, self.constrNum)))

    def reducedCost(self): #目标系数，系数矩阵，基变量索引
        cBase = []
        for i in self.base:
            cBase.append(self.c[i])
        cBase = np.array(cBase)#基变量的目标系数
        rc = self.c - np.dot(cBase, self.a)
        return rc   #返回检验数
    
    def delSlackVariable(self):
        self.a = self.a[:, :self.variableNum]#删去松弛变量列

    def baseChange(self, rc):   #更新基变量
        swapIn = np.argmax(rc)  #换入基变量的索引
        n = np.shape(self.a)[0]#约束个数
        theta = [self.b[i]/ self.a[i, swapIn] if self.a[i, swapIn]!=0 else np.Infinity for i in range(n)]#计算入基变量
        theta = [np.Infinity if theta[i] <= 0 else theta[i] for i in range(n)]
        if (np.array(theta)==np.Infinity).all():
            self.status == 'Unbounded'
            return
        swapOut = np.argmin(theta)#得到换出变量
        self.base[swapOut] = swapIn  #基变量索引更换
        for i in range(n):  #更新系数矩阵
            if i == swapOut:    #若为轴心所在行
                self.b[i] = self.b[i]/self.a[swapOut, swapIn]
                self.a[i, :] = self.a[i, :]/self.a[swapOut, swapIn]
            else:
                self.b[i] = self.b[i] - self.b[swapOut]*self.a[i, swapIn]/self.a[swapOut, swapIn]
                self.a[i, :] = self.a[i, :] - self.a[swapOut, :]*self.a[i, swapIn]/self.a[swapOut, swapIn]
    
    def optimize(self,equalConstr=0, secondStage=False): 
        s1, s2 = np.shape(self.a)#约束个数，变量个数
        if secondStage == False:    #两阶段法第二阶段
            self.base = np.arange(s2, s2+s1, dtype='int')#基变量索引
            self.normalize()
        if equalConstr > 0: #若存在人工变量，两阶段法第一阶段
            #第一阶段的人工变量目标系数为1，其余为0
            self.c = [-1 if i >= s1+s2-equalConstr else 0 for i in range(s1+s2)]
        rc = self.reducedCost()
        while np.max(rc) > 0:   #当存在正的检验数
            self.baseChange(rc)
            if self.status == 'Unbounded':#若为无界解，跳出循环
                break
            rc = self.reducedCost()
    
    def solve(self): #小于约束系数矩阵，目标系数，等于约束矩阵
        aTemp = self.a.copy()
        bTemp = self.b.copy()
        self.status = 'Optimal'
        if self.eqConstr > 0:    #若存在等式约束
            temp = self.c#记录原目标函数系数
            self.optimize(equalConstr=np.size(self.eqConstr))
            for i in range(self.unEqConstr+self.variableNum, self.constrNum+self.variableNum):
                if i in self.base:  #当人工变量出现在基变量内，则无可行解
                    self.status = 'Infeasible'
                    break
            if self.status == 'Optimal':    #第二阶段
                self.c = temp
                self.c = np.append(self.c, np.zeros((1, self.unEqConstr)))
                self.a = np.delete(self.a, [i for i in range(self.unEqConstr+self.variableNum, self.constrNum+self.variableNum)],axis = 1)#删除人工变量所在列
                self.optimize(secondStage=True)
        else:
            self.optimize()
        #输出结果
        if self.status == 'Infeasible':
            print('the model is infeasible')
        elif self.status == 'Unbounded':
            print('the model is unbounded')
        else:
            self.x = {}
            for i in range(self.variableNum):
                if i in self.base:   #如果变量x_{i}为基变量
                    self.x[i] = self.b[int(np.argwhere(self.base==i))]
                else:
                    self.x[i] = 0
            for i in range(self.constrNum):
                if np.dot((aTemp[i, :self.variableNum]),np.array(list(self.x.values()))) > bTemp[i]:
                    self.status = 'Infeasible'
                    print('the model is infeasible')
                    break
            self.obj = np.dot(self.c[:self.variableNum],np.array(list(self.x.values())))
     
    def isInteger(self): #是否所有整数约束都得到满足，如果不是则添加约束
        for i in self.x:
            if i in self.Int:    #如果i为整数变量
                if self.x[i]%1 != 0:
                    return i, self.x[i]
        return True
    
    def addLeftConstr(self, ind, num):    #插入变量索引及对应值
        newConstr = np.array([0.0 if i != ind else 1.0 for i in range(np.shape(self.a)[1])])#新的行
        self.a = np.vstack((newConstr, self.a))
        #为新的约束添加松弛变量
        self.b = np.append(num, self.b)
        self.unEqConstr += 1
        self.constrNum += 1
        
    def addRightConstr(self, ind, num):
        newConstr = np.array([0.0 if i != ind else -1.0 for i in range(np.shape(self.a)[1])])#新的行
        self.a = np.vstack((newConstr, self.a))
        #为新的约束添加松弛变量
        self.b = np.append(-num, self.b)
        self.unEqConstr += 1
        self.constrNum += 1
    
    def branchBound(self):
        self.layer = 1
        self.matrix = {}
        self.matrix[self.layer] = [self.a.copy(), self.b.copy()]
        self.matrix1 = {}
        self.branchAndBound()
           
    def branchAndBound(self):
        self.c = self.cOrigin.copy()
        self.solve()
        if self.status == 'Optimal':    #当有最优解，开始分支
            print('解为:',self.x)
            if self.obj > self.bestObjective:
                if self.isInteger() == True:    #若整数约束得到满足
                    self.bestObjective = self.obj
                    self.bestX = self.x
                    self.layer -= 1
                else:   #若变量不为整数
                    ind, num = self.isInteger()#非整数变量索引和对应的值
                    self.matrix1[self.layer] = [ind, num]
                    self.a = self.matrix[self.layer][0].copy()
                    self.b = self.matrix[self.layer][1].copy()
                    self.addLeftConstr(self.matrix1[self.layer][0], np.floor(self.matrix1[self.layer][1]))#添加小于约束
                    print('添加约束:x%d<=%d,'%(ind, np.floor(self.matrix1[self.layer][1])),end='')
                    self.layer += 1
                    self.matrix[self.layer] = [self.a.copy(), self.b.copy()]#更新下层的ab
                    #print(self.a, self.b)
                    self.branchAndBound()#递归
                    #恢复原矩阵
                    self.a = self.matrix[self.layer][0].copy()
                    self.b = self.matrix[self.layer][1].copy()
                    self.unEqConstr -= 1
                    self.constrNum -= 1
                    self.addRightConstr(self.matrix1[self.layer][0], np.ceil(self.matrix1[self.layer][1]))
                    print('添加约束:x%d>=%d,'%(ind, np.ceil(self.matrix1[self.layer][1])),end='')
                    self.layer += 1
                    self.matrix[self.layer] = [self.a.copy(), self.b.copy()]#更新下层的ab
                    self.branchAndBound()
                    self.unEqConstr -= 1
                    self.constrNum -= 1
                    self.layer -= 1
            else:
                self.layer -= 1
        else:
            self.layer -= 1
            

        
if __name__ == '__main__':
    a = np.array([[0.5, 0.25], [-1, -3]])
    b = np.array([4, -20]).transpose()
    c = np.array([-2, -3])
    aEq = np.array([[1, 1]])
    bEq = np.array([10])
    #a = np.array([[1, 0, 1], [1, 2, 0], [0, 1, 0]])
    #b = np.array([5, 10, 4]).transpose()
    #c = np.array([2, 3, 1])
    lp1 = simplex([a, b], c, mat2=[aEq, bEq])
    lp1.solve()#目标值
    print('变量值为:',lp1.x,'\n目标函数值为%g'%lp1.obj)
    
    a1 = np.array([[2.0, 1.0], [5.0, 7.0]])
    b1 = np.array([9.0, 35.0]).transpose()
    c1 = np.array([6.0, 5.0])
    lp2 = simplex([a1,b1], c1, Inte=[0, 1])
    lp2.branchBound()
    print('变量值为:',lp2.x,'\n目标函数值为%g'%lp2.obj)
    
    
 


    