import numpy as np

class simplex:
    def __init__(self) -> None:
        pass
    
    def normalize(self):  #将矩阵标准化
        n = np.shape(self.a)[0]#约束数量
        relaxedMat = np.identity(n, dtype='float')
        self.a = np.hstack((self.a, relaxedMat))#将单位矩阵从右侧合并到系数矩阵a
        self.c = np.append(self.c, np.zeros((1, n)))

    def reducedCost(self): #目标系数，系数矩阵，基变量索引
        n = np.shape(self.a)[0]
        cBase = []
        for i in self.base:
            cBase.append(self.c[i])
        cBase = np.array(cBase)#基变量的目标系数
        rc = self.c - np.dot(cBase, self.a)
        return rc   #返回检验数

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
        s1, s2 = np.shape(self.a)#约数个数，变量个数
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
    
    def solve(self, mat1, c, mat2=[]): #小于约束系数矩阵，目标系数，等于约束矩阵
        self.a = mat1[0]
        self.b = mat1[1]
        self.c = c
        s1, s2 = np.shape(self.a)#不等式约束个数，变量个数
        self.status = 'Optimal'
        if mat2:    #若存在等式约束
            s3 = np.size(mat2[1])#等式约束个数
            self.a = np.vstack((self.a, mat2[0]))#扩展矩阵
            self.b = np.hstack((self.b, mat2[1]))
            temp = self.c#记录原目标函数系数
            self.optimize(equalConstr=np.size(s3))
            for i in range(s1+s2, s1+s2+s3):
                if i in self.base:  #当人工变量出现在基变量内，则无可行解
                    self.status = 'Infeasible'
                    break
            if self.status == 'Optimal':    #第二阶段
                self.c = temp
                self.c = np.append(c, np.zeros((1, s1)))
                self.a = np.delete(self.a, [i for i in range(s1+s2,s1+s2+s3)],axis = 1)#删除人工变量所在列
                #print(self.a, self.b, self.c)
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
            for i in range(np.size(mat1[1])):
                if i in self.base:   #如果变量x_{i}为基变量
                    self.x['x%d'%(i+1)] = self.b[int(np.argwhere(self.base==i))]
                else:
                    self.x['x%d'%(i+1)] = 0
            self.obj = np.dot(c[:np.size(mat1[1])],np.array(list(self.x.values())))
        
        
if __name__ == '__main__':
    #a = np.array([[1, 0, 1], [1, 2, 0], [0, 1, 0]])
    a = np.array([[0.5, 0.25], [-1, -3]])
    b = np.array([4, -20]).transpose()
    c = np.array([-2, -3])
    aEq = np.array([[1, 1]])
    bEq = np.array([10])
    #b = np.array([5, 10, 4]).transpose()
    #c = np.array([2, 3, 1])
    simplex = simplex()
    simplex.solve([a, b], c, [aEq, bEq])#目标值
    print('变量值为:',simplex.x,'\n目标函数值为%g'%simplex.obj)
    
    
 


    