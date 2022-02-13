import numpy as np

class simplex:
    def __init__(self) -> None:
        pass
    
    def normalize(self, c, a):  #将矩阵标准化
        n = np.shape(a)[0]#约束数量
        relaxedMat = np.identity(n, dtype='float')
        newMatrix = np.hstack((a, relaxedMat))#将单位矩阵从右侧合并到系数矩阵a
        newObjective = np.append(c, np.zeros((1, n)))
        return newObjective, newMatrix

    def reducedCost(self, c,a, base): #目标系数，系数矩阵，基变量索引
        n = np.shape(a)[0]
        cBase = []
        for i in base:
            cBase.append(c[i])
        cBase = np.array(cBase)#基变量的目标系数
        rc = c - np.dot(cBase, a)
        return rc   #返回检验数

    def baseChange(self, a, b, rc, base):   #更新基变量
        swapIn = np.argmax(rc)  #换入基变量的索引
        n = np.shape(a)[0]#约束个数
        theta = [b[i]/ a[i, swapIn] if a[i, swapIn]>0 else np.Infinity for i in range(n)]#计算入基变量
        theta = [np.Infinity if theta[i] <= 0 else theta[i] for i in range(n)]
        swapOut = np.argmin(theta)#得到换出变量
        base[swapOut] = swapIn  #基变量索引更换
        for i in range(n):  #更新系数矩阵
            if i == swapOut:    #若为轴心所在行
                b[i] = b[i]/a[swapOut, swapIn]
                a[i, :] = a[i, :]/a[swapOut, swapIn]
            else:
                b[i] = b[i] - b[swapOut]*a[i, swapIn]/a[swapOut, swapIn]
                a[i, :] = a[i, :] - a[swapOut, :]*a[i, swapIn]/a[swapOut, swapIn]
        return a, b, base
    
    def optimize(self, a, b, c):    #系数矩阵，约束系数， 目标系数
        s1, s2 = np.shape(a)#约数个数，变量个数
        base = np.arange(s2, s2+s1, dtype='int')#基变量索引
        c, a =  self.normalize(c, a)
        rc = self.reducedCost(c, a, base)
        while np.max(rc) > 0:   #当存在正的检验数
            a, b, base = self.baseChange(a, b, rc, base)
            rc = self.reducedCost(c, a, base)
        self.a = a
        self.b = b
        self.c = c
        self.base = base
        #输出结果
        x = {}
        for i in range(s1):
            if i in base:   #如果变量x_{i}为基变量
                x['x%d'%(i+1)] = b[int(np.argwhere(base==i))]
            else:
                x['x%d'%(i+1)] = 0
        print(x)
        obj = np.dot(c[:s1],np.array(list(x.values())))
        self.x = x
        self.obj = obj
    

if __name__ == '__main__':
    a = np.array([[1, 0, 1], [1, 2, 0], [0, 1, 0]])
    b = np.array([5, 10, 4]).transpose()
    c = np.array([2, 3, 1])
    simplex = simplex()
    simplex.optimize(a, b, c)#目标值
    print('变量值为:',simplex.x,'\n目标函数值为%g'%simplex.obj)
    
    
 


    