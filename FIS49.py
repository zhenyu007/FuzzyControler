# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 09:47:39 2020

@author: zhenyu
"""
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
import math
import copy


class Integer():#用来计算积分的类
    def __init__(self,a1,a2,b2,b1,b0,step =math.pow(10,-4)):#math.pow(10,-3)
        
        self.a2 = a2#上界
        self.a1 = a1#下界
        self.step = step#步长
        
        #函数系数
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        
        self.n = math.floor((self.a2 - self.a1)/self.step)#n个完整的直角梯形
        self.s = 0#积分结果
        
    def func(self,x):#函数
        return self.b2*math.pow(x,2) + self.b1*x + self.b0
    
    def discreNum(self,k):#第k个离散数
        return self.a1 + self.step * k
    
    
    def calculate(self): 
        i = 0
        for i in range(self.n):#遍历所有完整的小直角梯形
            self.s += 0.5 * self.step * (self.func(self.discreNum(i)) + self.func(self.discreNum(i + 1)))
            #self.s += self.step * self.func(self.discreNum(i ))
        #self.s +=  0.5*(self.a2 - self.a1 - self.step * self.n) * (self.func(self.discreNum(i + 1)) + self.func(self.a2))#加上不完整的直角梯形
        
        
        #print("被积函数为：{0}x^2 + {1}x + {2}".format(self.b2,self.b1,self.b0))
        #print("积分下限为：{0}，积分上限为：{1}".format(self.a1,self.a2))
        #print("步长为：{0}，积分结果为：{1}".format(self.step,self.s))
        #print()
        return self.s

class Lv(Enum):
    NB = 0
    NM = 1
    NS = 2
    ZO = 3
    PS = 4
    PM = 5
    PB = 6


class Rule():#规则类
    def __init__(self,Lv_E,Lv_EC,Lv_dKp,Lv_dKi,Lv_dKd):#对规则的描述
        self.Lv_E = Lv_E
        self.Lv_EC = Lv_EC
        self.Lv_dKp = Lv_dKp
        self.Lv_dKi = Lv_dKi
        self.Lv_dKd = Lv_dKd
    
    
    def calcMsv(self,e,ec,E,EC):#calculate Membership Value计算隶属度值，并取最小值
        
         E_ms = E.msSet[self.Lv_E.name]#E的某个模糊子集
         EC_ms = EC.msSet[self.Lv_EC.name]#EC的某个模糊子集
         
         #输入对E和EC的适配度
         a1 = E_ms.calc(e)
         a2 = EC_ms.calc(ec)
         
         #返回最小值
         if (a1 >= a2):
             return a2
         else:
             return a1
         
    
         
    def compare(self,x,c1,c2,maxflag):#比较一个区间内两个一次函数的大小
        
        a1 = c1[0]*x + c1[1]
        a2 = c2[0]*x + c2[1]
        if(a1>a2):
            if (maxflag == True):
                return c1
            else:
                return c2
        elif(a1 == a2):
            return -1
            
        else:
            if (maxflag == True):
                return c2
            else:
                return c1
            
            
            
            
    def whichregion(self,x,points):#给定输入，返回输入所在的区间
        if (x < points[0] or x > points[-1]):
            return -1
        elif (x == points[0]):
            return 0
        else:
            for p in points:
                if(x <= p):
                    return points.index(p)-1
     
        
        
        
    def cutTop(self,ms,a):#削顶运算,参数为一个模糊子集和一个适配度
        
        
        cutTopPoints = []#削顶后的分段点
        cutTopCoefficient = []#削顶后的函数系数
        
        
        for p  in  ms.points:#加入隶属函数的分段点
            cutTopPoints.append(p)
        for i in range (len(ms.coefficient)):#加入切割线与隶属函数的交点
            k = ms.coefficient[i][0]
            b = ms.coefficient[i][1]
            cross = (a-b)/k
            cutTopPoints.append(cross)  
        cutTopPoints = list(set(cutTopPoints))
        cutTopPoints.sort()
        
        
        for i in range(len(cutTopPoints)-1):
            region = [cutTopPoints[i],cutTopPoints[i+1]]
            c1 = [0,a]
            x = (region[0] + region[1])/2
            num = self.whichregion(x,ms.points)
            c2 = ms.coefficient[num]
            
            coeff = self.compare(x,c1,c2,False)
            
            cutTopCoefficient.append(copy.deepcopy(coeff))
            
        return cutTopPoints,cutTopCoefficient
    
    
    
    
            
    def get_outms(self,dKp,dKi,dKd):#得到输出的模糊子集
        dKp_ms = dKp.msSet[self.Lv_dKp.name]#dKp的某个模糊子集
        dKi_ms = dKi.msSet[self.Lv_dKi.name]#dKi的某个模糊子集
        dKd_ms = dKd.msSet[self.Lv_dKd.name]#dKd的某个模糊子集
        return dKp_ms,dKi_ms,dKd_ms
            
            
    
    
    
    
    
    
    def  integrate(self,p1,p2,c1,c2):#将削顶后的图形合成
        
        
        integratePoints = []   #所有规则合成后的分段点
        integrateCoefficient = [] #所有规则合成后的函数参数
        tempPoints = [] #临时的分段点
        
        
        for p in p1:
            tempPoints.append(p)
        for p in p2:
            tempPoints.append(p)
        
        tempPoints = list(set(tempPoints))
        tempPoints.sort()
        
        for p in tempPoints:#将临时分段点加入
            integratePoints.append(p)
            
        
        for i in range(len(tempPoints) - 1):
            
            x1 = tempPoints[i]#第一个输入
            x2 = tempPoints[i+1]#第二个输入
            
            num2_1 = self.whichregion(x2,p1)#第二个输入在第一个图形的哪个区间
            if num2_1 == -1:
                c2_1 = [0,0]
            else:
                c2_1 = c1[num2_1]#对应区间的函数参数
                
            num2_2 = self.whichregion(x2,p2)#第二个输入在第二个图形的哪个区间
            if num2_2 == -1:
                c2_2 = [0,0]
            else:   
                c2_2 = c2[num2_2]#对应区间的函数参数
            
            result1 = self.compare(x1,c2_1,c2_2,True)
            result2 = self.compare(x2,c2_1,c2_2,True)
            
            if(result1 == result2):#在两点之间无交叉
                if result1 == -1:#两线重合
                    integrateCoefficient.append(copy.deepcopy(c2_1))
                else:#两线平行
                    integrateCoefficient.append(copy.deepcopy(result1))
            elif result1 == -1:#端点相等
                integrateCoefficient.append(copy.deepcopy(result2))
            elif result2 == -1:#端点相等
                integrateCoefficient.append(copy.deepcopy(result1))
                
            else:#两点之间有交叉
                cross = (c2_2[1] - c2_1[1])/(c2_1[0] - c2_2[0])
                integratePoints.append(cross)
                new_result1 =  self.compare((cross+tempPoints[i])/2,c2_1,c2_2,True)
                new_result2 =  self.compare((cross+tempPoints[i+1])/2,c2_1,c2_2,True)
                integrateCoefficient.append(copy.deepcopy(new_result1))
                integrateCoefficient.append(copy.deepcopy(new_result2))
        
        integratePoints = list(set(integratePoints))
        integratePoints.sort()
        return integratePoints,integrateCoefficient
    
    def draw(self,points,coefficient,color = 'b'):
        drawpoints = []
        for i in range(len(points)-1):
            drawpoints.append([points[i],points[i+1]])
            
        for i in range (len(drawpoints)):
            x = np.linspace(drawpoints[i][0],drawpoints[i][1],100)
            y = coefficient[i][0]*x + coefficient[i][1]
            plt.plot(x,y,color,linewidth=2, markersize=2)
        
            
            
        #plt.legend(loc = "best")
        plt.show()
        
    
                
                
                
    

   
        
class Var():#变量类
    def __init__(self,universe):
        self.universe = universe#变量的论域
        self.name = ["NB","NM","NS","ZO","PS","PM","PB"]#模糊子集的名称
        self.msSet ={}#模糊子集的集合
        for i in range (7):#创建七个模糊子集,并存入字典msSet中
            self.msSet[self.name[i]] = Ms(-universe + (i - 1) * universe/3,-universe 
                      + i * universe/3,-universe + (i + 1) * universe/3,universe)
            
           
        

class Ms():#模糊子集类
    def __init__(self,a,b,c,universe):
        
        #三角形三个参数
        self.a = a
        self.b = b
        self.c = c
        
        
        #论域范围
        self.universe = universe
        
        
        #函数分段点,及其函数解析式
        if(self.a < -self.universe):#NB
            self.points = [self.b,self.c]#分段点
            self.drawpoints =[[self.b,self.c]]#画图用的分段点
            self.coefficient =[[1/(b-c),-c/(b-c)]]#函数系数
        elif(self.c > self.universe):#PB
            self.points = [self.a,self.b]
            self.drawpoints =[[self.a,self.b]]
            self.coefficient =[[1/(b-a),-a/(b-a)]]
        else:#除去两端的一般情况
            self.points = [self.a,self.b,self.c]
            self.drawpoints =[[self.a,self.b],[self.b,self.c]]
            self.coefficient =[[1/(b-a),-a/(b-a)],[1/(b-c),-c/(b-c)]]
        
        
    def calc(self,x):#计算隶属度
        isflag = -1#标志输入变量属于哪个区间
        for i in range(len(self.drawpoints)):
            if (x>=self.drawpoints[i][0] and x<=self.drawpoints[i][1]):#如果在第i个区间
                isflag = i
                break
        
        if (isflag == -1):#不属于任何区间，即隶属度为零
            return 0
        else :
            return self.coefficient[isflag][0]*x + self.coefficient[isflag][1]
            
        
                
        
    def draw(self):#绘出隶属度图像
        for i in range (len(self.drawpoints)):
            x = np.linspace(self.drawpoints[i][0],self.drawpoints[i][1],100)
            y = self.coefficient[i][0]*x + self.coefficient[i][1]
            plt.plot(x,y,'bo',linewidth=2, markersize=2)
            
        plt.legend(loc = "best")
        
   
    
    
        
    
    
    


    
def fis(e,ec,E,EC,rule,dKp,dKi,dKd):#推理得到一条规则的三个输出
    
    a = rule.calcMsv(e,ec,E,EC)#最小适配度
    
    
    dKp_ms,dKi_ms,dKd_ms = rule.get_outms(dKp,dKi,dKd)#输出的模糊子集
    
    
    p1,c1 = rule.cutTop(dKp_ms,a)#dKp的削顶后的隶属函数
    #rule1.draw(p1,c1)
    
    p2,c2 = rule.cutTop(dKi_ms,a)#dKi的削顶后的隶属函数
    #rule1.draw(p2,c2)
    
    p3,c3 = rule.cutTop(dKd_ms,a)#dKd的削顶后的隶属函数
    #rule1.draw(p3,c3)
    
    return p1,c1,p2,c2,p3,c3


    
def defuzzy(points,coefficient):#解模糊
    myregions = []
    myup = copy.deepcopy(coefficient)
    mydown = copy.deepcopy(coefficient)
    length = len(points)-1
    for m in range(length):
        myregions.append([points[m],points[m+1]])#积分上下界
    for m in range(length):
        myup[m].append(0)#分子
        mydown[m].insert(0,0)#分母 
    #print("myup is:{}".format(myup))
    #print("mydown is:{}".format(mydown))
    #print("myregions is:{}".format(myregions))
    
    upvalue = 0#分子值
    downvalue = 0#分母值
    
    for i in range(len(points)-1):
        try:
            
            integer1 = Integer(myregions[i][0], myregions[i][1], myup[i][0],myup[i][1], myup[i][2])
            upvalue += integer1.calculate()
            
            integer2= Integer(myregions[i][0], myregions[i][1], mydown[i][0],mydown[i][1], mydown[i][2])
            downvalue += integer2.calculate()
        except:
            print("积分异常")
            
        
    return upvalue/downvalue#合成后图形重心的横坐标
    

def  fuzzyControl(e,ec,Kx,all_virs,rulesNum = 49):
    #得到单个规则的输出
    for i in range(1,rulesNum + 1):
        all_virs["output"+str(i)] = fis(e,ec,E,EC,all_virs["rule"+str(i)],dKp,dKi,dKd)
        
        
        
        
    stylelist = ["dKp","dKi","dKd"]
    Kx = Kx#0,1,2分别代表PID
    #所有规则的输出合成
    temp_p = output1[2*Kx]
    temp_c = output1[2*Kx + 1]
    
    for i in range(2,rulesNum + 1):
        temp_p,temp_c = rule1.integrate(temp_p,all_virs["output"+str(i)][2*Kx],temp_c,all_virs["output"+str(i)][2*Kx + 1])

    rule1.draw(temp_p,temp_c)  
    
    #print("temp_p is:{}".format(temp_p))
    #print("temp_c is:{}".format(temp_c))
    u = defuzzy(temp_p,temp_c)  #解模糊
    print("{0}={1}".format(stylelist[Kx],u))    

    
    
if __name__ =="__main__":
        
    #添加输入输出变量
    E = Var(6)
    EC = Var(6)
    dKp = Var(3)
    dKi = Var(3)
    dKd = Var(3)
    
    
    #建立模糊规则
    rule1 = Rule(Lv.NB,Lv.NB,Lv.PB,Lv.NB,Lv.PS)
    rule2 = Rule(Lv.NB,Lv.NM,Lv.PB,Lv.NB,Lv.NS)
    rule3 = Rule(Lv.NB,Lv.NS,Lv.PM,Lv.NM,Lv.NB)
    rule4 = Rule(Lv.NB,Lv.ZO,Lv.PM,Lv.NM,Lv.NB)
    rule5 = Rule(Lv.NB,Lv.PS,Lv.PS,Lv.NS,Lv.NB)
    rule6 = Rule(Lv.NB,Lv.PM,Lv.ZO,Lv.ZO,Lv.NM)
    rule7 = Rule(Lv.NB,Lv.PB,Lv.ZO,Lv.ZO,Lv.PS)
    
    
    rule8 = Rule(Lv.NM,Lv.NB,Lv.PB,Lv.NB,Lv.PS)
    rule9 = Rule(Lv.NM,Lv.NM,Lv.PB,Lv.NB,Lv.NS)
    rule10 = Rule(Lv.NM,Lv.NS,Lv.PM,Lv.NM,Lv.NB)
    rule11 = Rule(Lv.NM,Lv.ZO,Lv.PS,Lv.NS,Lv.NM)
    rule12 = Rule(Lv.NM,Lv.PS,Lv.PS,Lv.NS,Lv.NM)
    rule13 = Rule(Lv.NM,Lv.PM,Lv.ZO,Lv.ZO,Lv.NS)
    rule14 = Rule(Lv.NM,Lv.PB,Lv.NS,Lv.PS,Lv.ZO)
    
    
    rule15 = Rule(Lv.NS,Lv.NB,Lv.PM,Lv.NB,Lv.ZO)
    rule16 = Rule(Lv.NS,Lv.NM,Lv.PM,Lv.NM,Lv.NS)
    rule17 = Rule(Lv.NS,Lv.NS,Lv.PM,Lv.NS,Lv.NM)
    rule18 = Rule(Lv.NS,Lv.ZO,Lv.PS,Lv.NS,Lv.NM)
    rule19 = Rule(Lv.NS,Lv.PS,Lv.ZO,Lv.ZO,Lv.NS)
    rule20 = Rule(Lv.NS,Lv.PM,Lv.NS,Lv.PM,Lv.NS)
    rule21 = Rule(Lv.NS,Lv.PB,Lv.NS,Lv.PM,Lv.ZO)
    
    
    rule22 = Rule(Lv.ZO,Lv.NB,Lv.PM,Lv.NM,Lv.ZO)
    rule23 = Rule(Lv.ZO,Lv.NM,Lv.PM,Lv.NM,Lv.NS)
    rule24 = Rule(Lv.ZO,Lv.NS,Lv.PS,Lv.NS,Lv.NS)
    rule25 = Rule(Lv.ZO,Lv.ZO,Lv.ZO,Lv.ZO,Lv.NS)
    rule26 = Rule(Lv.ZO,Lv.PS,Lv.NS,Lv.PS,Lv.NS)
    rule27 = Rule(Lv.ZO,Lv.PM,Lv.NM,Lv.PM,Lv.NS)
    rule28 = Rule(Lv.ZO,Lv.PB,Lv.NM,Lv.PM,Lv.ZO)
    
    
    rule29 = Rule(Lv.PS,Lv.NB,Lv.PS,Lv.NM,Lv.ZO)
    rule30 = Rule(Lv.PS,Lv.NM,Lv.PS,Lv.NS,Lv.ZO)
    rule31 = Rule(Lv.PS,Lv.NS,Lv.ZO,Lv.ZO,Lv.ZO)
    rule32 = Rule(Lv.PS,Lv.ZO,Lv.NS,Lv.PS,Lv.ZO)
    rule33 = Rule(Lv.PS,Lv.PS,Lv.NS,Lv.PS,Lv.ZO)
    rule34 = Rule(Lv.PS,Lv.PM,Lv.NM,Lv.PM,Lv.ZO)
    rule35 = Rule(Lv.PS,Lv.PB,Lv.NM,Lv.PB,Lv.ZO)
    
    
    rule36 = Rule(Lv.PM,Lv.NB,Lv.PS,Lv.ZO,Lv.PB)
    rule37 = Rule(Lv.PM,Lv.NM,Lv.ZO,Lv.ZO,Lv.PM)
    rule38 = Rule(Lv.PM,Lv.NS,Lv.NS,Lv.PS,Lv.PM)
    rule39 = Rule(Lv.PM,Lv.ZO,Lv.NM,Lv.PS,Lv.PM)
    rule40 = Rule(Lv.PM,Lv.PS,Lv.NM,Lv.PM,Lv.PS)
    rule41 = Rule(Lv.PM,Lv.PM,Lv.NM,Lv.PB,Lv.PS)
    rule42 = Rule(Lv.PM,Lv.PB,Lv.NB,Lv.PB,Lv.PB)
    
    
    rule43 = Rule(Lv.PB,Lv.NB,Lv.ZO,Lv.ZO,Lv.PB)
    rule44 = Rule(Lv.PB,Lv.NM,Lv.ZO,Lv.ZO,Lv.PM)
    rule45 = Rule(Lv.PB,Lv.NS,Lv.NM,Lv.PS,Lv.PM)
    rule46 = Rule(Lv.PB,Lv.ZO,Lv.NM,Lv.PM,Lv.PM)
    rule47 = Rule(Lv.PB,Lv.PS,Lv.NM,Lv.PM,Lv.PS)
    rule48 = Rule(Lv.PB,Lv.PM,Lv.NB,Lv.PB,Lv.PS)
    rule49 = Rule(Lv.PB,Lv.PB,Lv.NB,Lv.PB,Lv.PB)
    
    
    #给定输入
    input1 = -5#E[-6,6]
    input2 = -3#EC[-6,6] 
    
    all_virs=locals()#所有变量的字典
    
    
    
    #进行模糊推理、解模糊，并打印出结果，参数为e,ec,Kx
    fuzzyControl(input1,input2,0,all_virs)#dKp
    fuzzyControl(input1,input2,1,all_virs)#dKi
    fuzzyControl(input1,input2,2,all_virs)#dKd
    
    
   
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    