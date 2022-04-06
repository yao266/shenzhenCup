#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random

plt.rcParams["font.sans-serif"]=["SimHei"] #设置字体
plt.rcParams["axes.unicode_minus"]=False #该语句解决图像中的“-”负号的乱码问题




#求两点间距
def distanceOfTwoPoints(x1,y1):
    sums = 0 
    for i in range(len(x1)):
        sums += (x1[i] -y1[i])**2
    return np.sqrt(sums)

def nextCoorOfDog(pre,tsh,rudis,V,shNext):
    #求解犬在羊一个运动时间段tsh后运动的两个坐标,并判断这两个坐标与羊的距离，选择小的一个返回。
    seg = (V*tsh*180)/(np.pi*rudis)
    temp = pre[0]**2+pre[1]**2-rudis**2+2*rudis**2*(np.cos(seg*np.pi/180))
    b = 4*pre[0]*temp
    a = 4*(pre[0]**2+pre[1]**2)
    c = temp**2-4*pre[1]**2*rudis**2
    delta = b**2-4*a*c
    if delta < 0:
        return pre
    if delta == 0 :
        x = -b/(2*a)
        y = np.sqrt(rudis**2-x**2)
        return [x,y]
    if delta > 0 :
        x1 = (-b-np.sqrt(delta))/(2*a)
        y1 = np.sqrt(rudis**2-x1**2)
        x2 = (-b+np.sqrt(delta))/(2*a)
        y2 = np.sqrt(rudis**2-x2**2)
        if distanceOfTwoPoints([x1,y1],shNext) > distanceOfTwoPoints([x2,y2],shNext):
            return [x2,y2]
        else:
            return [x1,y1]

def correspondingOfDog(pre,rudis,shNext):
    #返回羊运动到下一个点对应的犬的位置点
    k = pre[1]/pre[0]
    a = k+1
    b = 0
    c = (k+1)*rudis**2
    delta = b**2-4*a*c
    if delta < 0 :
        return pre
    if delta == 0:
        return [-(b/2/a),-k*(b/2/a)]
    if delta > 0:
        x1 = (-b-np.sqrt(delta))/(2*a)
        y1 = k*x1
        x2 = (-b+np.sqrt(delta))/(2*a)
        y2 = k*x2
        if distanceOfTwoPoints([x1,y1])


#求解羊下一次运动的坐标
def nextCoorOfSheep(pre,angel,preRudis,nextRudis):
    #求解羊下一次运动位置点，传入参数：羊上一次坐标，运动方向（随机值），运动到下一个点所在圆的半径
    k = np.cos(np.pi*angel/180)
    temp = pre[0]**2+pre[1]**2+4*k*preRudis*nextRudis-nextRudis**2
    a = 4*(pre[0]**2+pre[1]**2)
    b = -4*temp*pre[0]
    c = temp**2-4*pre[1]**2*nextRudis**2
    delta = b**2-4*a*c
    if delta < 0:
        return pre
    if delta == 0 :
        x = -(b/2/a)
        y = np.sqrt(nextRudis**2-x**2)
        return [x,y]
    if delta > 0 :
        x = max((-b+np.sqrt(delta))/2/a,(b+np.sqrt(delta))/2/a)
        y = np.sqrt(nextRudis**2-x**2)
        return [x,y]

def mainModel(coordinateDataOfSheep,coordinateDataOfDog,R,N,r,ang,V,v,num,x0,sh0,timeOfSheep,timeOfDog):
#绘制同心圆与犬羊的初始位置
    theta = np.linspace(0, 2*np.pi, 200)
    figure, axes = plt.subplots(1)
    for i in range(1,N+1):
        a = i*np.cos(theta)*(R/N)
        b = i*np.sin(theta)*(R/N)
        axes.plot(a, b,c='black',alpha=0.05)
        axes.set_aspect(1)
    axes.scatter(0,0,s=3,c='r')
    axes.scatter(x0[0],x0[1],c='orange',s=3,label='Dog')
    axes.scatter(sh0[0],sh0[1],c='blue',s=3,label='Sheep')
    plt.legend(loc = 'lower left',prop={'family':'SimHei','size':'x-small'})

    #循环处理
    space = np.arange(r,R+(R/N), R/N)
    # print(space)
    for j in range(1,len(space)):
        ang  = random.randrange(0,61)
        sh1 = nextCoorOfSheep(coordinateDataOfSheep[-1],ang,space[j-1],space[j])
        # print(sh1)
        while sh1 == coordinateDataOfSheep[-1]:
            ang  = random.randrange(0,61)
            sh1 = nextCoorOfSheep(coordinateDataOfSheep[-1],ang,space[j-1],space[j])
        tsh = distanceOfTwoPoints(sh1,sh0)/v
        xtemp = correspondingOfDog(coordinateDataOfDog[-1],R,sh1)
        x1 = nextCoorOfDog(coordinateDataOfDog[-1],tsh,R,V,sh1)
        print(sh1)
        # print(x1)
        # if x1 == coordinateDataOfDog[-1]:
        #     x1 = nextCoorOfDog(coordinateDataOfDog[-1],tsh,R,V,sh1)
        timeOfSheep.append(tsh)
        timeOfDog.append(distanceOfTwoPoints(x1,coordinateDataOfDog[-1])/V)
        coordinateDataOfDog.append(x1)
        coordinateDataOfSheep.append(sh1)
        axes.scatter(coordinateDataOfDog[-1][0],coordinateDataOfDog[-1][1],c='orange',s=3,label='Dog')
        axes.scatter(coordinateDataOfSheep[-1][0],coordinateDataOfSheep[-1][1],c='blue',s=3,label='Sheep')
        # plt.show()
    plt.title('第'+str(num+1)+'批犬羊模型位置')
    plt.savefig('./data/('+str(num)+').png')
    print(coordinateDataOfSheep)
    # print(coordinateDataOfDog)
    # plt.show()
    #数据导出为csv（t,d,x,y) (t,sh,x1,y1)
    outdata = []
    for i in range(len(coordinateDataOfDog)):
        outdata.append([i,'d',"("+str(coordinateDataOfDog[i][0])+","+str(coordinateDataOfDog[i][1])+")"])
    pd.DataFrame(outdata).to_csv('./data/dog_coordinate('+str(num)+').csv')
    outdata = []
    for i in range(len(coordinateDataOfSheep)):
        outdata.append([i,'sh',"("+str(coordinateDataOfDog[i][0])+","+str(coordinateDataOfDog[i][1])+")"])
    pd.DataFrame(outdata).to_csv('./data/sheep_coordinate('+str(num)+').csv')
if __name__ =='__main__':
    V = 2 #犬的速度单位m/s
    v = 1 #羊的速度单位m/s
    vMax = 2#狼与羊的速度比最大值
    # angel = 30 #羊与犬初始位置顺时针方向所成夹角
    R = 20 #圆的半径
    N = 20 #构建N个同心圆，N趋于无穷，由此得到可微的概念
    for j in range(0,100):
        ang  = random.randrange(0,61) #羊运动到下一个点与上一点与x轴夹角
        tempx = random.randrange(-R,R)
        x0 = [tempx,np.sqrt(R**2-tempx**2)] #犬的初始位置
        r = random.randrange(0,(N-1)*R/N,R/N)
        sh0 =[] 
        if int (r) != 0:
            tempx = random.randrange(-r,r,R/N)
            sh0 = [tempx,np.sqrt(r**2-tempx**2)] #羊的初始位置
        else:
            sh0 = [0,0]
        coordinateDataOfDog = [] #犬的坐标信息
        coordinateDataOfSheep = [] #样的坐标信息
        timeOfDog = [] #犬每次运动花费时间
        timeOfSheep = [] #羊每次运动花费的时间
        coordinateDataOfDog.append(x0)
        coordinateDataOfSheep.append(sh0)
        mainModel(coordinateDataOfSheep,coordinateDataOfDog,R,N,r,ang,V,v,j,x0,sh0,timeOfDog,timeOfSheep)
# %%
