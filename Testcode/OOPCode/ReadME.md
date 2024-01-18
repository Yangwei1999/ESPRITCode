
# 2023-12-19

**实验1** ：在不相关widely spaced DoA 情况下，ESPRIT GESPRIT MUSIC GMUSIC CRB 的曲线($N = 80,T =160$  SNR $=0$)
[![123](./Figure/1.jpg)  ](https://github.com/Yangwei1999/ESPRITCode/tree/DebugTest/Testcode/OOPCode/Figure/2.jpg)
此时GESPRI\ ESPRIT一致 MUSIC 和GMUSIC仍然保持一致

---

**实验2** ：在相关widely spaced DoA 情况下，ESPRIT GESPRIT MUSIC GMUSIC CRB 的曲线($N = 80,T =160$ Speration condition SNR $=0$)

![123](./Figure/2.jpg)  
此时GESPRI算法应该优于ESPRIT算法  MUSIC 和GMUSIC保持一致  

---

**实验3** ：验证N,T增大特征值和特征值的极限的谱范数距离（还有N、T较大的单独实验） 期望和方差  
3.1 先跑SNR=2的情况下的特征值极限值和理论值  
    ![Pais](./metting/3.1.png){:height="10%" width="10%"}
    theta_true = [0,pi/3];  
    coeff =10;  
    N = 40 * coeff;  
    T = 80 * coeff;  


3.2 $E[\hat{\theta} - \theta] , Var[\hat{\theta} - \theta]$

![123](./metting/3.2.png)  
方差再降低，但是偏差增高，所以可以做一个实验
看下随着N增大特征值的极限（明天做）
---
**实验4** ：ESPRIT方法与子阵列的取值数量的关系 考虑在不同的子阵列下  
貌似阵列取少后会出现 在pi/ pi 相位模糊的情况

---


<!-- **Quesitons**
* Q2 : 重根条件下无法算，但是结果表明没有问题，实际上可以忽略重根影响
* Q3 ：当不满足分离条件的时候，舍弃特征值还是不修复    
talk 2023-12-21；
1. 重根问题先不考虑
2. 不满足分离条件时候应该舍去特征值
--- -->

<!-- * 画不同$\Delta$估计器的方差 偏差
其实我们用的是TAM
尝试两种--
1. 一种是先验条件 1/ds < u < 1/ds   
2. 一种是尝试进行矫正
廖老师好，我查阅了一些资料，我发现最早的ESPRIT论文里面，也是仅考虑了阵列的偏移为ds=1（deleta），然后我找了下资料，一本书上，他是这么说的，当ds=1的时候，解是唯一的，当ds>1的时候，需要一个先验条件，就是入射角度需要在-pi/ds ~1pi/ds之间才是唯一的解，也就是随着ds增大，可获得角度的区间越小

然后我跑了两组实验，一组是两个角度，满足-pi/ds < theta_i< 1pi/ds,然后计算MSE,BIAS,var，发现确实随着ds的增大，Var MSE 是降低的，然后我看这本书也是这样

我跑了n=N/2的情况（单角度，因为多角度的话没办法区分哪个对应的真实角度是哪个），采用的是估计出来后根据2pi周期去找跟真实角度最接近的，发现他的MSE是最小的，我发现不管n怎么取，bias都是很小的（几乎无偏），然后随着n的增大，MSE和方差都减小。

我的感觉是，这个delta取1 的时候，是不需要先验条件或者预处理，可以直接用，delta>1的时候，是需要先验条件（-pi/ds < theta_i< pi/ds，或者真实角度的范围区间），这个先验条件导致随着delta的增大，估计器的方差降低 -->'


## 随着N增大偏差增大的图


## 2024-1-8 
跑了ESPRIT -  GESPRIT 偏差的相关实验   
应该把图放在Figure 文件夹中 然后按照日期分类

Q1: 貌似当theta 取某些值的时候，ESPRIT算法的理论根有错误，理论两个根的相位方向相同了。。。。。。。