# abl


Reliability in erasure-coded storage systems( Average bytes lost per usable TB)

USAGE:

关注微信公共号 写个定理

存储一栏中的可靠性计算有用法介绍


在Greenan的HFRS的基础上修改得到abl

HFRS相关资料：

http://www.kaymgee.com/Kevin_Greenan/Software.html

由于HFRS门槛较高，代码中的一些参数也需要修改。

随着模型不同，参数的变化也不是能轻松掌握的

我根据实践经验，将HFRS原本要计算的erasure-coded storage systems的范围压缩到了reed solomon codes

删去了 XOR LDPC等的入口

又根据至少对三块盘故障的容忍度的要求修改了代码中的默认模拟模型

详细内容请关注微信公共号： 写个定理

或者我的简书： http://www.jianshu.com/users/d68483ca45f6/latest_articles

PS:

用python编写(不适用python3)   请先安装mpmath, scipy , numpy


仅适用于非并行修复的场景

关于并行修复的好处和负面影响如何平衡的问题虽然重要，但不是这个程序要解决的

在上面公共号中提到的文章中有一定的篇幅提到了并行修复的问题。

并行修复最直观的好处是大大缩短了修复过程

最直接的问题是文件的scope变大，坏单位物理单位的情况下丢失数据的可能性也随之上升

隐性的问题是维护元数据的开销增长的很快，最后修复的瓶颈一方面在网络另一方面在元数据上

元数据的压力可能会导致诸多问题，甚至是直接丢失数据

————————————————————————————

个人意见，在可靠性够看的情况下（根据计算不用并行修复系统已经足够可靠了）

尝试去RS ---> LRC  才是比较正确的选择

正如facebook azure他们所做的一样

安全 可靠 而不单单追求性能

可靠的信号传输是纠删码的使命

——————————————————————————————

总之，策略这种事，要实事求是



