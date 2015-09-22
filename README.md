# abl
Reliability in erasure-coded storage systems( Average bytes lost per usable TB)


在Greenan的HFRS的基础上修改得到abl
HFRS相关资料：
http://www.kaymgee.com/Kevin_Greenan/Software.html

由于HFRS门槛较高，代码中的一些参数也需要修改。
随着模型不同，参数的变化也不是能轻松掌握的

我根据实践经验，将HFRS原本要计算的erasure-coded storage systems的范围压缩到了reed solomon codes
删去了 XOR LDPC等系统模型

又根据至少对三块盘故障的容忍度的要求修改了代码中的默认模拟模型

详细内容请关注微信公共号  写个定理

或者我的简书： http://www.jianshu.com/users/d68483ca45f6/latest_articles

PS:
用python编写(不适用python3)   请先安装mpmath, scipy and numpy

