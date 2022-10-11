# SCENIC
SCENIC是一种同时重建基因调控网络并从单细胞RNA-seq数据中鉴定stable cell states的工具。基于共表达和DNA模基序（motif）分析推断基因调控网络 ，然后在每个细胞中分析网络活性以鉴定细胞状态。

SCENIC的分析步骤如下

第一步：GENIC3/GRNBoost，通过共表达分析找出某转录因子（TF）可能的靶基因

第二步：RcisTarget，通过分析转录因子锌指结构与靶基因的motif基序的结合情况，筛选可能的靶基因（不直接进行第二步的原因是，一个TF在不同种类的细胞中可能调控不同的基因）得到regulon

第三步：AUCell对转录因子与对应靶基因构成的基因集进行活性打分


