# Optimal Neighborhood Multiple Kernel Clustering with Adaptive Local Kernels (ON-ALK)

Matalb implementation for IEEE TKDE paper:
- Jiyuan Liu, Xinwang Liu, Jian Xiong, Qing Liao, Sihang Zhou, Siwei Wang and Yuexiang Yang:  
[Optimal Neighborhood Multiple Kernel Clustering with Adaptive Local Kernels](https://liujiyuan13.github.io/pubs/Optimal_Neighborhood_Multiple_Kernel_Clustering_with_Adaptive_Local_Kernels.pdf), IEEE Transactions on Knowledge and Data Engineering (TKDE), 2020, doi: 10.1109/TKDE.2020.3014104.

## Introduction
**Abstract**

Multiple kernel clustering (MKC) algorithm aims to group data into different categories by optimally integrating information from a group of pre-specified kernels. Though demonstrating superiorities in various applications, we observe that existing MKC algorithms usually do not sufficiently consider the local density around individual data samples and excessively limit the representation capacity of the learned optimal kernel, leading to unsatisfying performance. In this paper, we propose an algorithm, called optimal neighborhood MKC with adaptive local kernels (ON-ALK), to address the two issues. In specific, we construct adaptive local kernels to sufficiently consider the local density around individual data samples, where different numbers of neighbors are discriminatingly selected on each sample. Further, the proposed ON-ALK algorithm boosts the representation of the learned optimal kernel via relaxing it into the neighborhood area of weighted combination of the pre-specified kernels. To solve the resultant optimization problem, a three-step iterative algorithm is designed and theoretically proven to be convergent. After that, we also study the generalization bound of the proposed algorithm. Extensive experiments have been conducted to evaluate the clustering performance. As indicated, the algorithm significantly outperforms state-of-the-art methods in recent literatures on six challenging benchmark datasets, verifying its advantages and effectiveness.

## Citation

If you find our code useful, please cite:

	@article{journals/tkde/LiuONALK,
	  author    = {Jiyuan Liu and
	               Xinwang Liu and
	               Jian Xiong and
	               Qing Liao and
	               Sihang Zhou and 
	               Siwei Wang and
	               Yuexiang Yang},
	  title     = {Optimal Neighborhood Multiple Kernel Clustering with Adaptive Local Kernels},
	  journal   = {IEEE Transactions on Knowledge and Data Engineering (TKDE)},
	  volume    = {},
	  pages     = {},
	  year      = {2020},
	  url       = {https://doi.org/10.1109/TKDE.2020.3014104},
	  doi       = {10.1109/TKDE.2020.3014104}
	}

## More
- For more related researches, please visit my homepage: [https://liujiyuan13.github.io](https://liujiyuan13.github.io).
- For data and discussion, please message me: liujiyuan13@nudt.edu.cn

