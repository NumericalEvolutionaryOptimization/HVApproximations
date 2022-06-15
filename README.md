## SMSEMOA-GP
SMSEMOA-GP is an EMOA based on mathematical models that approximate the hypervolume (HV) value of the Pareto set. The modeling strategy used to generate the models is [Genetic Programming](http://www.genetic-programming.com/) (GP), which produces models expressed as symbolic mathematical expressions. The [GP system](https://github.com/domsearson/gptips-2-0) was configured to obtain efficient and simple models, avoiding unnecessarily large or complex structures. The models use statistical descriptors as inputs, and through an efficient implementation are extremely fast and efficient. The main advantage of the resulting expressions is their low computational complexity, which significantly reduces computational times. This work represents the first time that a machine learning approach (GP) was used to generate models that approximate the HV indicator.

# Motivation

Hypervolume (HV) has become one of the most popular indicators to assess the quality of Pareto front approximations. However, the best algorithm for computing these values has a computational complexity of O(N^(k/3)polylog(N)) for N candidate solutions and k objectives. In this study, we propose a regression-based approach to learn new mathematical expressions to approximate the HV value and improve at the same time their computational efficiency.

# Structure
The content of this repository can be divided into three parts
- source code 
    - 3 objetives   
    - 4 objetives
    - 5 objetives
- data base 
    - 3 objetives   
    - 4 objetives
    - 5 objetives
- models
    - 3 objetives   
    - 4 objetives
    - 5 objetives
# Requirements
- The source codes which can be used to reproduce the experimental results in the paper are provided in the source code folder. The source codes use [PlatEMO](https://github.com/BIMK/PlatEMO/) v2.0 or higher.
- Databases can be used to train new models. In this work we use [GPTIPS2](https://github.com/domsearson/gptips-2-0) to evolve the models.
- The models can be used as a guide for the SMSEMOA-GP algorithm. It is necessary to train models for different number of objectives. We are working to overcome this limitation. 
- The models can be used as a guide for the SMSEMOA-GP algorithm. These models were trained with samples of fronts N = 300, for this population size they obtain their best results. It is necessary to train models for different number of objectives. We are working to overcome this limitation.
- The models can be used in DTLZ and WFG benchmark test suite.
# References
-  C. Sandoval, O. Cuate, L.C. González et al., Towards fast approximations for the hypervolume indicator for multi-objective optimization problems
by Genetic Programming, Applied Soft Computing (2022) 109103, https://doi.org/10.1016/j.asoc.2022.109103.
# Version
Current version: v0.1
# License
GNU GPLv3, see [LICENSE]()