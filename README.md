# SEDRepair
SEDRepair is the source code repository for the paper:

**Sector Error-Oriented Durability-Aware Fast Repair in Erasure-Coded Cloud Storage Systems**  
Yifei Xiao, Shijie Zhou, Linpeng Zhong, UESTC, China

SEDRepair is based on Zerasure(<https://github.com/zhoutl1106/zerasure>).

# Usage
## Software reference
- Jerasure 1.2A: 
  http://web.eecs.utk.edu/~plank/plank/www/software.html  
  Source code of this library is partially ported into the repo for convinient compiling purpose. However, neither the authors nor Texas A&M University own this source code. All rights are reserved to Dr. Plank, University of Tennessee, Knoxville.

## Compile
Zerasue use **qmake** to organize the project. No *QT* module, however, is used in the source code.
A **qmake** *.pro* file is included in the repo. Compilation could be simply:
~~~~
qmake
make
~~~~
in the root directory of the source code. An executable file *zerasure* will be generated.
A *makefile* with default compilation flags is also provided named **mfile**, you can compile using
~~~
make -f mfile
~~~

## Pre-optimized Cauchy Matrix
Several pre-optimzed X,Y arrays to define *Cauchy matrix* are provided in **PreOpt/ge_100_03_06_01_1000_weighted_s13.txt**.  
Each row corresponds to one specific (k,m,w) parameters obtained by *genetic algorithm, weighted cost function, strategy-(1,3)* with

- initial population = 100
- select rate = 0.3
- crossover rate = 0.6
- mutation rate = 0.1
- maximum population = 1000

Within one row, the first 3 elements are (k,m,w), then the k+m other numbers for the X,Y array (see the definition of Cauchy matrix). An example of using these numbers will be shown in next section.
This will also be the output format of Simulated Annealing and Genetic Algorithm.

## Usage
### Start  
~~~
single 0 8
~~~
### Init(n,k)
Input K(the amount of data nodes),M(the amount of parity nodes),W(how many bits you use to store every single value, e.g., W=4, means 2^4=16, values from 0 ~ 15)

~~~
6 3 4
~~~
will read K,M,W,X,Y from stdin (one line in txt file), and compute the simulation result: the total repair time.

### Change Value of Factor

If you want to change the factor, like the size of packet or time of T1, you should modify the source code.



## Code structure
- Algorithm/
    - SEDRepair: implement of SEDRepair
    - ERT:reconstruction-only, which is the conventional method based on RS codes, ERT only use reconstruction operations
    - ER: ER preserves recovered data in the parity nodes in BN, which means ER is support for data durability 
- Jerasure-1.2A/
    - Partial source code ported from Jerasure-1.2A, minor modification for compile purpose.
- utils: Some utility functions such as memory management and timing.
