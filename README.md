# PBRT-Learning
I change it to a Qt version just because I want to see the images directly, and maybe it needs to be compiled by Visual Studio. You need to change the Qt version after opening the .sln file.

Since the second part continues the first part, so there is only a final project of the second part. And since the third part refactors the second part a lot, I separate them.

Besides, I found there are errors when I put some functions just in header files, so I added many sources files.

Next, I want to try to use this frame to open more complicated models like triangle-based meshes, and try to use OpenMP to speed up.



**2021.12.13：导入三角网格**

1.使用的三角网格结构和PBRT有些区别。PBRT里面三角形是与Mesh绑定的，即每个三角形结构里面都有一个Mesh指针。但是基于每个三角形也可以作为单独的形状的考虑，把三角网格结构改写成了类似 vector<Triangle> 的形式。并且在网格结构里加了一个accelerator指针，对三角网格里的三角形先做一次加速结构划分，然后把三角网格作为一个整体放入场景，再对场景做一次加速结构划分。采用这样方法的原因是发现如果不把所有三角形作为一个整体的话，把所有三角形放入场景时有一些面片会显示不出来，而且对三角网格先做加速划分获得的速度提升非常大（如果代码没写错的话）。在内存上目前还没有做考虑，应该有可以提升的部分。

2.使用了OpenMesh。如果不使用的话可以直接把trianglenew.h里的#define USE_OPENMESH注释掉。自己写的读obj文件的程序只能读最简单的，如Cat_head.obj。

3.引入了PBRT的transform。

4.上周说的使用PBRT的加速结构有些参数结果有问题，发现主要产生影响的是maxPrims这个参数，即叶子结点里面所能包含的最大Prim数，这个参数如果比较小的话结果会有问题。猜测了一下原因，比如使用kdtree的话，对空间做划分，maxPrims是每个空间结点里面能包含的Prim数。假设maxPrims=1，如果场景中有两个bounding box相交的话，那么其实无论怎么划分空间都没有办法把这两个bounding box分开分成每个node只含其中一个。

5.其实现在也不大清楚PBRT的加速结构搬得有没有问题23333，虽然加速效果很明显。