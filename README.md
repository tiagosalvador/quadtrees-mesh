# Quadtrees Mesh
This repo contains code for the paper [Higher-order Adaptive Finite Difference Methods for Fully Nonlinear Elliptic Equations](https://doi.org/10.1007/s10915-017-0586-5).

## Code

The script `example_meshes.m` illustrates the wide range of different meshes that can be created. It displays the following meshes.
  
<p align="center">
<img align="middle" src="http://www.math.lsa.umich.edu/~saldanha/circle.png" width="200" height="150" />
<img align="middle" src="http://www.math.lsa.umich.edu/~saldanha/ellipse.png" width="200" height="150" />
<img align="middle" src="http://www.math.lsa.umich.edu/~saldanha/diammond_stretched.png" width="200" height="150" />
<img align="middle" src="http://www.math.lsa.umich.edu/~saldanha/clover.png" width="200" height="150" />
</p>

The script `example_MA.m` contains an example on how to call the solvers for the Monge-Ampere equation with both the monotone scheme and the filtered scheme.

The script `example_ConEnv.m` contains an example on how to call the solvers for the Convex Envelope equation with both the standard mesh and and a refined mesh.

## Citation
If you use this paper or code in your scientific work, please cite as
```
@article{MeshFreeHamfeldtSalvador,
  title     = {Higher-order adaptive finite difference methods for fully nonlinear elliptic equations},
  author    = {Hamfeldt, Brittany Froese and Salvador, Tiago},
  journal   = {Journal of Scientific Computing},
  volume    = {75},
  number    = {3},
  pages     = {1282--1306},
  year      = {2018},
  publisher = {Springer}
}
```
