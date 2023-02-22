A basic implementation via *[VCGlib](http://vcg.isti.cnr.it/vcglib)* of the developable surface modeling approach for triangle meshes proposed by Oded Stein, Eitan Grinspun, and Keenan Crane in the paper *[Developability of Triangle Meshes](https://doi.org/10.1145/3197517.3201303)*.

**This tool is also available in [*Meshlab*](http://www.meshlab.net/) as a [filter plugin](https://github.com/cnr-isti-vclab/meshlab/pull/1308)**, allowing for a more user-friendly usage within an environment packed with many tools for processing, editing and visualizing mesh properties.

![samples](https://user-images.githubusercontent.com/65720263/220656110-688194ec-e269-4e8b-a020-1e207ef1ac2d.png)

## Building
```
mkdir build
cd build
cmake ..
make
./develop
```

## Usage
```
exec [path to mesh] [optimizer] [max no. fun evals] [stop gradient threshold] [step size] {tau}
```

`optimizer` can be:

* `f`: gradient method opt with fixed step size
* `b`: gradient method opt with backtracking line search (Armijo condition)

Optimization stops when the max number of function evaluations is reached, or the squared 2-norm gradient goes below the given threshold.

When using backtracking opt, step size is progressively scaled down by factors of `tau` according to the Armijo condition.

Example usage:

```
./build/develop ./assets/bunny.obj b 400 0 0.0001 0.8
```

## About discrete developability

In a smooth and continuous setting, a surface is said to be *developable* if it can be isometrically mapped to a plane while also being a twice differentiable immersion. The first condition, called *flattenability*, requires a zero Gaussian curvature and implies the surface can be obtained from a plane without introducing any distortion; this property alone however admits surfaces that may not be easily manufactured, like a piece of paper that needs to be crumpled in a very specific way. To avoid this, it is also required the flattenable surface be a $C^2$ immersion, inducing the presence of straight ruling lines passing through points on the surface.

When dealing with triangle meshes, developability is commonly acknowledged just as *local flattenability*, which is met when the vertices have zero angle defect and thus zero Gaussian curvature over a small neighborhood; however, as in the smooth case, this condition alone may admit behaviors that are not desirable in a manufacturing scenario.

The authors of the paper suggest a new definition of developability for manifold triangle meshes that captures both the property of flattenability and presence of straight ruling lines: if all the vertex stars are flat or form an hinge, than the mesh is *discrete developable*; in both cases, the local Gaussian curvature will be zero, and any non-flat region will have parallel ruling lines running along the hinges.

Based on this definition, the authors derive an approach based on variational principles that alters the geometry of manifold triangle meshes in order to make them more developable: the approach simply amounts to a gradient descent over an energy that quantifies the developability of a mesh, given by how flat or hinge-like each of its vertex stars is; the result is then a mesh that is similar to the initial, but instead comprised of one or more developable pieces held toghether by highly regular seam curves, i.e. path of edges which vertex stars did not form an hinge or a flat spot.
