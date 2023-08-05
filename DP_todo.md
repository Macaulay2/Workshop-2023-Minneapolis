# Recent changes to SimplicialModule type:
## I made an option to specify whether you want the simplicialModule command to compute the degeneracy maps. The reason for this is as follows: the normalization actually does not use the information of the degeneracy maps, so computing these is unnecessary for many purposes of computing derived functors. Currently the option is set to NOT compute degeneracy maps by default, but you can add the option Degeneracy => true to compute the degeneracy maps. 

# Derived Functors and General Compatibility
## Need to write a symmetric/divided power functor that works for arbitrary matrices
## The command simplicialModuleMap converts a morphism of complexes into a SimplicialModuleMap.
## I updated the simplicialModule command to assume that the topDegree of the outputted simplicial module is equal to the length of the input complex (but probably we should add a warning that the user can specify the top degree)
## I now have 3 types of nonlinear functors implemented: Exterior powers, Schur functors (this includes symmetric powers), and tensor products. Note: I am aware that exterior powers are an example of a Schur functor, but the exterior power code runs significantly faster than the Schur functor code.


# Normalization
## I aded a command that actually has the proper normalization.
## Normalization now applies to morphisms of simplicial modules and outputs a morphism of the respective normalizations.

# Construction Optimization
## Caching A matrix
## Rewriting so only column needed is pulled by facemap builders

# Implement basic functions/binary ops
[ ] - == comparison





# General changes from 8/5 push:
## Simplicial module maps have essentially all the same functionality that complex maps have: we can add, subtract, multiply, take direct sums, tensor products, etc.
## Many basic comparisons operators such as == are implemented
## Data from the direct summands of the input data is cached when applying simplicialModule to a direct sum of objects
## There is a forgetComplex and forgetDegeneracy command if you want to forget certain pieces of data from a simplicial module (in case you want to run a computation that might take too long if this extra data is taken into account, you don't have to entirely redefine the object to get rid of it)
## When taking direct sums, more information is cached and can be recovered (particularly useful when wanting the induced map upon restricting to a particular summand)
## Can now take kernels/cokernels/images/coimages of morphisms of simplicial modules
## Can now prune simplicial modules and morphisms of simplicial modules (very useful when combined with taking kernels/etc as above)
## Can use "inducedMap" to try to get a naively induced map between simplicial modules. For instance, if phi is a morphism of simplicial modules, inducedMap(source phi,ker phi) will return the inclusion Ker phi --> source phi. 
## simplicialTensor has been updated to cache the indices of the direct summands, so morphisms between tensor products and the differentials/degeneracies can now be restricted to summands of the tensor product
## the normalize command has been altered slightly in a few ways: it now no longer prunes the output, the user should use "prune" on the output if they want it pruned (I started getting some seriously annoying bugs when the terms of the simplicial module were not just free modules; the previous data pruned the individual maps and then naively glued them back together, but this does not work if the ambient data is not free). Also: when normalizing something that has the key "complex" cached, it doesn't bother performing the normalization, it instead just outputs the value S.complex (where S is the simplicial module). If you don't want it to do that, you can use the option CheckComplex => false OR use the forgetComplex command. Likewise, similar functionality when it sees that the simplicial module is a direct sum, it instead normalizes each component then outputs the direct sum of the normalizations (this can be avoided by using the option CheckSum => false)


