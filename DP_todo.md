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
