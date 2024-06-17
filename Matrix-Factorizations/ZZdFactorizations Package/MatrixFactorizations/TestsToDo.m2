Some tests to perform (in words):

1. Check periodicity of the differentials --DONE
2. Taking direct sums: can access the inclusions/projections onto the summands, should check that
the cokernel of the inclusion is equal to the other factorization --DONE
3. Basic morphisms such as the identity and 0 maps should be well-defined morphisms of factorizations
that are also commutative --DONE modulo one thing
4. Given a matrix factorization F, the endomorphisms End(F) should be a complex --DONE
5. Basic checks on all of the things we are allowed to do to ZZdfactorizations and their maps:
     for instance, we can tensor with rings/modules, apply ring maps
     there is functoriality, so all operations such as tensor/direct sum/etc should work for maps
     a basic test such as id_F ++ id_F == id_(F++F) == 1 for example (hopefully that is true)
6. A good source of ideas for tests is to look at the Complexes tests, we want similar types of things
7. Make sure that pruning/id/isdFactorization works for factorizations
   that do not necessarily have free modules
8. Basic tests with randomFactorizationMap: we should verify that any
   such map satisfies "isWellDefined", and if we do the option
   Cycle => true then the output should satisfy "isFactorizationMorphism"
   and also "isCommutative"
9. 
