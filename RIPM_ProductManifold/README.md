# README

- [RIPM_ProductManifold.m](RIPM_ProductManifold.m)

  Modified version of RIPM that works for product manifolds composed of two components, the original RIPM.m cannot handle it.

- [test_OnlyIneq_RIPM_ProductManifold.m](test_OnlyIneq_RIPM_ProductManifold.m) (P2) Only inequality constraints.

- [test_ProductManifold.m](test_ProductManifold.m) (P1) Equality constraints.

- [Notes_17th_September_2023.pdf](Notes_17th_September_2023.pdf) Records examples in P1 and P2 forms.

## modified_supporting_functions

Extended versions of the original files with the same name.

- [ProdInner.m](modified_supporting_functions\ProdInner.m)
- [ProdLincomb.m](modified_supporting_functions\ProdLincomb.m)
- [RIPM_checkupto2ndorder.m](modified_supporting_functions\RIPM_checkupto2ndorder.m)

## new_supporting_functions

Newly added files, some functions for operating on structs by fields.

- [addStructs.m](new_supporting_functions\addStructs.m)
- [divideStructs.m](new_supporting_functions\divideStructs.m)
- [findMinInStruct.m](new_supporting_functions\findMinInStruct.m)
- [innerProductStructs.m](new_supporting_functions\innerProductStructs.m)
- [maxZeroStruct.m](new_supporting_functions\maxZeroStruct.m)
- [multiplyStructWithScalar.m](new_supporting_functions\multiplyStructWithScalar.m)
- [multiplyStructs.m](new_supporting_functions\multiplyStructs.m)
- [normOfStruct.m](new_supporting_functions\normOfStruct.m)
- [subtractStructs.m](new_supporting_functions\subtractStructs.m)
