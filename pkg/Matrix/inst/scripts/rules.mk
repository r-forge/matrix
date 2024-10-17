## Generated by ./rules.sh :
Csparse.o: Csparse.c Mdefines.h \
  Msymbols.h utils.h cs-etc.h \
  SuiteSparse/CXSparse/Include/cs.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h cholmod-etc.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h Csparse.h
attrib.o: attrib.c Mdefines.h \
  Msymbols.h utils.h attrib.h
bind.o: bind.c Mdefines.h \
  Msymbols.h utils.h M5.h coerce.h \
  bind.h
cholmod-api.o: cholmod-api.c Mdefines.h \
  Msymbols.h utils.h cholmod-api.h \
  cholmod-etc.h SuiteSparse/SuiteSparse_config/SuiteSparse_config.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h
cholmod-etc.o: cholmod-etc.c Mdefines.h \
  Msymbols.h utils.h idz.h \
  cholmod-etc.h SuiteSparse/SuiteSparse_config/SuiteSparse_config.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h
coerce.o: coerce.c Mdefines.h \
  Msymbols.h utils.h M5.h idz.h \
  coerce.h
cs-etc.o: cs-etc.c Mdefines.h \
  Msymbols.h utils.h cs-etc.h \
  SuiteSparse/CXSparse/Include/cs.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h
dense.o: dense.c Mdefines.h \
  Msymbols.h utils.h M5.h idz.h \
  dense.h
determinant.o: determinant.c Mdefines.h \
  Msymbols.h utils.h \
  cholmod-etc.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h determinant.h
expm.o: expm.c Lapack-etc.h \
  Mdefines.h \
  Msymbols.h utils.h expm.h
factor.o: factor.c Lapack-etc.h \
  cs-etc.h \
  SuiteSparse/CXSparse/Include/cs.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h cholmod-etc.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h Mdefines.h \
  Msymbols.h utils.h M5.h factor.h
idz.o: idz.c Mdefines.h \
  Msymbols.h utils.h M5.h idz.h
init.o: init.c Mdefines.h \
  Msymbols.h utils.h Csparse.h \
  attrib.h bind.h cholmod-api.h cholmod-etc.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h coerce.h dense.h determinant.h \
  expm.h factor.h kappa.h matmult.h objects.h perm.h solve.h sparse.h \
  subassign.h subscript.h utils-R.h validity.h \
  
kappa.o: kappa.c Lapack-etc.h \
  Mdefines.h \
  Msymbols.h utils.h kappa.h
matmult.o: matmult.c Lapack-etc.h \
  cholmod-etc.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h Mdefines.h \
  Msymbols.h utils.h M5.h idz.h coerce.h dense.h sparse.h matmult.h
objects.o: objects.c Mdefines.h \
  Msymbols.h utils.h objects.h
perm.o: perm.c Mdefines.h \
  Msymbols.h utils.h perm.h
solve.o: solve.c Lapack-etc.h \
  cs-etc.h \
  SuiteSparse/CXSparse/Include/cs.h \
  SuiteSparse/SuiteSparse_config/SuiteSparse_config.h cholmod-etc.h \
  SuiteSparse/CHOLMOD/Include/cholmod.h Mdefines.h \
  Msymbols.h utils.h M5.h idz.h solve.h
sparse.o: sparse.c Mdefines.h \
  Msymbols.h utils.h M5.h idz.h \
  sparse.h
subassign.o: subassign.c Mdefines.h \
  Msymbols.h utils.h subassign.h \
  t_subassign.c
subscript.o: subscript.c Mdefines.h \
  Msymbols.h utils.h M5.h idz.h \
  subscript.h
utils-R.o: utils-R.c Mdefines.h \
  Msymbols.h utils.h M5.h version.h \
  utils-R.h t_rle.c
utils.o: utils.c Mdefines.h \
  Msymbols.h utils.h
validity.o: validity.c Mdefines.h \
  Msymbols.h utils.h validity.h
