sources = \
	Check/cholmod_check.c Check/cholmod_l_check.c \
	Check/cholmod_read.c  Check/cholmod_l_read.c \
	Check/cholmod_write.c Check/cholmod_l_write.c \
	Cholesky/cholmod_amd.c          Cholesky/cholmod_l_amd.c \
	Cholesky/cholmod_analyze.c      Cholesky/cholmod_l_analyze.c \
	Cholesky/cholmod_colamd.c       Cholesky/cholmod_l_colamd.c \
	Cholesky/cholmod_etree.c        Cholesky/cholmod_l_etree.c \
	Cholesky/cholmod_factorize.c    Cholesky/cholmod_l_factorize.c \
	Cholesky/cholmod_postorder.c    Cholesky/cholmod_l_postorder.c \
	Cholesky/cholmod_rcond.c        Cholesky/cholmod_l_rcond.c \
	Cholesky/cholmod_resymbol.c     Cholesky/cholmod_l_resymbol.c \
	Cholesky/cholmod_rowcolcounts.c Cholesky/cholmod_l_rowcolcounts.c \
	Cholesky/cholmod_rowfac.c       Cholesky/cholmod_l_rowfac.c \
	Cholesky/cholmod_solve.c        Cholesky/cholmod_l_solve.c \
	Cholesky/cholmod_spsolve.c      Cholesky/cholmod_l_spsolve.c \
	Utility/cholmod_aat.c           Utility/cholmod_l_aat.c \
	Utility/cholmod_add.c           Utility/cholmod_l_add.c \
	Utility/cholmod_band.c          Utility/cholmod_l_band.c \
	Utility/cholmod_change_factor.c Utility/cholmod_l_change_factor.c \
	Utility/cholmod_common.c        Utility/cholmod_l_common.c \
	Utility/cholmod_complex.c       Utility/cholmod_l_complex.c \
	Utility/cholmod_copy.c          Utility/cholmod_l_copy.c \
	Utility/cholmod_dense.c         Utility/cholmod_l_dense.c \
	Utility/cholmod_error.c         Utility/cholmod_l_error.c \
	Utility/cholmod_factor.c        Utility/cholmod_l_factor.c \
	Utility/cholmod_memory.c        Utility/cholmod_l_memory.c \
	Utility/cholmod_sparse.c        Utility/cholmod_l_sparse.c \
	Utility/cholmod_transpose.c     Utility/cholmod_l_transpose.c \
	Utility/cholmod_triplet.c       Utility/cholmod_l_triplet.c \
	Utility/cholmod_version.c       Utility/cholmod_l_version.c \
	MatrixOps/cholmod_drop.c      MatrixOps/cholmod_l_drop.c \
	MatrixOps/cholmod_horzcat.c   MatrixOps/cholmod_l_horzcat.c \
	MatrixOps/cholmod_norm.c      MatrixOps/cholmod_l_norm.c \
	MatrixOps/cholmod_scale.c     MatrixOps/cholmod_l_scale.c \
	MatrixOps/cholmod_sdmult.c    MatrixOps/cholmod_l_sdmult.c \
	MatrixOps/cholmod_ssmult.c    MatrixOps/cholmod_l_ssmult.c \
	MatrixOps/cholmod_submatrix.c MatrixOps/cholmod_l_submatrix.c \
	MatrixOps/cholmod_symmetry.c  MatrixOps/cholmod_l_symmetry.c \
	MatrixOps/cholmod_vertcat.c   MatrixOps/cholmod_l_vertcat.c \
	Modify/cholmod_rowadd.c Modify/cholmod_l_rowadd.c \
	Modify/cholmod_rowdel.c Modify/cholmod_l_rowdel.c \
	Modify/cholmod_updown.c Modify/cholmod_l_updown.c \
	Partition/cholmod_camd.c          Partition/cholmod_l_camd.c \
	Partition/cholmod_ccolamd.c       Partition/cholmod_l_ccolamd.c \
	Partition/cholmod_csymamd.c       Partition/cholmod_l_csymamd.c \
	Partition/cholmod_metis.c         Partition/cholmod_l_metis.c \
	Partition/cholmod_metis_wrapper.c \
	Partition/cholmod_nesdis.c        Partition/cholmod_l_nesdis.c \
	Supernodal/cholmod_super_numeric.c  Supernodal/cholmod_l_super_numeric.c \
	Supernodal/cholmod_super_solve.c    Supernodal/cholmod_l_super_solve.c \
	Supernodal/cholmod_super_symbolic.c Supernodal/cholmod_l_super_symbolic.c
objects = $(sources:.c=.o)
archive = CHOLMOD.a

PKG_CPPFLAGS = -I./Include -I../SuiteSparse_config \
	-I../AMD/Include -I../COLAMD/Include \
	-I../CAMD/Include -I../CCOLAMD/Include \
	-I. -I./SuiteSparse_metis/include \
	-I./SuiteSparse_metis/GKlib -I./SuiteSparse_metis/libmetis

all : $(archive)

$(archive) : $(objects)
	rm -f $@
	$(AR) -cr $@ $(objects)
	$(RANLIB) $@

.c.o :
	$(CC) $(PKG_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $*.o

clean :
	@rm -f $(objects) $(archive)
