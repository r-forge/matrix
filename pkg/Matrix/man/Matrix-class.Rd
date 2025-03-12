\name{Matrix-class}
\title{Virtual Class \sQuote{Matrix} of Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{Matrix-class}
%
\description{
  \code{Matrix} is a virtual class representing matrices.  All classes
  of matrix defined in package \pkg{Matrix} are extensions of
  \code{Matrix}.  Nonvirtual class \code{\linkS4class{matrix}} defined
  in package \pkg{methods}, from which traditional matrices inherit
  \emph{implicitly}, does \emph{not} extend \code{Matrix}, but methods
  are defined so that \code{Matrix} objects and \code{matrix} objects
  are interoperable.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{an integer vector of length 2 giving the
      dimensions of the matrix.  It is typically accessed using
      \code{dim}; see \code{\link{dim-methods}}.}
    \item{\code{Dimnames}}{a list of length 2.  Each element must be
      \code{NULL} or a character vector of length equal to the
      corresponding element of \code{Dim}.  It is typically accessed
      using \code{dimnames}; see \code{\link{dimnames-methods}}.}
  }
}
\section{Subclasses}{
  Objects inheriting from virtual class \code{Matrix} can be categorized
  in three ways: by data type, by structure, and by storage format.
  Accordingly, it can be helpful to think of nonvirtual subclasses of
  \code{Matrix} as triples
\preformatted{	(<data type>, <structure>, <storage format>)}
  and to think of virtual subclasses of \code{Matrix} as sets of triples
  satisfying constraints.

  The data type of a \code{Matrix} object refers to the set of possible
  values of entries in the represented matrix.  Package \pkg{Matrix}
  supports five data types: boolean (or \dQuote{nonzero pattern}),
  logical, integer, double, and complex.  How boolean data are stored
  varies by format, whereas logical, integer, double, and complex data
  are always stored in atomic vectors of the corresponding type.
  Virtual classes \code{\linkS4class{nMatrix}},
  \code{\linkS4class{lMatrix}}, \code{\linkS4class{iMatrix}},
  \code{\linkS4class{dMatrix}}, and \code{\linkS4class{zMatrix}}
  designate \code{Matrix} objects whose data type is one of the
  supported five (n=nonzero pattern, l=logical, i=integer, d=double,
  z=complex, following the LAPACK naming scheme).

  The structure of a \code{Matrix} object refers to the mathematical
  properties of the represented matrix.  Virtual classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{triangularMatrix}},
  \code{\linkS4class{diagonalMatrix}}, and
  \code{\linkS4class{indexMatrix}} represent general (not having any
  specific properties), Hermitian or symmetric, triangular, diagonal,
  and index matrices.  \code{\linkS4class{posdefMatrix}} is a virtual
  subclass of \code{symmetricMatrix} representing positive semidefinite
  matrices.

  The storage format of a \code{Matrix} object describes the relation
  between what is stored in memory and entries in the represented
  matrix.  A matrix can be represented in multiple formats, some more
  efficient than others depending on, e.g., the sparsity of the matrix
  and the order in which entries are accessed in computation.  Virtual
  class \code{\linkS4class{sparseMatrix}} designates \code{Matrix}
  objects whose format is such that the space required to represent an
  \eqn{m \times n}{m-by-n} matrix with \eqn{p} nonzero entries is
  \eqn{O(\max\{m, n, p\})}{O(max(m, n, p))}.  Its \dQuote{complement} is
  virtual class \code{\linkS4class{denseMatrix}}.
  \code{\linkS4class{CsparseMatrix}},
  \code{\linkS4class{RsparseMatrix}}, and
  \code{\linkS4class{TsparseMatrix}} are virtual subclasses of
  \code{sparseMatrix} corresponding to the compressed sparse column,
  compressed sparse row, and triplet formats.
  \code{\linkS4class{unpackedMatrix}} and
  \code{\linkS4class{packedMatrix}} are virtual subclasses of
  \code{denseMatrix} corresponding to the conventional and packed
  formats defined by LAPACK, storing data of length \eqn{m n}{m*n} and
  \eqn{n (n+1)/2}{n*(n+1)/2} and listing matrix entries in
  column-major order: \url{https://netlib.org/lapack/lug/node121.html}

  A canonical example of a nonvirtual subclass of \code{Matrix} is
  \code{\linkS4class{dgCMatrix}}.  Objects of this class have double
  precision data, general structure, and compressed sparse column
  format.  Accordingly, \code{dgCMatrix} extends virtual classes
  \code{dMatrix}, \code{generalMatrix}, and \code{CsparseMatrix}.
  Indeed, nonvirtual subclasses of \code{Matrix} defined in package
  \pkg{Matrix} tend to have names of the form
  \code{\var{x}\var{yy}Matrix}, where the first character indicates the
  data type and the second and third characters together indicate the
  structure and format.  This naming scheme is inspired by LAPACK:
  \url{https://netlib.org/lapack/lug/node24.html}
}
\section{Methods}{
  S4 methods are defined for the S4 generic functions listed by
  \code{\link{.S4methods}(class = "Matrix")}.

  S3 methods are defined for the S3 generic functions listed by
  \code{\link{.S3methods}(class = "Matrix")}.  They are defined so that
  \code{Matrix} objects are handled correctly even in environments where
  the S4 generic functions are not available, such as the evaluation
  frame of a function from a package that does not import the S4 generic
  function.  Notably, the S3 method defined for \code{\link{as.matrix}}
  ensures that \pkg{base} functions like \code{\link{apply}},
  \code{\link{eigen}}, and \code{\link{svd}} behave as documented even
  if their first argument is a \code{Matrix} (rather than \code{matrix})
  object.  Note that this mechanism is not needed for generic functions
  that dispatch \dQuote{internally}; see \code{\link{InternalGenerics}}.

  Methods are documented by generic function rather than by class,
  hence see:
\preformatted{	help.search(keyword = "methods", package = "Matrix")
	help("<generic>-methods", package = "Matrix")}

  Methods defined in package \pkg{Matrix} for virtual class
  \code{Matrix} or one of its virtual subclasses can be dispatched for
  subclasses of \code{Matrix} defined in other packages.  The behaviour
  is undefined if the \dQuote{foreign} subclass does not extend any
  nonvirtual \dQuote{local} subclass.  Authors of packages that define
  such foreign subclasses may need to override local methods with
  methods of their own.  An example is class \code{MultiCompanion} in
  package \CRANpkg{mcompanion} (version 0.6 at the time of writing).
}
\seealso{
  Nonvirtual subclasses such as \code{\linkS4class{dgCMatrix}} and
  \code{\linkS4class{dgeMatrix}}.  Function \code{\link{Matrix}}, a
  widely used (but not the only) constructor for \code{Matrix} objects.
}
\examples{
cl <- getClass("Matrix")
showClass(cl) # output with all slots, superclasses, subclasses

subcl <- cl@subclasses
subcl <- subcl[vapply(subcl, packageSlot, "") == "Matrix"] # exclude foreign ones
(nms <- names(subcl))
stopifnot(all(endsWith(nms, "Matrix"))) # our convention

iv <- vapply(nms, isVirtualClass, FALSE)
( v <- nms[ iv]) #    virtual subclasses
(nv <- nms[!iv]) # nonvirtual subclasses
(nv.prefix <- sub("Matrix$", "", nv))
tabulate(nchar(nv.prefix), 6L) # 3-letter prefix, usually ...
}
