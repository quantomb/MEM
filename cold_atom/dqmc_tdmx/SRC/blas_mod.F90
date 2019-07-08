module blas_mod

  ! 
  !  This module is designed for interface checking
  !

#define DB kind(1.0d0)
  
  interface

     ! y = x
     subroutine dcopy(n, x, incx, y, incy)
       integer   :: n, incx, incy
       real(DB)  :: x(*), y(*)
     end subroutine dcopy

     ! y = a*x + y
     subroutine daxpy(n, a, x, incx, y, incy)
       integer   :: n, incx, incy
       real(DB)  :: a, x(*), y(*)
     end subroutine daxpy
          
     ! dx = da*dx
     subroutine dscal(n, da, dx, incx)
       integer   :: n, incx
       real(DB)  :: da, dx(*) 
     end subroutine dscal

     ! matrix-vector multiplication: y = alpha*A*x + beta*y
     ! DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     subroutine dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
       character :: trans
       integer   :: m, n, lda, incx, incy
       real(DB)  :: alpha, beta, A(lda,*), x(*), y(*)
     end subroutine dgemv


     ! matrix-matrix multiplication
     subroutine dgemm(transA, transB, m, n, k, alpha, A, lda, &
          B, ldb, beta, c, ldc)
       character :: transA, transB
       integer   :: m, n, k, lda, ldb, ldc
       real(DB)  :: alpha, beta, A(lda,*), B(ldb,*), C(ldc,*)
     end subroutine dgemm

     ! triangular matrix-matrix multiplication
     subroutine dtrmm (side, uplo, transA, diag, m, n, alpha, &
          A, lda, B, ldb)
       character :: side, uplo, transA, diag
       integer   :: m, n, lda, ldb
       real(DB)  :: alpha, A(lda,*), B(ldb,*)
     end subroutine dtrmm
     
  end interface
  
end module blas_mod
