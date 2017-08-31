#include "math.h"
#include <stdlib.h>
#include "R.h"
#include "Rmath.h"
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <assert.h>
#include <stdio.h>
#include <R_ext/Utils.h>

/******* common functions *******/
void set_equal(int *length_x, double *x, double *copy_x) {
     /********************************************************************
     # sets copy_x equal to x.
     # Args:
     #   length_x = the length of array x (integer).
     #   x = array of length length_x.
     #   copy_x = array of length length_x that is the copy of x.
     **********************************************************************/
     int i;
     for(i = 0; i < *length_x; i++) {
         copy_x[i] = x[i];
     }
}

void mean(int *length_x, double *x, double *mean_x) {
  /********************************************************************
    # calculates the mean of x.
    # Args:
    #   length_x = the length of vector x (scalar).
    #   x = vector (length_x x 1).
    #   mean_x = scalar for storing the mean.
    **********************************************************************/
    int i;
    *mean_x = 0.0;
    for(i = 0; i < *length_x; i++) {
        *mean_x += x[i];
    }
    *mean_x /= (double) (*length_x);
}

void variance(int *length_x, double *x, double *var_x) {
  /********************************************************************
    # calculates the variance of x.
    # Args:
    #   length_x = the length of vector x (scalar).
    #   x = vector (length_x x 1).
    #   var_x = scalar for storing the variance.
    **********************************************************************/
    int i;
    double x1 = 0;
    double x2 = 0;
  
    if(*length_x == 1) {
        *var_x = 0.0;
    } else {
        for(i = 0; i < *length_x; i++) {
            x1 += x[i];
            x2 += x[i] * x[i];
    }
    *var_x = (x2 - x1 * x1 / (double) (*length_x)) / ((double) (*length_x) - 1.0);
  }
}

void find_integer_max(int *length_x, int *x, int *max_x) {
    /********************************************************************
    # finds the maximum value in a vector of integers.
    # Args:
    #   length_x: positive integer indicating the length of the vector
    #      to search through.
    #   x: the vector of integers to search through.
    #   max_x: the maximum of x.
    **********************************************************************/
    int i;
    *max_x = x[0];
    if(*length_x > 1) {
        for(i = 1; i < *length_x; i++) {
            if(x[i] > *max_x) {
                *max_x = x[i];
            }
        }
    }
}

void find_interval(int *length_x, double *x, double *tau, int *bin) {
    /********************************************************************
    # finds which bin of vector x that scalar tau is in
    # Args:
    #   length_x = length of x (scalar)
    #   x = vector of unique values that is sorted in increasing order (length_x x 1)
    #   tau = value to find in the interval (scalar)
    #   bin = updated so that x[bin] <= tau < x[bin + 1], in {0, 1, .., length_x = 1}. (scalar)
    #   Note the program will search indefinitely if tau < x[0].
    #   I used the source code of the findInterval function in R to create this code.
    **********************************************************************/
    int i_step;
    int middle = 0;
    int i_hi = *bin + 1;
  
    if(*tau < x[0]) {
        Rprintf("Lowest value of x is greater than tau.\n");
    }
  
    if (*tau <  x[ihi]) {
        if (*tau >= x[*bin]) {} /* `lucky': same interval as last time */
    else {
        /* **** now *x < x[*bin] . decrease *bin to capture *x */
        for (i_step = 1; ; i_step *= 2) {
            i_hi = *bin;
            *bin = i_hi - i_step;
        if (*bin <= 1)
            break;
        if (*tau >=  x[*bin])
            goto L50;
        }
        *bin = 1;
    }
    }
    else {
    /* **** now *tau>= * xt[ihi] . increase i_hi to capture *x*/
        for (i_step = 1; ; i_step *= 2) {
            *bin = ihi;
            i_hi = *bin + i_step;
            if (i_hi >= *length_x)
                break;
            if (*tau < x[ihi])
                goto L50;
        }
        i_hi = *length_x;
    }
    L50:
    /* **** now * x[*bin] <= *tau < *x[ihi]. narrow the interval. */
    while(middle != *bin) {
    /* note. it is assumed that middle = * bin in case i_hi = * bin + 1. */
        if (*tau >= x[middle]) {
            *bin = middle;
        } else {
            i_hi = middle;
        }
        middle = (*bin + ihi) / 2;
    }
}
    
void matrix_multiply (int *t_a, int *t_b, int *nrow_a, int *ncol_a, int *nrow_b,
    int *ncol_b, double *a, double *b, double *c){
    /********************************************************************
    # returns c = a * b
    # Args:
    #   t_a = indicator if a should be transposed
    #   t_b = indicator if b should be transposed
    #   nrow_a = number of rows in matrix a
    #   ncol_a = number of columns in matrix a
    #   nrow_b = number of rows in matrix b
    #   ncol_b = number of columns in matrix b
    #   a is matrix (nrow_a x ncol_a)
    #   b is matrix (nrow_b x ncol_b)
    #   c is matrix of dimension:
    nrow_a x ncol_b if t_a = 0 and t_b = 0
    ncol_a x ncol_b if t_a = 1 and t_b = 0
    nrow_a x nrow_b if t_a = 0 and t_b = 1
    ncol_a x nrow_b if t_a = 1 and t_b = 1
    **********************************************************************/
    int i;
    double dummy_one = 1.0;
    
    int nrow_c = *nrow_a;
    int ncol_c = *ncol_b;
    int mid_c = *ncol_a;
    
    char transa = 'N';
    if (*t_a == 1){
        nrow_c = *ncol_a;
        mid_c = *nrow_a;
        transa = 'T';
    }
    char transb = 'N';
    if (*t_b == 1){
        ncol_c = *nrow_b;
        transb = 'T';
    }
    int dim_c = nrow_c * ncol_c;
    for(i = 0; i < dim_c; i++) {
        c[i] = 0.0;
    }
    /*
    lda informs BLAS how large a "stride" to take inside the matrix
    C language is row major so stride is the number of columns
    */
    F77_CALL(dgemm)(&transa, &transb, &nrow_c, &ncol_c, &mid_c, &dummy_one, a, nrow_a, b, nrow_b, &dummy_one, c, &nrow_c);
}
    
    void Kronecker (int *nrow_a, int *ncol_a, int *nrow_b, int *ncol_b, double *a,
    double *b, double *c){
    /********************************************************************
    # produces c = kronecker(a, b)
    # Args:
    #   nrow_a = number of rows of a
    #   ncol_a = number of columns of a
    #   nrow_b = number of rows of b
    #   ncol_b = number of columns of b
    #   a = first matrix (nrow_a x ncol_a)
    #   b = second matrix (nrow_b x ncol_b)
    #   c = matrix of dimension nrow_a * nrow_b x ncol_a * ncol_b
    **********************************************************************/
    int a1, a2, b1, b2;
    for(a1 = 0; a1 < *nrow_a; a1++){
    for(a2 = 0; a2 < *ncol_a; a2++){
    for(b1 = 0; b1 < *nrow_b; b1++){
    for(b2 = 0; b2 < *ncol_b; b2++){
    c[(*nrow_b * a1 + b1)  + (*nrow_a * *nrow_b) * (*ncol_b * a2 + b2)] =  b[b1 + *nrow_b * b2] *  a[a1 + *nrow_a * a2];
    }
    }
    }
    }
    }
    
    void QuadraticForm (int *length_x, double *x, double *a, double *qf){
    /********************************************************************
    # calculates qf = x'a x, where a is symmetric
    # Args:
    #   length_x = the length of vector x (scalar).
    #   x = vector (length_x x 1).
    #   a = matrix (length_x x length_x).
    #   qf = quadratic form (scalar).
    **********************************************************************/
      int i, j;
    
    if(*length_x == 1) {
      *qf = x[0] * x[0] * a[0];
    } else {
      *qf = 0.0;
      for(i = 1; i < *length_x; i++) {
        for(j = 0; j < i; j++) {
          *qf += x[i] * x[j] * a[i + *length_x * j];
        }
      }
      *qf *= 2.0;
      for(i = 0; i < *length_x; i++) {
        *qf += x[i] * x[i] * a[i + *length_x * i];
      }
    }
    }
  
  void QuadraticFormKronecker (int *nrow_d, int *nrow_e, double *d, double *e,
                               double *x, double *qf){
    /********************************************************************
      # This function returns x ' (kronecker(d, e)) x
      #   without calculating the kronecker product
      #
      # Args:
      #  nrow_d: number of rows of square matrix d.
      #  nrow_e: number of rows of square matrix e.
      #  d: square matrix with nrow_d rows.
      #  e: square matrix with nrow_e rows.
      #  x: vector of length nrow_d * nrow_e.
      #  qf: scalar for storing the quadratic form.
      **********************************************************************/
      int i;
    int t_a = 0;
    int t_b = 0;
    int length_exd = *nrow_e * *nrow_d;
    
    double *ex = (double *)Calloc(length_exd, double);
    double *exd = (double *)Calloc(length_exd, double);
    
    for(i = 0; i < length_exd; i++) {
      ex[i] = 0.0;
      exd[i] = 0.0;
    }
    
    matrix_multiply(&t_a, &t_b, nrow_e, nrow_e, nrow_e, nrow_d, e, x, ex);
    matrix_multiply(&t_a, &t_b, nrow_e, nrow_d, nrow_d, nrow_d, ex, d, exd);
    
    *qf = 0.0;
    for(i = 0; i < length_exd; i++) {
      *qf += x[i] * exd[i];
    }
    
    Free(ex);
    Free(exd);
  }
  
  void CholMatrix(int *nrow_v, double *v, double *v_chol, double *log_det){
    /********************************************************************
      # finds the lower triangular Cholesky decomposition of V
      # Args:
      #   nrow_v = number of rows of v (scalar)
      #   v = positive definite symmetric matrix (nrow_v x nrow_v)
      #   v_chol = updated to the array with chol(v) (nrow_v x nrow_v)
      #   log_det = updated to the log determinant of v (scalar)
      **********************************************************************/
      char uplo = 'L';
      int i, j, info;
      int nrow_v_2 = *nrow_v * *nrow_v;
      for(i = 0; i < nrow_v_2; i++){
        v_chol[i] = v[i];
      }
      /*cholesky v*/
        F77_CALL(dpotrf)(&uplo, nrow_v, v_chol, nrow_v, &info);
      if (info) {
        Rprintf("Error with chol(mat_chol): info = %d\n", info);
      }
      /*calculate the log determinant*/
        *log_det = log(v_chol[0]);
        for(i = 1; i < *nrow_v; i++){
          *log_det += log(v_chol[i + *nrow_v * i]);
        }
        *log_det *= 2;
        /*ensure the upper portion of the triangle is 0*/
          for(i = 0; i < *nrow_v; i++){
            for(j = i + 1; j < *nrow_v; j++){
              v_chol[i + *nrow_v * j] = 0.0;
            }
          }
  }
  
  void InvertSymmetricMatrix(int *nrow_v, double *v, double *v_inv,
                             double *log_det){
    /********************************************************************
      # finds inverse of positive definite symmetric matrix and log determinant of the inverse
      #
      # Args:
      #   nrow_v = number of rows of v (scalar)
      #   v = positive definite symmetric matrix (nrow_v x nrow_v)
      #   v_inv = updated to the nrow_v x nrow_v array with the inverted matrix (nrow_v x nrow_v)
      #   log_det = updated to the log determinant of v_inv (scalar)
      **********************************************************************/
      char uplo = 'L';
      int i, j, info;
      /*cholesky dummy_inv*/
        CholMatrix(nrow_v, v, v_inv, log_det);
      *log_det *= -1;
      /*solve for the upper part of the triangle*/
        F77_CALL(dpotri)(&uplo, nrow_v, v_inv, nrow_v, &info);
      for(i = 0; i < *nrow_v; i++){
        for(j = i + 1; j < *nrow_v; j++){
          v_inv[i + *nrow_v * j] = v_inv[j + *nrow_v * i];
        }
      }
  }
  
  void MakeAbsoluteDistance(int *n_timepoints, double *timepoints,
                            double *distance_matrix) {
    /******************************************************
      # Construct an absolute distance matrix.
      # Args:
      #   n_timepoints: the number of timepoints.
      #   timepoints: vector of unique timepoints.
      #   distance_matrix:
      n_timepoints x n_timepoints array for storing the distances.
    ******************************************************/
      int i, j;
    for (i = 0; i < *n_timepoints; i++) {
      for (j = 0; j < *n_timepoints; j++) {
        distance_matrix[i + *n_timepoints * j] =
          fabs(timepoints[i] - timepoints[j]);
      }
    }
  }
  
  void MakeDiagonalPrecision(int *n_timepoints, double *rho, double *dist,
                             double *prec_ar, double *log_det) {
  }
  
  void MakeAutoregressivePrecision(int *n_timepoints, double *rho, double *dist,
                                   double *prec_ar, double *log_det) {
    /******************************************************
      # Construct an AR-1 precision matrix.
      # Args:
      #   n_timepoints: the number of timepoints.
      #   rho: AR-1 correlation matrix.
      #   dist: matrix indicating the distance between timepoints.
      #   prec_ar: n_timepoints x n_timepoints array for storing the precision.
      #   log_det: scalar for storing the log determinant of the precision matrix.
      ******************************************************/
      int i, j;
    double *corr_mat = (double *)Calloc(*n_timepoints * *n_timepoints, double);
    for (i = 0; i < *n_timepoints; i++) {
      for (j = 0; j < *n_timepoints; j++) {
        corr_mat[i + *n_timepoints * j] = pow(*rho, dist[i + *n_timepoints * j]);
      }
    }
    InvertSymmetricMatrix(n_timepoints, corr_mat, prec_ar, log_det);
    Free(corr_mat);
  }
  
  void MakeSpatialPrecision(int *n_locs, double *rho, double *dist,
                            double *prec_sp, double *log_det) {
    /******************************************************
      # Construct a spatial exponential precision matrix.
      # Args:
      #   n_locs: the number of locations.  Locations must be unique.
      #   rho: range parameter (positive scalar).
      #   dist: euclidean distance matrix that is n_locs x n_locs.
      #   prec_sp: n_locs x n_locs array for storing the precision.
      #   log_det: scalar for storing the log determinant of the precision matrix.
      ******************************************************/
      int i, j;
    double *corr_mat = (double *)Calloc(*n_locs * *n_locs, double);
    
    for (i = 0; i < *n_locs; i++) {
      for (j = 0; j < *n_locs; j++) {
        corr_mat[i + *n_locs * j] = exp(-1.0 * dist[i + *n_locs * j] / *rho);
      }
    }
    InvertSymmetricMatrix(n_locs, corr_mat, prec_sp, log_det);
    Free(corr_mat);
  }
  
  void MakeSplineKnots(int *interior_knots_length, int *spline_df,
                       double *interior_knots, double *spline_knots) {
    /******************************************************
      # Constructs knot sequence for cubic spline basis.
      # Args:
      # interior_knots_length: positive integer indicating the length of the interior knots.
      # spline_df: positive integer indicating the spline degrees of freedom.
      # interior_knots: vector of knots that are in (0, 1).
      # spline_knots: vector of length (2 * spline_df + interior_knots_length) for storing the spline_knots.
      ******************************************************/
      int i;
    
    for(i = 0; i < *spline_df; i++) {
      spline_knots[i] = 0.0;
    }
    for(i = *spline_df; i < (*spline_df + *interior_knots_length); i++) {
      spline_knots[i] = interior_knots[i - *spline_df];
    }
    for(i = (*spline_df + *interior_knots_length); i < (2 * *spline_df + *interior_knots_length + 2); i++) {
      spline_knots[i] = 1.0;
    }
  }
  
  void MSpline(double *tau, int *spline_df, double *spline_knots, int *m,
               double *v) {
    /******************************************************
      # Given knots spline.knots, evaluates the mth m-spline of
      #   spline_df degrees of freedom at tau.
      # Args:
      #   tau: scalar in (0, 1).
      #   spline_df: degrees of freedom for splines.
      #   spline_knots: knot sequence.
      #   m: index for the spline.
      #   v: scalar for updating.
      ******************************************************/
      
      int spline_df_1 = *spline_df - 1;
      double v1 = 0.0;
      double v2 = 0.0;
      double d1 = (*tau - spline_knots[*m]);
      double d2 = (spline_knots[*m + *spline_df] - *tau);
      int m_1 = *m + 1;
      
      if(*tau < spline_knots[*m] || *tau >= spline_knots[*m + *spline_df]) {
        *v = 0.0;
      } else if (*spline_df == 1) {
        *v = 1.0 / (spline_knots[*m + 1] - spline_knots[*m]);
      } else {
        MSpline(tau, &spline_df_1, spline_knots, m, &v1);
        MSpline(tau, &spline_df_1, spline_knots, &m_1, &v2);
        *v = d1 * v1 + d2 * v2;
        *v *= *spline_df / ((*spline_df - 1) * (spline_knots[*m + *spline_df] - spline_knots[*m]));
      }
  }
  
  void ISpline(double *tau, int *spline_df, int *spline_knots_length,
               double *spline_knots, int *m,  double *v) {
    /******************************************************
      # Given knots spline.knots, evaluates the mth i-spline of
      #   spline_df degrees of freedom at tau.
      # Args:
      #   tau: scalar in (0, 1).
      #   spline_df: degrees of freedom for splines.
      #   spline_knots_length: length of the knot sequence.
      #   spline_knots: knot sequence.
      #   m: index for the spline.
      #   v: scalar for updating.
      ******************************************************/
      
      int spline_df_1 = *spline_df + 1;
      double v1 = 0.0;
      int l;
      
      int bin = 2;
      find_interval(spline_knots_length, spline_knots, tau, &bin);
      
      if(bin < *m) {
        *v = 0.0;
      } else if (bin - *spline_df + 1 > *m) {
        *v = 1.0;
      } else {
        *v = 0.0;
        for(l = *m; l < bin + 1; l++) {
          MSpline(tau, &spline_df_1, spline_knots, &l, &v1);
          *v += v1 * (spline_knots[l + *spline_df + 1] - spline_knots[l]) / spline_df_1;
        }
      }
  }
  
  void CubicMSpline(double *tau, int *spline_knots_length, double *spline_knots,
                    int *m, int *bin, double *v){
    /******************************************************
      # Given knots spline_knots, evaluates the mth cubic m-spline of
      #   3 degrees of freedom at tau.
      # Args:
      #   tau: scalar in (0, 1).
      #   m: indicator for element of the basis.
      #   spline_knots_length: length of the vector of spline_knots.
      #   spline_knots: knots for the cubic spline basis.
      #   bin: index for the interval in which tau resides.
      #   v: Given spline_knots and 3 degrees of freedom,
      #       the mth m-spline evaluated at tau.
      ******************************************************/
      find_interval(spline_knots_length, spline_knots, tau, bin);
    *v = 0.0;
    if (*bin == *m){
      *v = 3 * pow(*tau - spline_knots[*m], 2) /
        ((spline_knots[*m + 1] - spline_knots[*m]) *
           (spline_knots[*m + 2] - spline_knots[*m]) *
           (spline_knots[*m + 3] - spline_knots[*m])
        );
    }
    else if(*bin == (*m + 1)){
      *v = 3 / ((spline_knots[*m + 2] - spline_knots[*m + 1]))  -
        3 *(pow(spline_knots[*m + 3] - *tau, 2) / ((spline_knots[*m + 3] - spline_knots[*m + 1]) * (spline_knots[*m + 3] - spline_knots[*m]) * (spline_knots[*m + 2] - spline_knots[*m + 1]))
            + pow(*tau - spline_knots[*m], 2) / ((spline_knots[*m + 3] - spline_knots[*m]) * (spline_knots[*m + 2] - spline_knots[*m]) * (spline_knots[*m + 2] - spline_knots[*m + 1]))
        );
    }
    else if(*bin == (*m + 2)){
      *v = 3 * pow(spline_knots[*m + 3] - *tau,2) / ((spline_knots[*m + 3] - spline_knots[*m + 2]) * (spline_knots[*m + 3] - spline_knots[*m + 1]) * (spline_knots[*m + 3] - spline_knots[*m]));
    }
  }
  
  void CubicISpline(double *tau, int *spline_knots_length, double *spline_knots,
                    int *m, int *bin, double *v) {
    /******************************************************
      # Given knots spline_knots, evaluates the mth cubic i-spline of
      #   3 degrees of freedom at tau.
      # Args:
      #   tau: scalar in (0, 1).
      #   m: indicator for element of the basis.
      #   spline_knots_length: length of the vector of spline_knots.
      #   spline_knots: knots for the cubic spline basis.
      #   bin: index for the interval in which tau resides.
      #   v: Given spline_knots and 3 degrees of freedom,
      #       the mth m-spline evaluated at tau.
      ******************************************************/
      *v = 0.0;
      find_interval(spline_knots_length, spline_knots, tau, bin);
      if ((*bin - 2) > *m) {
        *v = 1.0;
      }
      else if (*bin == *m) {
        *v = pow((*tau - spline_knots[*m]), 3) /
          ((spline_knots[*m + 1]  - spline_knots[*m]) *
             (spline_knots[*m + 2]  - spline_knots[*m]) *
             (spline_knots[*m + 3]  - spline_knots[*m]));
      }
      else if (*bin == (*m + 1)) {
        double i1 =  (*tau - spline_knots[*m]) / (spline_knots[*m + 2] - spline_knots[*m + 1]);
        double i2 =  (pow(*tau - spline_knots[*m + 1], 2) - pow(spline_knots[*m + 3] - *tau, 2) ) /
          ((spline_knots[*m + 3] - spline_knots[*m + 1]) * (spline_knots[*m + 2] - spline_knots[*m + 1]));
        double i3 =   pow(spline_knots[*m + 3] - *tau,3) / ((spline_knots[*m + 3] - spline_knots[*m + 1]) *
                                                              (spline_knots[*m + 3] - spline_knots[*m]) *
                                                              (spline_knots[*m + 2] - spline_knots[*m + 1])) -
          pow(*tau - spline_knots[*m] , 3) /
          ((spline_knots[*m + 3] - spline_knots[*m]) *
             (spline_knots[*m + 2] - spline_knots[*m]) *
             (spline_knots[*m + 2] - spline_knots[*m + 1]));
        *v = i1 + i2 + i3;
      }
      else if (*bin == (*m + 2)) {
        *v = 1 - pow(spline_knots[*m + 3]- *tau,3) / ((spline_knots[*m + 3] - spline_knots[*m + 2]) *
                                                        (spline_knots[*m + 3] - spline_knots[*m + 1]) *
                                                        (spline_knots[*m + 3]- spline_knots[*m]));
      }
  }
  
  void MakeCubicSplineCoefficients (int *n_bf, double *spline_knots,
                                    double *spline_coefs) {
    /******************************************************
      # Makes the coefficients used in CubicISpline2.
      #   3 degrees of freedom at tau using precalculated spline coefficients.
      # Args:
      #   n_bf: integer number of basis functions.
      #   spline_knots: knots for the cubic spline basis.
      #   spline_coefs: coefficients for the cubic spline basis.
      ******************************************************/
      
      int dummy_bin, l, m;
    int spline_knots_length = *n_bf + 4;
    int M_1 = *n_bf - 1;
    
    double *d_a = (double *)Calloc(spline_knots_length, double);
    double *d_b = (double *)Calloc(spline_knots_length * 4, double);
    double *d_c = (double *)Calloc(spline_knots_length, double);
    
    
    for (m = 0; m < spline_knots_length; m++) {
      d_a[m] = 0.0; d_c[m] = 0.0;
    }
    for (m = 0; m < spline_knots_length * 4; m++) {
      d_b[m] = 0.0;
    }
    
    for (m = 2; m < M_1; m++) {
      d_a[m] = 1 / ((spline_knots[m + 1] - spline_knots[m]) * (spline_knots[m + 2] - spline_knots[m]) * (spline_knots[m + 3] - spline_knots[m]));
    }
    for (dummy_bin = 2; dummy_bin < M_1; dummy_bin++) {
      m = dummy_bin;
      d_b[dummy_bin]           = 1.0 / (spline_knots[m + 1] - spline_knots[m]);
      d_b[dummy_bin + spline_knots_length * 1] = 1.0 / ((spline_knots[m + 2] - spline_knots[m]) * (spline_knots[m + 1] - spline_knots[m]));
    }
    for (dummy_bin = 2; dummy_bin < M_1; dummy_bin++) {
      m = dummy_bin - 1;
      d_b[dummy_bin + spline_knots_length * 2] = -1.0 / ((spline_knots[m + 3] - spline_knots[m + 1]) * (spline_knots[m + 3] - spline_knots[m]) * (spline_knots[m + 2] - spline_knots[m + 1]));
      d_b[dummy_bin + spline_knots_length * 3] = -1.0 / ((spline_knots[m + 3] - spline_knots[m])     * (spline_knots[m + 2] - spline_knots[m]) * (spline_knots[m + 2] - spline_knots[m + 1]));
    }
    for (dummy_bin = 2; dummy_bin < M_1; dummy_bin++) {
      m = dummy_bin - 2;
      d_c[dummy_bin] = 1.0 / ((spline_knots[m + 3] - spline_knots[m + 2]) * (spline_knots[m + 3] - spline_knots[m + 1]) * (spline_knots[m + 3] - spline_knots[m]));
    }
    
    for (m = 0; m < 4 * *n_bf * M_1; m++) {
      spline_coefs[m] = 0.0;
    }
    for (l = 1; l < *n_bf; l++) {
      for (m = 2; m < M_1; m++) {
        if (m > (l + 1)) {
          spline_coefs[3 + 4 * (l + *n_bf * m)] = 1.0;
        }
        if (l == m + 1) {
          spline_coefs[0 + 4 * (l + *n_bf * m)] =      d_a[m];
          spline_coefs[1 + 4 * (l + *n_bf * m)] = -3.0 * d_a[m] * spline_knots[m];
          spline_coefs[2 + 4 * (l + *n_bf * m)] =  3.0 * d_a[m] * spline_knots[m] * spline_knots[m];
          spline_coefs[3 + 4 * (l + *n_bf * m)] = -1.0 * d_a[m] * spline_knots[m] * spline_knots[m] * spline_knots[m];
        }
        if (l == (m)) {
          spline_coefs[0 + 4 * (l + *n_bf * m)] =      d_b[m + spline_knots_length * 2] + d_b[m + spline_knots_length * 3];
          spline_coefs[1 + 4 * (l + *n_bf * m)] = -3.0 * d_b[m + spline_knots_length * 2] * spline_knots[m + 2]                       - 3 * d_b[m + spline_knots_length * 3] * spline_knots[m - 1];
          spline_coefs[2 + 4 * (l + *n_bf * m)] =  3.0 * d_b[m + spline_knots_length * 2] * spline_knots[m + 2] * spline_knots[m + 2] + 3 * d_b[m + spline_knots_length * 3] * spline_knots[m - 1] * spline_knots[m - 1] + d_b[m + spline_knots_length * 0] + 2.0 * d_b[m + spline_knots_length * 1] * (spline_knots[m + 2] - spline_knots[m]);
          spline_coefs[3 + 4 * (l + *n_bf * m)] = -1.0 * d_b[m + spline_knots_length * 2] * pow(spline_knots[m + 2], 3) - d_b[m + spline_knots_length * 3] * pow(spline_knots[m - 1], 3) - d_b[m + spline_knots_length * 0] * spline_knots[m - 1] + d_b[m + spline_knots_length * 1] * (pow(spline_knots[m],2) - pow(spline_knots[m + 2], 2));
        }
        if (l == (m - 1)) {
          spline_coefs[0 + 4 * (l + *n_bf * m)] =      d_c[m];
          spline_coefs[1 + 4 * (l + *n_bf * m)] = -3 * d_c[m] * spline_knots[m + 1];
          spline_coefs[2 + 4 * (l + *n_bf * m)] =  3 * d_c[m] * spline_knots[m + 1] * spline_knots[m + 1];
          spline_coefs[3 + 4 * (l + *n_bf * m)] = -1 * d_c[m] * pow(spline_knots[m + 1], 3) + 1;
        }
      }
    }
    
    for(m = 0; m < M_1; m++) {
      spline_coefs[3 + 4 * *n_bf * m] = 1.0;
    }
    Free(d_a);
    Free(d_b);
    Free(d_c);
  }
  
  void CubicMSpline2(double *tau, int *spline_knots_length, double *spline_knots,
                     double *spline_coefs, int *m, int *bin, double *v) {
    /******************************************************
      # Given knots spline_knots, evaluates the mth cubic m-spline of
      #   3 degrees of freedom at tau using precalculated spline coefficients.
      # Args:
      #   tau: scalar in (0, 1).
      #   spline_knots_length: length of the vector of spline_knots.
      #   spline_knots: knots for the cubic spline basis.
      #   spline_coefs: coefficients for the cubic spline basis.
      #   m: indicator for element of the basis.
      #   bin: index for the interval in which tau resides.
      #   v: Given spline_knots and 3 degrees of freedom,
      #       the mth m-spline evaluated at tau.
    ******************************************************/
      int n_bf = *spline_knots_length - 4;
      find_interval(spline_knots_length, spline_knots, tau, bin);
      int m_1 = *m + 1;
      
      *v = 3.0 * spline_coefs[0 + 4 * (m_1 + n_bf * *bin)] *  *tau * *tau +
        2.0 * spline_coefs[1 + 4 * (m_1  + n_bf * *bin)] * *tau +
        spline_coefs[2 + 4 * (m_1  + n_bf * *bin)];
      
  }
  
  void CubicISpline2(double *tau, int *spline_knots_length, double *spline_knots,
                     double *spline_coefs, int *m, int *bin, double *v) {
    /******************************************************
      # Given knots spline_knots, evaluates the mth cubic i-spline of
      #   3 degrees of freedom at tau using precalculated spline coefficients.
      # Args:
      #   tau: scalar in (0, 1).
      #   spline_knots_length: length of the vector of spline_knots.
      #   spline_knots: knots for the cubic spline basis.
      #   spline_coefs: coefficients for the cubic spline basis.
      #   m: indicator for element of the basis.
      #   bin: index for the interval in which tau resides.
      #   v: Given spline_knots and 3 degrees of freedom,
      #       the mth i-spline evaluated at tau.
    ******************************************************/
      int n_bf = *spline_knots_length - 4;
      find_interval(spline_knots_length, spline_knots, tau, bin);
      int m_1 = *m + 1;
      
      *v = spline_coefs[0 + 4 * (m_1 + n_bf * *bin)] * pow(*tau, 3) +
        spline_coefs[1 + 4 * (m_1  + n_bf * *bin)] * pow(*tau, 2) +
        spline_coefs[2 + 4 * (m_1  + n_bf * *bin)] * *tau +
        spline_coefs[3 + 4 * (m_1  + n_bf * *bin)];
      
  }
  
  void MakeStickBreakingWeights(int *length_v, double *v, double *nu){
    /**************************************************
      # Transforms a vector of beta random variables and a 1 at the end into
      #  latent class weights.
      #
      # Args:
      #   length_v: integer indicating the number of latent classes.
      #   v: vector containing beta random variables,
      #      except for the final element which is 1.
      #
      # Updates:
      #   nu : Latent class weights.
      ******************************************************/
      int i;
    double cumprobs = (1 - v[0]);
    nu[0] = v[0];
    if(*length_v > 1) {
      for(i = 1; i < *length_v; i++){
        nu[i] = v[i] * cumprobs;
        cumprobs *= (1 - v[i]);
      }
    }
  }
  
  void LogDNorm (double *y, double *loc, double *scale, double *shape,
                 double *d_y){
    /********************************************************************
      # finds the log density of the normal distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape (shape is 0 for normal)
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter (scalar)
      #   d_y = density (scalar).
      **********************************************************************/
      *d_y = dnorm(*y, *loc, *scale, 1);
  }
  
  void LogDT (double *y, double *loc, double *scale, double *shape,
              double *d_y){
    /********************************************************************
      # finds the log density of the t distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *d_y =  dt((*y - *loc) / *scale, *shape, 1) - log(*scale);
  }
  
  void LogDAsymmetricLaplace (double *y, double *loc, double *scale,
                              double *shape, double *d_y){
    /********************************************************************
      # finds the log density of the asymmetric Laplace distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter in (0, 1) (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double abs_diff = *y - *loc;
      double kappa = sqrt(*shape / (1.0 - *shape));
      double kappa2 = kappa * kappa;
      *d_y = 0.5 * log(2) - log(*scale) + log(kappa) - log1p(kappa2);
      if (abs_diff < 0.0){
        kappa = 1.0 / kappa;
        abs_diff *= -1;
      }
      *d_y -= sqrt(2) / *scale * abs_diff * kappa;
  }
  
  void LogDLogistic (double *y, double *loc ,double *scale, double *shape,
                     double *d_y){
    /********************************************************************
      # finds the log density of the asymmetric logistic distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape (shape is ignored)
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double z = (*y - *loc)/ *scale;
      *d_y = (-z - 2 * log(1 + exp(-z)) - log(*scale));
  }
  
  void LogDWeibull (double *y, double *loc ,double *scale, double *shape,
                    double *d_y){
    /********************************************************************
      # finds the log density of the Weibull distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *d_y = dweibull(*y - *loc, *shape, *scale, 1);
  }
  
  void LogDGamma (double *y, double *loc ,double *scale, double *shape,
                  double *d_y){
    /********************************************************************
      # finds the log density of the gamma distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *d_y = dgamma(*y - *loc, *shape, *scale, 1);
  }
  
  void LogDBeta(double *y, double *shape1, double *shape2, double *d_y){
    /********************************************************************
      # finds the log density of the beta distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for the density (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *d_y = dbeta(*y, *shape1, *shape2, 1);
  }
  
  void LogDInverseWishart(int *p, double *df,
                          double *sai, double *x, double *d_y) {
    /********************************************************************
      # evaluates the log density of the inverse wishart distribution evaluated at x
      #   with scale matrix sai and df degrees of freedom.
      # Args:
      #   p = number of rows of x.
      #   df = degrees of freedom.
      #   sai = scale matrix.
      #   x = positive definite matrix argument.
      #   d_y = log density.
      **********************************************************************/
      int i;
    int p2 = *p * *p;
    int int_zero = 0;
    double log_det = 0.0;
    double log_det_sai = 0.0;
    double trace = 0.0;
    double *x_inv = (double *) Calloc(p2, double);
    double *product_matrix = (double *) Calloc(p2, double);
    double double_p = (double) (*p);
    double double_i;
    double this_double;
    double ff = 0.0;
    for(i = 0; i < p2; i++) {
      x_inv[i] = 0.0;
      product_matrix[i] = 0.0;
    }
    InvertSymmetricMatrix(p, x, x_inv, &log_det);
    log_det *= -1.0;
    matrix_multiply(&int_zero, &int_zero, p, p, p, p, sai, x_inv, product_matrix);
    for(i = 0; i < *p; i++) {
      trace += product_matrix[i + *p * i];
      double_i = (double) (i);
      this_double = 0.5 * (*df - double_i);
      ff += lgamma(this_double);
    }
    ff *= -1.0;
    /*calculate log determinant of sai*/
      InvertSymmetricMatrix(p, sai, x_inv, &log_det_sai);
    log_det_sai *= -1.0;
    ff += 0.5 * *df * log_det_sai;
    ff -= 0.5 * *df * double_p * log(2);
    ff -= 0.25 * double_p * (double_p - 1.0) * log(M_PI);
    *d_y = -0.5 * ((double_p + *df + 1.0) * log_det + trace) + ff;
    Free(x_inv);
    Free(product_matrix);
  }
  
  void PNorm (double *y, double *loc, double *scale, double *shape,
              double *p_y){
    /********************************************************************
      # finds the cdf of the normal distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      #   (shape is 0 for normal).
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *p_y = pnorm(*y, *loc, *scale, 1, 0);
  }
  
  void PT (double *y, double *loc, double *scale, double *shape, double *p_y){
    /********************************************************************
      # finds the cdf of the t distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double z = (*y - *loc)/ *scale;
      *p_y = pt(z, *shape, 1, 0);
  }
  
  void PLogistic (double *y, double *loc, double *scale, double *shape,
                  double *p_y){
    /********************************************************************
      # finds the cdf of the logistic distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape (shape is ignored)
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      *p_y = pow((1.0 + exp((*loc - *y) / *scale)), -1);
  }
  
  void PAsymmetricLaplace(double *y, double *loc, double *scale, double *shape,
                          double *p_y){
    /********************************************************************
      # finds the cdf of the asymmetric Laplace distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter in (0, 1) (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double temp5;
    double kappa = sqrt(*shape / (1.0 - *shape));
    double kappa2 = kappa * kappa;
    temp5 = (kappa2) / (1.0 + kappa2);
    double exponent = -1.0 * (sqrt(2.0) / *scale) * (*y - *loc) * kappa;
    if(*y < *loc){
      exponent /= -1.0 * kappa2;
    }
    temp5 = exp(exponent) / (1 + kappa * kappa);
    if(*y < *loc){
      *p_y = kappa2 * temp5;
    }
    else{
      *p_y = 1.0 - temp5;
    }
  }
  
  void PWeibull (double *y, double *loc ,double *scale, double *shape,
                 double *p_y){
    /********************************************************************
      # finds the cdf of the Weibull distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double z = (*y - *loc);
      if(z <= 0) {
        *p_y = 0.0;
      } else {
        *p_y = pweibull(z, *shape, *scale, 1, 0);
      }
  }
  
  void PGamma (double *y, double *loc ,double *scale, double *shape,
               double *p_y){
    /********************************************************************
      # finds the cdf of the gamma distribution evaluated at y
      #   with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   y = argument for cdf (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   p_y = cdf (scalar).
      **********************************************************************/
      double z = (*y - *loc);
      if(z <= 0) {
        *p_y = 0.0;
      } else{
        *p_y = pgamma(z, *shape, *scale, 1, 0);
      }
  }
  
  void QNorm(double *tau, double *loc, double *scale, double *shape,
             double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the normal distribution
      #     with mean/scale/shape parameters loc/scale/shape (shape is 0 for normal)
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      *q_tau = *loc + *scale * qnorm(*tau, 0, 1, 1, 0);
  }
  
  void QT(double *tau, double *loc, double *scale, double *shape,
          double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the t distribution
      #     with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      *q_tau = *loc + *scale * qt(*tau, *shape, 1, 0);
  }
  
  void QLogistic(double *tau, double *loc, double *scale, double *shape,
                 double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the standard logistic distribution
      #     with mean/scale/shape parameters loc/scale/shape (shape is 0 for logistic)
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      *q_tau = *loc + *scale * (log(*tau) - log(1 - *tau));
  }
  
  void QAsymmetricLaplace (double *tau, double *loc, double *scale,
                           double *shape, double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the asymmetric Laplace distribution
      #     with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = shape parameter in (0, 1) (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      double kappa = sqrt(*shape / (1.0 - *shape));
      double kappa2 = kappa * kappa;
      double temp5 = (kappa2) / (1.0 + kappa2);
      if(*tau <= temp5){
        *q_tau = *loc + (*scale * kappa) * log(*tau / temp5) / pow(2, .5);
      }
      else{
        *q_tau = *loc - (*scale / kappa) *
          (log1p(kappa2) + log1p(-1.0 * *tau)) / sqrt(2.0);
      }
  }
  
  void QWeibull(double *tau, double *loc, double *scale, double *shape,
                double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the Weibull distribution
      #     with mean/scale/shape parameters loc/scale/shape
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      *q_tau = *loc + qweibull(*tau, *shape, *scale, 1, 0);
  }
  
  void QGamma(double *tau, double *loc ,double *scale, double *shape,
              double *q_tau){
    /********************************************************************
      # finds the tauth quantile of the gamma distribution
      #     with mean/scale/shape parameters loc/scale/shape (loc parameter is not used)
      # Args:
      #   tau = quantile level in (0, 1) (scalar)
      #   loc = location parameter (scalar)
      #   scale = positive scale parameter (scalar)
      #   shape = positive shape parameter (scalar)
      #   q_tau = quantile (scalar).
      **********************************************************************/
      *q_tau = *loc + qgamma(*tau, *shape, *scale, 1, 0);
  }
  
  void TruncatedGaussianMoments(int *type, double *mu, double *sigma, double *a,
                                double *b, double *mn, double *vr){
    /**************************************************
      # Finds the mean and variance of a truncated normal random variable.
      #
      # Args:
      #   type: indicator for truncation type, with 1/2/3 corresponding to left/right/interval.
      #   mu: the mean parameter.
      #   sigma: the scale parameter.
      #   a: the left truncation point.  -99 corresponds to -Inf.
      #   b: the right truncation point.  -99 corresponds to Inf.
      #   mn: mean.
      #   vr: variance.
      #
    ******************************************************/
      double a_centered = (*a - *mu) / *sigma;
      double b_centered = (*b - *mu) / *sigma;
      double sigma2 = *sigma * *sigma;
      double dn = 0.0;
      double dn2 = 0.0;
      double pn = 0.0;
      double pn2 = 0.0;
      double lambda = 0.0;
      double delta = 0.0;
      double r;
      double dummy_zero = 0.0;
      double dummy_one = 1.0;
      double dummy_shape = 0.0;
      
      
      
      if(*type == 1) {
        LogDNorm(&a_centered, &dummy_zero, &dummy_one, &dummy_shape, &dn);
        dn = exp(dn);
        PNorm(&a_centered, &dummy_zero, &dummy_one, &dummy_shape, &pn);
        if(pn == 1.0) {
          Rprintf("pn is numerically 1.\n");
        }
        lambda = dn / (1.0 - pn);
        delta = lambda * (lambda - a_centered);
        *mn = *mu + *sigma * lambda;
        *vr = sigma2 * (1 - delta);
      } else if(*type == 2) {
        LogDNorm(&b_centered, &dummy_zero, &dummy_one, &dummy_shape, &dn);
        dn = exp(dn);
        PNorm(&b_centered, &dummy_zero, &dummy_one, &dummy_shape, &pn);
        lambda = dn / pn;
        *mn = *mu - *sigma * lambda;
        *vr = sigma2 * (1 - b_centered * lambda - lambda * lambda);
      } else {
        LogDNorm(&a_centered, &dummy_zero, &dummy_one, &dummy_shape, &dn);
        dn = exp(dn);
        LogDNorm(&b_centered, &dummy_zero, &dummy_one, &dummy_shape, &dn2);
        dn2 = exp(dn2);
        PNorm(&a_centered, &dummy_zero, &dummy_one, &dummy_shape, &pn);
        PNorm(&b_centered, &dummy_zero, &dummy_one, &dummy_shape, &pn2);
        r = (dn - dn2) / (pn2 - pn);
        *mn = *mu + *sigma * r;
        *vr = sigma2 * (1.0 +
                          (a_centered * dn - b_centered * dn2) /
                          (pn2 - pn) -
                          r * r);
      }
  }
  
  void RandomNormal(int *n_samples, int *dim_y, double *mu, double *sigma,
                    double *y) {
    /**************************************************
      # Draw a multivariate normal sample.
      # Args:
      #   n_samples: the number of random draws (scalar).
      #   dim_y: the dimension of the random variable (scalar).
      #   mu: the mean (dim_y x 1).
      #   sigma: the covariance (dim_y x dim_y).
      #   y: the random draw (dim_y x 1).
      ******************************************************/
      int i, j, n, this_y;
    double dummy_log_det = 0.0;
    int dim_y_2 = *dim_y * *dim_y;
    int length_z = *n_samples * *dim_y;
    double *sigma_chol  = (double *)Calloc(dim_y_2, double);
    double *z           = (double *)Calloc(length_z, double);
    
    for(n = 0; n <  dim_y_2; n++) {
      sigma_chol[n] = 0.0;
    }
    CholMatrix(dim_y, sigma, sigma_chol, &dummy_log_det);
    GetRNGstate();
    for(n = 0; n < length_z; n++) {
      z[n] = rnorm(0, 1);
    }
    PutRNGstate();
    
    for(n = 0; n < *n_samples; n++) {
      for(i = 0; i < *dim_y; i++) {
        this_y = n + *n_samples * i;
        y[this_y] = mu[i];
        for(j = 0; j < (i + 1); j++) {
          y[this_y] += sigma_chol[i + *dim_y * j] * z[n + *n_samples * j];
        }
      }
    }
    Free(sigma_chol);
    Free(z);
  }
  
  void RandomTruncatedNormal(int *n_samples, double *mu, double *sigma,
                             double *a, double *b, double *y){
    /**************************************************
      # Generates a truncated normal random variable.
      #
      # Args:
      #   n_samples: number of samples.
      #   mu: the mean parameter.
      #   sigma: the scale parameter.
      #   a: the left truncation point.  -99 corresponds to -Inf.
      #   b: the right truncation point.  -99 corresponds to Inf.
      #   y: truncated normal realization of length n_samples.
      #
      ******************************************************/
      int n;
    double u;
    double z = 0.0;
    double tau = 0.0;
    double mn_0 = 0.0;
    double scale_1 = 1.0;
    double shape_0 = 0.0;
    double pnorm_alpha = 0.0;
    double pnorm_beta = 0.0;
    
    GetRNGstate();
    for(n = 0; n < *n_samples; n++){
      /* -99 indicates -Inf for lower limit*/
        if(a[n] == -99.0){
          pnorm_alpha = 0.0;
        } else {
          z = (a[n] - mu[n]) / sigma[n];
          PNorm(&z, &mn_0, &scale_1, &shape_0, &pnorm_alpha);
        }
      /* -99 indicates Inf for upper limit*/
        if(b[n] == -99.0){
          pnorm_beta = 1.0;
        } else {
          z = (b[n] - mu[n]) / sigma[n];
          PNorm(&z, &mn_0, &scale_1, &shape_0, &pnorm_beta);
        }
      u = runif(0, 1);
      tau = pnorm_alpha + u * (pnorm_beta - pnorm_alpha);
      QNorm(&tau, &mn_0, &scale_1, &shape_0, &z);
      y[n] = z * sigma[n] + mu[n];
    }
    PutRNGstate();
  }
  
  void RandomInverseWishart(int *n_row, double *df, double *sai, double *w) {
    /**************************************************
      # Draw an inverse Wishart sample.
      # Args:
      #   n_row: the dimension of the random variable (scalar).
      #   df: degrees of freedom (scalar).
      #   sai: the scale matrix (n_row x n_row).
      #   w: the random draw (n_row x n_row).
      ******************************************************/
      int i, j, s, info;
    int n_row_2 = *n_row * *n_row;
    
    int T_A = 0;
    int T_B = 0;
    
    char uplo = 'L';
    double this_df, dummy_log_det;
    double *sai_inv = (double *)Calloc(*n_row * *n_row, double);
    double *sai_inv_chol = (double *)Calloc(*n_row * *n_row, double);
    double *z_mat = (double *)Calloc(*n_row * *n_row, double);
    
    for (s = 0; s < n_row_2; s++) {
      z_mat[s] = 0;
    }
    /*invert sai*/
      InvertSymmetricMatrix(n_row, sai, sai_inv, &dummy_log_det);
    /*cholesky inverse sai*/
      CholMatrix(n_row, sai_inv, sai_inv_chol, &dummy_log_det);
    GetRNGstate();
    /*construct Z*/
      for (s = 0; s < *n_row; s++) {
        this_df = (*df - s) / 2;
        z_mat[s + *n_row * s] = sqrt(rgamma(this_df, 2));
      }
    if (*n_row > 1) {
      for (i = 1; i < *n_row; i++) {
        for (j = 0; j < i; j++) {
          z_mat[i + *n_row * j] = rnorm(0,1);
        }
      }
    }
    PutRNGstate();
    matrix_multiply (&T_A, &T_B,
                    n_row, n_row, n_row, n_row, sai_inv_chol, z_mat, w);
    F77_CALL(dpotri)(&uplo, n_row, w, n_row, &info);
    for (i = 0; i < *n_row; i++) {
      for (j = i + 1; j < *n_row; j++) {
        w[i + *n_row * j] = w[j + *n_row * i];
      }
    }
    Free(z_mat);
    Free(sai_inv);
    Free(sai_inv_chol);
  }
  
  void RandomDirichlet(int *n_samples, int *n_cells, double *probs,  double *x){
    /**************************************************
      # Generates random Dirichlet samples.
      #
      # Args:
      #   n_samples: number of samples.
      #   n_cells: number of categories.
      #   probs: n_cells x 1 vector of probabilities for each cell.  Should sum to 1.
      #   x: n_cells x 1 vector for storing random Dirichlet realization.
      #
      ******************************************************/
      int i, h;
    double gammasum;
    GetRNGstate();
    for(i = 0; i < *n_samples; i++){
      gammasum = 0.0;
      for(h = 0; h < *n_cells; h++){
        x[i + *n_samples * h] = rgamma(probs[h], 1);
        gammasum += x[i + *n_samples * h];
      }
      for(h = 0; h < *n_cells; h++){
        x[i + *n_samples * h] /= gammasum;
      }
    }
    PutRNGstate();
  }
  
  void RandomMultinomial(int *n_samples, int *n_cells, double *probs, int *x){
    /**************************************************
      # Generates random multinomial samples.
      #
      # Args:
      #   n_samples: number of samples.
      #   n_cells: number of categories.
      #   probs: n_cells x 1 vector of probabilities for each cell.
      #     Should sum to 1.
      #   x: n_cells x 1 vector for storing random multinomial realization.
      #
      ******************************************************/
      GetRNGstate();
    rmultinom(*n_samples, probs, *n_cells, x);
    PutRNGstate();
  }
  
  void RandomUniform(double *a, double *b, double *x){
    /**************************************************
      # Generates a unif(a, b) random variable.
      #
      # Args:
      #   a: lower limit in the support.
      #   b: upper limit in the support.
      #   x: scalar for storing random beta realization.
      #
      ******************************************************/
      //  double double_zero = 0.0;
      //  double double_one = 1.0;
      GetRNGstate();
      *x = runif(*a, *b);
      PutRNGstate();
  }
  
  void RandomBeta(double *shape1, double *shape2, double *x){
    /**************************************************
      # Generates a random beta sample.
      #
      # Args:
      #   shape1: first shape parameter.
      #   shape2: second shape parameter.
      #   x: scalar for storing random beta realization.
      #
      ******************************************************/
      GetRNGstate();
    *x = rbeta(*shape1, *shape2);
    PutRNGstate();
  }
  
  void RandomGamma(double *shape, double *scale, double *x){
    /**************************************************
      # Generates a random gamma sample.
      #
      # Args:
      #   shape: shape parameter.
      #   scale: scale parameter.
      #   x: scalar for storing random gamma realization.
      #
      ******************************************************/
      GetRNGstate();
    *x = rgamma(*shape, *scale);
    PutRNGstate();
  }
  
  void RandomExponential(double *scale, double *x){
    /**************************************************
      # Generates a random exponential sample.
      #
      # Args:
      #   scale: scale parameter.
      #   x: scalar for storing random gamma realization.
      #
      ******************************************************/
      GetRNGstate();
    *x = rexp(*scale);
    PutRNGstate();
  }
  
  void UpdateKroneckerInverseWishart(
    int *nrow_omega,
    int *nrow_sigma,
    double *nu_0, double *sigma_0,
    double *e, double *omega,
    double *sigma) {
    /**************************************************
      # Performs a Gibbs update of the posterior inverse Wishart
      #   random variable sigma, where
      #   e|omega, sigma ~ N(0, Kronecker(omega, sigma))
      #   and sigma ~ IW(nu.0, sigma.0).
      #
      # Args:
      #   nrow_omega: the number of rows of omega.
      #   nrow_sigma: the number of rows of sigma.
      #   nu.0: the prior scale for sigma.
      #   sigma.0: the prior location for sigma.
      #   e: a multivariate normal r.v. that has mean 0
    #     and covariance Kronecker(omega ^ -1, sigma)).
    #   omega: a valid precision matrix.
    #
    # Returns:
    #   sigma: the updated covariance.
    **************************************************/
      int i, j, k, l;
    double df = *nrow_omega + *nu_0;
    int nrow_sigma_2 = *nrow_sigma * *nrow_sigma;
    double *s = (double *) Calloc(nrow_sigma_2, double);
    for(i = 0; i < nrow_sigma_2; i++) {
      s[i] = sigma_0[i];
    }
    for(i = 0; i < *nrow_omega; i++) {
      for(j = 0; j < *nrow_omega; j++) {
        for(k = 0; k < *nrow_sigma; k++) {
          for(l = 0; l < *nrow_sigma; l++) {
            s[k + *nrow_sigma * l] +=
              omega[i + *nrow_omega * j] *
              e[j * *nrow_sigma + k] *
              e[i * *nrow_sigma + l];
          }
        }
      }
    }
    RandomInverseWishart(nrow_sigma, &df, s, sigma);
    Free(s);
  }
  
  /**********************  prior functions *****************************/
    void MakeBasis (int *basis_ind, int *n_bf, double *basis_knots, int *n_tau,
                    double *tau, double *shape, double *basis_matrix) {
      /**********************************************
        # Calculates the basis matrix of quantile functions.
        # Args:
        #  basis_ind: integer indicator where (-1, 0, 1, 2, 3, 4, 5) maps to
        #      ("spline", "Gaussian", "t", "logistic", "ALAP", "Weibull", "gamma").
        #  n_bf: number of basis functions.
        #  basis_knots: vector of knots.
        #  n_tau: number of quantile levels.
        #  tau: vector of quantile levels in (0,1).
        #  shape: optional parameter, needed for t, ALAP, Weibull and gamma distributions.
        #  basis_matrix: n_tau x n_bf array that holds the basis.
        *************************************************/
        int i, j, k, m;
      int knots_length = *n_bf + 4;
      int dummy_bin = 3;
      double scale = 1; /*for identifiability mean is 0 and scale is 1*/
        double mn = 0;
        
        double *q_k  = (double *)Calloc(*n_bf, double); /*quantile function evaluated at the knots*/
          
          void (*q_ptr)(double *, double *, double *, double *, double*); /*pointer to the quantile function*/
          q_ptr = &QNorm;
          switch(*basis_ind) {
            case 0:
              q_ptr = &QNorm;
              break;
              case 1:
                q_ptr = &QT;
                break;
                case 2:
                  q_ptr = &QLogistic;
                  break;
                  case 3:
                    q_ptr = &QAsymmetricLaplace;
                    break;
                    case 4:
                      q_ptr = &QWeibull;
                      break;
                      case 5:
                        q_ptr = &QGamma;
                        break;
          }
          
          if (*basis_ind == -1) {
            for (i = 0; i < *n_tau; i++) {
              basis_matrix[i] = 1.0;
            }
            for (i = 0; i < *n_tau; i++) {
              find_interval(&knots_length, basis_knots, &tau[i], &dummy_bin);
              for (j = 1; j < *n_bf; j++) {
                k = j - 1;
                CubicISpline(&tau[i], &knots_length, basis_knots, &k, &dummy_bin, &basis_matrix[i + *n_tau * j]);
              }
            }
          } else {
            q_k[0] = -999.999;
            if (*basis_ind > 3) {
              q_k[0] = 0;
            }
            q_k[*n_bf - 1] = 999.999;
            for (m = 1; m < (*n_bf - 1); m++) {
              q_ptr(&basis_knots[m], &mn, &scale, shape, &q_k[m]);
            }
            for (j = 0; j < *n_tau; j++) {
              basis_matrix[j] = 1.0;
            }
            /*second column of basis_matrix*/
              for (j = 0; j < *n_tau; j++) {
                if (tau[j] <= basis_knots[1]) {
                  q_ptr(&tau[j], &mn, &scale, shape, &basis_matrix[j + *n_tau * 1]);
                } else{
                  q_ptr(&basis_knots[1], &mn, &scale, shape, &basis_matrix[j + *n_tau * 1]);
                }
              }
            /*rest of columns*/
              if (*n_bf > 2) {
                for (m = 2; m < *n_bf; m++) {
                  for (j = 0; j < *n_tau; j++) {
                    if (tau[j] <= basis_knots[m - 1]) {
                      basis_matrix[j + *n_tau * m] = 0;
                    }else if (tau[j] > basis_knots[m]) {
                      basis_matrix[j + *n_tau * m] = q_k[m] - q_k[m - 1];
                    } else {
                      q_ptr(&tau[j], &mn, &scale, shape, &basis_matrix[j + *n_tau * m]);
                      basis_matrix[j + *n_tau * m] -= q_k[m - 1];
                    }
                  }
                }
              }
          }
          Free(q_k);
    }
  
  void MakeITIGaussian (int *n_bf, double *basis_knots, double *iti) {
    /**********************************************
      # Calculates i'i for all quantile levels, where i is the matrix
      #  of Gaussian basis functions.
      # Args:
      #   n_bf: number of basis functions.
      #   c.indicator: an indicator if the computation should be performed in C.
      # Returns:
      #   iti:  i'i, where i is the Gaussian quantile basis.
      *************************************************/
      
      int i, j, k, l;
    int n_bf_1 = *n_bf - 1;
    //  int n_bf_2 = *n_bf * *n_bf;
    int int_one = 1;
    int int_two = 2;
    int int_three = 3;
    //  int t_a = 0;
    //  int t_b = 0;
    double a = 0.0;
    double b = 0.0;
    double loc = 0.0;
    double scale = 1.0;
    double shape = 1.0;
    double *weight      = (double *) Calloc(n_bf_1, double);
    double *starters    = (double *) Calloc(n_bf_1, double);
    double *finishers   = (double *) Calloc(n_bf_1, double);
    double *bf_means    = (double *) Calloc(n_bf_1, double);
    double *bf_vars     = (double *) Calloc(n_bf_1, double);
    double *tn_means    = (double *) Calloc(n_bf_1, double);
    double *tn_vars     = (double *) Calloc(n_bf_1, double);
    
    if(*n_bf == 2) {
      iti[0] = 1.0;
      iti[1] = 0.0;
      iti[2] = 0.0;
      iti[3] = 1.0;
    } else {
      for(i = 0; i < n_bf_1; i++) {
        weight[i] = basis_knots[i + 1] - basis_knots[i];
        starters[i] = 0.0;
        finishers[i] = 0.0;
        bf_means[i] = 0.0;
        bf_vars[i] = 0.0;
        tn_means[i] = 0.0;
        tn_vars[i] = 0.0;
      }
      
      QNorm(&basis_knots[1], &loc, &scale, &shape, &finishers[0]);
      for(i = 1; i < (*n_bf - 2); i++) {
        QNorm(&basis_knots[i], &loc, &scale, &shape, &starters[i]);
        QNorm(&basis_knots[i + 1], &loc, &scale, &shape, &finishers[i]);
        finishers[i] -= starters[i];
      }
      QNorm(&basis_knots[*n_bf - 2], &loc, &scale, &shape, &starters[*n_bf - 2]);
      
      a = -99.0;
      QNorm(&basis_knots[1], &loc, &scale, &shape, &b);
      TruncatedGaussianMoments(&int_two, &loc, &scale, &a, &b, &tn_means[0], &tn_vars[0]);
      for(i = 1; i < (*n_bf - 2); i++) {
        QNorm(&basis_knots[i], &loc, &scale, &shape, &a);
        QNorm(&basis_knots[i + 1], &loc, &scale, &shape, &b);
        TruncatedGaussianMoments(&int_three, &loc, &scale, &a, &b, &tn_means[i], &tn_vars[i]);
      }
      QNorm(&basis_knots[*n_bf - 2], &loc, &scale, &shape, &a);
      b = 99.0;
      TruncatedGaussianMoments(&int_one, &loc, &scale, &a, &b, &tn_means[*n_bf - 2], &tn_vars[*n_bf - 2]);
      
      bf_means[0] = weight[0] * tn_means[0];
      for(i = 1; i < n_bf_1; i++) {
        bf_means[0] += weight[i] * finishers[0];
      }
      for(i = 1; i < *n_bf - 2; i++) {
        bf_means[i] = weight[i] * (tn_means[i] - starters[i]);
        for(j = i + 1; j < n_bf_1; j++) {
          bf_means[i] += weight[j] * finishers[i];
        }
      }
      bf_means[n_bf_1 - 1] = weight[n_bf_1 - 1] * (tn_means[n_bf_1 - 1] - starters[n_bf_1 - 1]);
      
      for(i = 0; i < n_bf_1; i++) {
        for(j = 0; j < n_bf_1; j++) {
          k = i - j;
          if(k < -1) {
            k = -1;
          }
          if(k > 1) {
            k = 1;
          }
          switch(k) {
            case 1:
              bf_vars[i] +=  weight[j] * (0 - bf_means[i]) * (0 - bf_means[i]);
              break;
              case 0:
                a = tn_means[i] - bf_means[i] - starters[i];
                bf_vars[i] += weight[j] * (tn_vars[i] + a * a);
                break;
                case -1:
                  a = (finishers[i] - bf_means[i]);
                  bf_vars[i] += weight[j] * a * a;
                  break;
          }
        }
      }
      
      for(i = 0; i < n_bf_1; i++) {
        iti[(i + 1) + *n_bf * 0] = bf_means[i];
        iti[0 + *n_bf * (i + 1)] = bf_means[i];
        iti[(i + 1) + *n_bf * (i + 1)] = bf_means[i] * bf_means[i] + bf_vars[i];
      }
      iti[0] = 1.0;
      for(i = 0; i < *n_bf - 2; i++) {
        for(j = (i + 1); j < n_bf_1; j++) {
          l = (i + 1) + *n_bf * (j + 1);
          iti[l] = weight[j] * (tn_means[j] - starters[j]) * finishers[i];
          if(j < *n_bf - 2) {
            for(k = j + 1; k < n_bf_1; k++) {
              iti[l] += weight[k] * finishers[i] * finishers[j];
            }
          }
          iti[(j + 1) + *n_bf * (i + 1)] = iti[l];
        }
      }
      
    }
    Free(weight);
    Free(starters);
    Free(finishers);
    Free(bf_means);
    Free(bf_vars);
    Free(tn_means);
    Free(tn_vars);
  }
  
  void MakeITISpline (int *n_bf, double *basis_knots, double *iti) {
    /**********************************************
      # Calculates i'i for all quantile levels, where i is the matrix
      #  of spline basis functions.
      # Args:
      #   n_bf: number of basis functions.
      #   basis_knots: a sequence of knots for the parametric basis functions.
      # Updates:
      #   iti:  i'i, where i is the Gaussian quantile basis.
      *************************************************/
      int n_bf_2 = *n_bf * *n_bf;
      int bin, k, l, m;
      int n_bf_1 = *n_bf - 1;
      double *spline_coefs = (double *)Calloc((4 * *n_bf  * n_bf_1), double);
      double s1, s2, s3, s4, s5, s6, s7;
      for (m = 0; m < 4 * *n_bf * n_bf_1; m++) {
        spline_coefs[m] = 0.0;
      }
      MakeCubicSplineCoefficients (n_bf, basis_knots, spline_coefs);
      for(m = 0; m < n_bf_2; m++) {
        iti[m] = 0.0;
      }
      for(l = 0; l < *n_bf; l++) {
        for(m = 0 ; m < *n_bf; m++) {
          k = l + *n_bf * m;
          for(bin = 2; bin < n_bf_1; bin++) {
            s7 = (spline_coefs[0 + 4 * l + 4 * *n_bf * bin] * spline_coefs[0 + 4 * m + 4 * *n_bf * bin]);
            s7 *= (pow(basis_knots[bin + 1], 7.0) - pow(basis_knots[bin], 7.0));
            s7 /= 7.0;
            s6 = (spline_coefs[0 + 4 * l + 4 * *n_bf * bin] * spline_coefs[1 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[1 + 4 * l + 4 * *n_bf * bin] * spline_coefs[0 + 4 * m + 4 * *n_bf * bin]);
            s6 *= (pow(basis_knots[bin + 1], 6.0) - pow(basis_knots[bin], 6.0));
            s6 /= 6.0;
            s5 = (spline_coefs[0 + 4 * l + 4 * *n_bf * bin] * spline_coefs[2 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[1 + 4 * l + 4 * *n_bf * bin] * spline_coefs[1 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[2 + 4 * l + 4 * *n_bf * bin] * spline_coefs[0 + 4 * m + 4 * *n_bf * bin]);
            s5 *= (pow(basis_knots[bin + 1], 5.0) - pow(basis_knots[bin], 5.0));
            s5 /= 5.0;
            s4 = (spline_coefs[0 + 4 * l + 4 * *n_bf * bin] * spline_coefs[3 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[1 + 4 * l + 4 * *n_bf * bin] * spline_coefs[2 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[2 + 4 * l + 4 * *n_bf * bin] * spline_coefs[1 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[3 + 4 * l + 4 * *n_bf * bin] * spline_coefs[0 + 4 * m + 4 * *n_bf * bin]);
            s4 *= (pow(basis_knots[bin + 1], 4.0) - pow(basis_knots[bin], 4.0));
            s4 /= 4.0;
            s3 = (spline_coefs[1 + 4 * l + 4 * *n_bf * bin] * spline_coefs[3 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[2 + 4 * l + 4 * *n_bf * bin] * spline_coefs[2 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[3 + 4 * l + 4 * *n_bf * bin] * spline_coefs[1 + 4 * m + 4 * *n_bf * bin]);
            s3 *= (pow(basis_knots[bin + 1], 3.0) - pow(basis_knots[bin], 3.0));
            s3 /= 3.0;
            s2 = (spline_coefs[2 + 4 * l + 4 * *n_bf * bin] * spline_coefs[3 + 4 * m + 4 * *n_bf * bin] +
                    spline_coefs[3 + 4 * l + 4 * *n_bf * bin] * spline_coefs[2 + 4 * m + 4 * *n_bf * bin]);
            s2 *= (pow(basis_knots[bin + 1], 2.0) - pow(basis_knots[bin], 2.0));
            s2 /= 2.0;
            s1 = (spline_coefs[3 + 4 * l + 4 * *n_bf * bin] * spline_coefs[3 + 4 * m + 4 * *n_bf * bin]);
            s1 *= (pow(basis_knots[bin + 1], 1.0) - pow(basis_knots[bin], 1.0));
            iti[k] += s7 + s6 + s5 + s4 + s3 + s2 + s1;
          }
        }
      }
      Free(spline_coefs);
  }
  
  void MakeMuSpline (int *n_bf, double *basis_knots, double *theta_int) {
    /**********************************************
      # Returns the vector for theta that results in a mean 0,
      # unit variance process for the spline basis.
      # Args:
      #   n_bf: number of basis functions.
      #   basis_knots: a sequence of knots for the parametric basis functions.
      # Updates:
      #   theta_int: the weights for the intercept process that result
      #     in a mean 0, unit variance random variable.
      *************************************************/
      int bin, m;
    double a, theta1, theta2;
    int n_bf_1 = *n_bf - 1;
    double *spline_coefs = (double *)Calloc((4 * *n_bf  * n_bf_1), double);
    double *iti = (double *)Calloc((*n_bf  * *n_bf), double);
    double *w = (double *)Calloc((*n_bf), double);
    
    for (m = 0; m < 4 * *n_bf * n_bf_1; m++) {
      spline_coefs[m] = 0.0;
    }
    for (m = 0; m < *n_bf * *n_bf; m++) {
      iti[m] = 0.0;
    }
    for (m = 0; m < *n_bf; m++) {
      w[m] = 0.0;
    }
    
    MakeCubicSplineCoefficients (n_bf, basis_knots, spline_coefs);
    MakeITIGaussian (n_bf, basis_knots, iti);
    w[0] = 1.0;
    for(m = 1; m < *n_bf; m++) {
      for(bin = 2; bin < n_bf_1; bin++) {
        w[m] +=
          (1 / 4) * spline_coefs[0 + 4 * (m + *n_bf * bin)] *
          (pow(basis_knots[bin + 1], 4) - pow(basis_knots[bin], 4)) +
          (1 / 3) * spline_coefs[1 + 4 * (m + *n_bf * bin)] *
          (pow(basis_knots[bin + 1], 3) - pow(basis_knots[bin], 3)) +
          (1 / 2) * spline_coefs[2 + 4 * (m + *n_bf * bin)] *
          (pow(basis_knots[bin + 1], 2) - pow(basis_knots[bin], 2)) +
          (1    ) * spline_coefs[3 + 4 * (m + *n_bf * bin)] *
          (basis_knots[bin + 1] - basis_knots[bin]);
      }
    }
    a = w[1];
    for(m = 2; m < *n_bf; m++) {
      a += w[m];
    }
    w[0] = a;
    for(m = 1; m < *n_bf; m++) {
      w[m] = 1.0;
    }
    theta2 = 0.0;
    QuadraticForm(n_bf, w, iti, &theta2);
    theta2 = sqrt(theta2);
    theta2 = 1 / theta2;
    theta1 = a * theta2;
    theta_int[0] = theta1;
    for(m = 1; m < *n_bf; m++) {
      theta_int[m] = theta2;
    }
    Free(spline_coefs);
    Free(iti);
    Free(w);
  }
  
  void MakePrecTimepoints (int *n_timepoints, double *timepoints,
                           double *prec_timepoints) {
    /**********************************************
      # Creates 2 x 2 precision matrix from timepoints.
      # Args:
      #   n_timepoints: the number of timepoints.
      #   timepoints: a sequence of knots in (0,1).
      # Updates:
      #   prec_timepoints: 2 x 2 inverted precision matrix of intercept
      #     and timepoints.
      *************************************************/
      int t;
    int int_two = 2;
    double dummy_log_det = 0.0;
    double *cov_timepoints = (double *)Calloc(4, double);
    cov_timepoints[0] = *n_timepoints;
    cov_timepoints[1] = timepoints[0];
    cov_timepoints[3] = timepoints[0] * timepoints[0];
    for(t = 1; t < *n_timepoints; t++) {
      cov_timepoints[1] += timepoints[t];
      cov_timepoints[3] += timepoints[t] * timepoints[t];
    }
    cov_timepoints[2] = cov_timepoints[1];
    for(t = 0; t < 4; t++) {
      cov_timepoints[t] /= (double) (*n_timepoints);
    }
    InvertSymmetricMatrix(&int_two, cov_timepoints,
                          prec_timepoints, &dummy_log_det);
    Free(cov_timepoints);
  }
  
  void MakePrecisionThetastar(int *dim_thetastar,
                              double *prec_bf, double *prec_preds, double *prec_time,
                              double *prec_clusters, double *prec_responses,
                              double *prec_thetastar) {
    /**************************************************
      # This function creates the precision matrix of thetastar.
      #
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   prec_bf: The precision of the basis functions.
      #   prec_preds: The precision across predictors.
      #   prec.time: The precision across time.
      #   prec_clusters: The precision across clusters.
      #   prec_responses: The precision across responses.
      #   prec_thetastar: the full precision matrix.
    **************************************************/
      int i;
    
    int n_bf = dim_thetastar[0];
    int n_preds = dim_thetastar[1];
    int dim_time = dim_thetastar[2];
    int n_clusters = dim_thetastar[3];
    int n_responses = dim_thetastar[4];
    
    int n_bf_preds = (n_bf * n_preds);
    int n_bf_preds_2 = n_bf_preds * n_bf_preds;
    
    int multivariate_indicator = dim_time * n_clusters * n_responses;
    
    if(multivariate_indicator == 1) {
      Kronecker(&n_bf, &n_bf, &n_preds, &n_preds, prec_bf,
                prec_preds, prec_thetastar);
      for(i = 0; i < n_bf_preds_2; i++) {
        prec_thetastar[i] *= prec_responses[0];
      }
    } else {
      int length_omega1 =  n_bf_preds_2;
      int length_omega2 = (n_bf_preds * dim_time) * (n_bf_preds * dim_time);
      int length_omega3 = (n_bf_preds * dim_time * n_clusters) *
        (n_bf_preds * dim_time * n_clusters);
      
      int n_bf_preds_time = (n_bf * n_preds * dim_time);
      int n_bf_preds_time_clusters = (n_bf * n_preds * dim_time * n_clusters);
      
      double *omega1 = (double *)Calloc(length_omega1, double);
      double *omega2 = (double *)Calloc(length_omega2, double);
      double *omega3 = (double *)Calloc(length_omega3, double);
      
      for(i = 0; i < length_omega1; i++) {
        omega1[i] = 0.0;
      }
      for(i = 0; i < length_omega2; i++) {
        omega2[i] = 0.0;
      }
      for(i = 0; i < length_omega3; i++) {
        omega3[i] = 0.0;
      }
      Kronecker(&n_bf, &n_bf, &n_preds, &n_preds, prec_bf, prec_preds, omega1);
      Kronecker(&n_bf_preds, &n_bf_preds, &dim_time, &dim_time,
                omega1, prec_time, omega2);
      Kronecker(&n_bf_preds_time, &n_bf_preds_time, &n_clusters, &n_clusters,
                omega2, prec_clusters, omega3);
      Kronecker(&n_bf_preds_time_clusters, &n_bf_preds_time_clusters,
                &n_responses, &n_responses, omega3,
                prec_responses, prec_thetastar);
      Free(omega1);
      Free(omega2);
      Free(omega3);
    }
  }
  
  void LogPriorThetastar(int *dim_thetastar, double *log_dets,
                         double *e_thetastar, double *prec_thetastar, double *lp_thetastar) {
    /**************************************************
      # Returns the log prior of the regression parameters.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: A vector of length 5 containing the log determinants
      #     of prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   e_thetastar: The 5 dimensional array
      #     of the residuals for the regression coefficients.
      #   prec_thetastar: The precision matrix of thetastar.
      #   lp: the log prior probability.
      **************************************************/
      int i;
    int length_thetastar = dim_thetastar[0];
    for(i = 1; i < 5; i++) {
      length_thetastar *= dim_thetastar[i];
    }
    QuadraticForm(&length_thetastar, e_thetastar, prec_thetastar, lp_thetastar);
    *lp_thetastar *= -1.0;
    
    *lp_thetastar += log_dets[0] * (double) (dim_thetastar[1]) *
      (double) (dim_thetastar[2]) * (double) (dim_thetastar[3]) *
      (double) (dim_thetastar[4]);
    *lp_thetastar += log_dets[1] * (double) (dim_thetastar[0]) *
      (double) (dim_thetastar[2]) * (double) (dim_thetastar[3]) *
      (double) (dim_thetastar[4]);
    *lp_thetastar += log_dets[2] * (double) (dim_thetastar[0]) *
      (double) (dim_thetastar[1]) * (double) (dim_thetastar[3]) *
      (double) (dim_thetastar[4]);
    *lp_thetastar += log_dets[3] * (double) (dim_thetastar[0]) *
      (double) (dim_thetastar[1]) * (double) (dim_thetastar[2]) *
      (double) (dim_thetastar[4]);
    *lp_thetastar += log_dets[4] * (double) (dim_thetastar[0]) *
      (double) (dim_thetastar[1]) * (double) (dim_thetastar[2]) *
      (double) (dim_thetastar[3]);
    *lp_thetastar *= 0.5;
  }
  
  void UpdateMu(int *dim_thetastar, double *thetastar, double *e_thetastar,
                double *prec_thetastar, int *indicator_mu, double *mu_0,
                double *prec_mu_0, double *mu, double *mu_array) {
    
    /**************************************************
      # Performs a Gibbs update of mu, updates mu_array and e_thetastar.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   thetastar: a dim_thetastar array of regression coefficients.
      #   e_thetastar: a dim_thetastar array of thetastar - mu_array.
      #   prec_thetastar: the precision matrix of the regression coefficients.
      #   indicator_mu: an indicator if mu
      #     is constant for basis functions greater than 1.
      #   mu_0: the prior mean for mu.
      #   prec_mu_0: the prior precision for mu.
    #   mu: an array of means of dimension num_bf x num_preds x dim_time
    #     or 2 x num_preds x dim.time if indicator_mu = 1.
    #   mu_array: an array of mu of dimension dim_thetastar.
    ******************************************************/
      
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int dim_prod = n_bf * n_preds * dim_time * n_clusters * n_responses;
      
      int int_one = 1;
      int nrow_ata = n_bf * n_preds * dim_time;
      int nrow_v_mat = nrow_ata;
      if(*indicator_mu == 1) {
        nrow_v_mat = 2 * n_preds * dim_time;
      }
      int nrow_v_mat_2 = nrow_v_mat * nrow_v_mat;
      int nrow_ata_2   = nrow_ata * nrow_ata;
      
      double *ata              = (double *)Calloc(nrow_ata_2, double); /* a' prec.thetastar a*/
      double *v_mat            = (double *)Calloc(nrow_v_mat_2, double);
      double *v_mat_inv        = (double *)Calloc(nrow_v_mat_2, double);
      double *v_mat_inv_chol   = (double *)Calloc(nrow_v_mat_2, double);
      double *dummy_v_mat      = (double *)Calloc(nrow_ata_2, double);
      double *v_init           = (double *)Calloc(nrow_ata, double);
      double *v0               = (double *)Calloc(nrow_v_mat, double);
      double *v                = (double *)Calloc(nrow_v_mat, double);
      double *m22              = (double *)Calloc(4, double); /*2 x 2 dummy matrix*/
      
      int d, g, h, i, j, k, l, m, p, row_starter, col_starter;
      
      int gh = n_clusters * n_responses;
      int pd = n_preds * dim_time;
      
      double dummy_log_det = 0.0;
      
      for(i = 0; i < nrow_ata_2; i++) {
      ata[i] = 0.0;
      dummy_v_mat[i] = 0.0;
      }
      for(i = 0; i < nrow_v_mat_2; i++) {
      v_mat[i] = prec_mu_0[i];
      v_mat_inv[i] = 0.0;
      v_mat_inv_chol[i] = 0.0;
      }
      for(i = 0; i < nrow_ata; i++) {
      v_init[i] = 0.0;
      }
      for(i = 0; i < nrow_v_mat; i++) {
      v0[i] = 0.0;
      v[i] = 0.0;
      }
      for(i = 0; i < 4; i++) {
      m22[i] = 0.0;
      }
      
      if(*indicator_mu == 1) {
      for(i = 0; i < gh; i++) {
      row_starter = i * nrow_ata;
      for(j = 0; j < gh; j++) {
      col_starter = j * nrow_ata;
      for(k = 0; k < nrow_ata; k++) {
      for(l = 0; l < nrow_ata; l++) {
      ata[k + nrow_ata * l] +=
      prec_thetastar[(row_starter + k) +  dim_prod * (col_starter + l)];
      }
      }
      }
      }
      for(i = 0; i < pd; i++) {
      row_starter = i * n_bf;
      for(j = 0; j < pd; j++) {
      col_starter = j * n_bf;
      m22[0] = ata[row_starter + nrow_ata * col_starter];
      m22[1] = 0;
      m22[2] = 0;
      m22[3] = 0;
      for(k = 1; k < n_bf; k++) {
      m22[1] += ata[(row_starter + k) + nrow_ata * col_starter];
      m22[2] += ata[row_starter + nrow_ata * (col_starter + k)];
      for(l = 1; l < n_bf; l++) {
      m22[3] += ata[(row_starter + k) + nrow_ata * (col_starter + l)];
      }
      }
      v_mat[i * 2 + 0 + nrow_v_mat * (j * 2 + 0)] += m22[0];
      v_mat[i * 2 + 1 + nrow_v_mat * (j * 2 + 0)] += m22[1];
      v_mat[i * 2 + 0 + nrow_v_mat * (j * 2 + 1)] += m22[2];
      v_mat[i * 2 + 1 + nrow_v_mat * (j * 2 + 1)] += m22[3];
      }
      }
      
      for(l = 0; l < gh; l++) {
      col_starter = l * nrow_ata;
      for(m = 0; m < nrow_ata_2; m++) {
      dummy_v_mat[m] = 0.0;
      }
      for(k = 0; k < gh; k++) {
      row_starter = k * nrow_ata;
      for(i = 0; i < nrow_ata; i++) {
      for(j = 0; j < nrow_ata; j++) {
      dummy_v_mat[i + nrow_ata * j] +=
      prec_thetastar[row_starter + i + dim_prod * (col_starter + j)];
      }
      }
      }
      for(i = 0; i < nrow_ata; i++) {
      for(j = 0; j < nrow_ata; j++) {
      v_init[i] +=
      dummy_v_mat[i + nrow_ata * j] * thetastar[col_starter + j];
      }
      }
      }
      for(i = 0; i < nrow_v_mat; i++) {
      for(j = 0; j < nrow_v_mat; j++) {
      v0[i] += prec_mu_0[i + nrow_v_mat * j] * mu_0[j];
      }
      }
      for(i = 0; i < pd; i++) {
      v0[2 * i]      += v_init[n_bf * i];
      for(m = 1; m < n_bf; m++) {
      v0[2 * i + 1]  += v_init[n_bf * i + m];
      }
      }
      } else {
      for(i = 0; i < gh; i++) {
      row_starter = i * nrow_v_mat;
      for(j = 0; j < gh; j++) {
      col_starter = j * nrow_v_mat;
      for(k = 0; k < nrow_v_mat; k++) {
      for(l = 0; l < nrow_v_mat; l++) {
      v_mat[k + nrow_v_mat * l] +=
      prec_thetastar[(row_starter + k) + dim_prod * (col_starter + l)];
      }
      }
      }
      }
      
      for(i = 0; i < nrow_v_mat; i++) {
      for(j = 0; j < nrow_v_mat; j++) {
      v0[i] += prec_mu_0[i + nrow_v_mat * j] * mu_0[j];
      }
      }
      
      for(j = 0; j < gh; j++) {
      col_starter = j * nrow_v_mat;
      for(m = 0; m < nrow_v_mat_2; m++) {
      dummy_v_mat[m] = 0.0;
      }
      for(i = 0; i < gh; i++) {
      row_starter = i * nrow_v_mat;
      for(k = 0; k < nrow_v_mat; k++) {
      for(l = 0; l < nrow_v_mat; l++) {
      dummy_v_mat[k + nrow_v_mat * l] +=
      prec_thetastar[(row_starter + k) + dim_prod * (col_starter + l)];
      }
      }
      }
      for(k = 0; k < nrow_v_mat; k++) {
      for(l = 0; l < nrow_v_mat; l++) {
      v0[k] +=
      dummy_v_mat[k + nrow_v_mat * l] * thetastar[col_starter + l];
      }
      }
      }
      }
      InvertSymmetricMatrix(&nrow_v_mat, v_mat, v_mat_inv, &dummy_log_det);
      for(i = 0; i < nrow_v_mat; i++) {
      for(j = 0; j < nrow_v_mat; j++) {
      v[i] += v_mat_inv[i + nrow_v_mat * j] * v0[j];
      }
      }
      dummy_log_det = 0.0;
      for(i = 0; i < nrow_v_mat; i++) {
      dummy_log_det += v[i];
      }
      dummy_log_det = 0.0;
      for(i = 0; i < nrow_v_mat * nrow_v_mat; i++) {
      dummy_log_det += v_mat[i];
      }
      RandomNormal(&int_one, &nrow_v_mat, v, v_mat_inv, mu);
      if(*indicator_mu == 1) {
      for(m = 0; m < n_bf; m++) {
      if(m == 0) {
      i = 0;
      } else {
      i = 1;
      }
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = i + 2 * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      } else {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = m + n_bf * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      }
      Free(ata);
      Free(v_mat);
      Free(v_mat_inv);
      Free(dummy_v_mat);
      Free(v_init);
      Free(v0);
      Free(v);
      Free(m22);
  }
      
      void UpdateCovResponses(int *dim_thetastar, double *log_dets,
      double *e_thetastar,
      double *prec_bf, double *prec_preds,
      double *prec_time, double *prec_clusters, double *prec_responses,
      double *cov_responses,
      double *nu_0, double *sigma_0,
      double *prec_thetastar, double *lp_thetastar) {
      
      /**************************************************
      # Gibbs updates cov_responses, the covariance of the responses, as well
      #  as prec_responses, prec_thetastar, log_dets, and lp_thetastar.
      #
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: a vector of length 5 containing the log determinants of
      #      prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   e_thetastar: dim_thetastar array of residuals of thetastar.
      #   prec_bf: The precision of the basis functions.
      #   prec_preds: The precision across predictors.
      #   prec_time: The precision across time.
      #   prec_clusters: The precision across clusters.
      #   prec_responses: The precision across responses.
      #   cov_responses: The covariance across responses.
      #   nu_0: the prior scale for cov_responses.
      #   sigma_0: the prior location for cov_responses.
      #   prec_thetastar: the updated precision.
      #   lp_thetastar: the log prior of thetastar.
      **************************************************/
      int nrow_omega1 = dim_thetastar[0] * dim_thetastar[1];
      int nrow_omega2 = dim_thetastar[2] * dim_thetastar[3];
      int nrow_omega = nrow_omega1 * nrow_omega2;
      
      int length_omega1 = nrow_omega1 * nrow_omega1;
      int length_omega2 = nrow_omega2 * nrow_omega2;
      int length_omega = nrow_omega * nrow_omega;
      
      int i;
      
      double *omega1 = (double *)Calloc(length_omega1, double);
      double *omega2 = (double *)Calloc(length_omega2, double);
      double *omega = (double *)Calloc(length_omega, double);
      
      for(i = 0; i < length_omega1; i++) {
      omega1[i] = 0.0;
      }
      for(i = 0; i < length_omega2; i++) {
      omega2[i] = 0.0;
      }
      for(i = 0; i < length_omega; i++) {
      omega[i] = 0.0;
      }
      Kronecker(&dim_thetastar[0], &dim_thetastar[0], &dim_thetastar[1],
      &dim_thetastar[1], prec_bf, prec_preds, omega1);
      Kronecker(&dim_thetastar[2], &dim_thetastar[2], &dim_thetastar[3],
      &dim_thetastar[3], prec_time, prec_clusters, omega2);
      
      double log_det = 0.0;
      for(i = 0; i < length_omega1; i++) {
      log_det += omega1[i];
      }
      log_det = 0.0;
      for(i = 0; i < length_omega2; i++) {
      log_det += omega2[i];
      }
      
      Kronecker(&nrow_omega1, &nrow_omega1, &nrow_omega2, &nrow_omega2, omega1,
      omega2, omega);
      UpdateKroneckerInverseWishart(&nrow_omega, &dim_thetastar[4], nu_0, sigma_0,
      e_thetastar, omega, cov_responses);
      InvertSymmetricMatrix(&dim_thetastar[4], cov_responses, prec_responses,
      &log_dets[4]);
      MakePrecisionThetastar(dim_thetastar, prec_bf, prec_preds, prec_time,
      prec_clusters, prec_responses, prec_thetastar);
      LogPriorThetastar(dim_thetastar, log_dets, e_thetastar,
      prec_thetastar, lp_thetastar);
      
      Free(omega1);
      Free(omega2);
      Free(omega);
      }
      
      void LogPriorXi(int *dim_thetastar, double *log_dets,
      double *e_xi, double *prec_clusters, double *tau2_xi, double *lp_xi) {
      /**************************************************
      # Returns the log prior of the tail parameters.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: A vector of length 5 containing the log determinants.
      #     of prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   e_xi: The n_clusters x n_responses array.
      #     of the residuals for xi.
      #   prec_clusters: The n_clusters x n_clusters correlation matrix.
      #   tau2_xi: The precision of xi.
      #   lp_xi: the log prior probability.
      **************************************************/
      int i;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_responses_2 = n_responses * n_responses;
      int length_xi = dim_thetastar[3] * dim_thetastar[4];
      int length_xi_2 = length_xi * length_xi;
      /*tau2_xi * kronecker(prec_clusters, diag(n_responses))*/
      double *prec_xi = (double *)Calloc(length_xi_2, double);
      double *diag_responses = (double *)Calloc(n_responses_2, double);
      
      for(i = 0; i < length_xi_2; i++) {
      prec_xi[i] = 0.0;
      }
      for(i = 0; i < n_responses_2; i++) {
      diag_responses[i] = 0.0;
      }
      for(i = 0; i < n_responses; i++) {
      diag_responses[i + n_responses * i] = 1.0;
      }
      Kronecker(&n_clusters, &n_clusters, &n_responses, &n_responses,
      prec_clusters, diag_responses, prec_xi);
      for(i = 0; i < length_xi_2; i++) {
      prec_xi[i] *= *tau2_xi;
      }
      QuadraticForm(&length_xi, e_xi, prec_xi, lp_xi);
      *lp_xi *= -1.0;
      *lp_xi += log_dets[3] * (double) (n_responses) + log(*tau2_xi) *
      (double) (n_clusters) * (double) (n_responses);
      *lp_xi *= 0.5;
      }
      
      void UpdatePriorXi (int *dim_thetastar, double *log_dets,
      double *xi, double *mu_xi, double *e_xi, double *mu_0_xi,
      double *prec_clusters, double *kappa_xi,
      double *tau2_xi,
      double *a_xi, double *b_xi, double *lp_xi) {
      /**********************************************
      # Gibbs updates mu_xi and tau2_xi, and updates e_xi and lp_xi.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: A vector of length 5 containing the log determinants
      #     of prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   xi: an n_clusters x n_responses array of tail parameters.
      #   mu_xi: an n_responses x 1 array of tail means.
      #   mu_0_xi: an n_responses x 1 array of prior means for mu_xi.
      #   prec_clusters: precision matrix of the clusters.
      #   kappa_xi: prior sample size for the mean.
      #   a_xi: prior shape for the precision.
      #   b_xi: prior scale for the precision.
      # Updates:
      #   mu_xi:  updated mean.
      #   e_xi:  updated residuals.
      #   tau2_xi:  updated precision.
      #   lp: updated prior.
      *************************************************/
      int i, j, k;
      int int_one = 1;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_responses_2 = n_responses * n_responses;
      int dim_xi = n_clusters * n_responses;
      int dim_lambda = dim_xi * dim_xi;
      
      double shape_xi = 0.5 * (double) (n_clusters) *
      (double) (n_responses) + *a_xi;
      double scale_1 = 0.0;
      double scale_2 = 0.0;
      double scale_3 = 0.0;
      double sum_prec_clusters = 0.0;
      double scale_xi, scaling_factor;
      
      
      double *dummy_e = (double *)Calloc(dim_xi, double);
      double *lambda = (double *)Calloc(dim_lambda, double);
      double *diag_n_responses = (double *)Calloc(n_responses_2, double);
      double *colsums_clusters = (double *)Calloc(n_clusters, double);
      double *weighted_sum = (double *)Calloc(n_responses, double);
      double *little_omega = (double *)Calloc(n_responses, double);
      double *big_omega = (double *)Calloc(n_responses_2, double);
      double *final_mu = (double *)Calloc(n_responses, double);
      double *final_sigma = (double *)Calloc(n_responses_2, double);
      
      
      for(i = 0; i < dim_xi; i++) {
      dummy_e[i] = 1.0;
      }
      for(i = 0; i < dim_lambda; i++) {
      lambda[i] = 0.0;
      }
      for(i = 0; i < n_responses_2; i++) {
      big_omega[i] = 0.0;
      diag_n_responses[i] = 0.0;
      final_sigma[i] = 0.0;
      }
      for(i = 0; i < n_responses; i++) {
      little_omega[i] = 0.0;
      diag_n_responses[i + n_responses * i] = 1.0;
      final_mu[i] = 0.0;
      }
      for(j = 0; j < n_clusters; j++) {
      colsums_clusters[j] = 0.0;
      for(i = 0; i < n_clusters; i++) {
      colsums_clusters[j] += prec_clusters[i + n_clusters * j];
      sum_prec_clusters += prec_clusters[i + n_clusters * j];
      }
      }
      scaling_factor = sum_prec_clusters + *kappa_xi;
      for(i = 0; i < n_responses; i++) {
      weighted_sum[i] = *kappa_xi * mu_0_xi[i];
      }
      
      for(j = 0; j < n_clusters; j++) {
      for(i = 0; i < n_responses; i++) {
      weighted_sum[i] += colsums_clusters[j] * xi[i + n_responses * j];
      }
      }
      
      Kronecker(&n_clusters, &n_clusters, &n_responses, &n_responses, prec_clusters,
      diag_n_responses, lambda);
      QuadraticForm(&dim_xi, xi, lambda, &scale_1);
      for(i = 0; i < n_responses;  i++) {
      scale_2 += mu_0_xi[i] * mu_0_xi[i];
      }
      scale_2 *= *kappa_xi;
      for(i = 0; i < n_responses; i++) {
      scale_3 += weighted_sum[i] * weighted_sum[i];
      }
      scale_3 *= -1.0 / scaling_factor;
      scale_xi = 1.0 / *b_xi + 0.5 * (scale_1 + scale_2 + scale_3);
      scale_xi = 1.0 / scale_xi;
      RandomGamma(&shape_xi, &scale_xi, tau2_xi);
      for(i = 0; i < n_responses; i++) {
      weighted_sum[i] /= scaling_factor;
      diag_n_responses[i + n_responses * i] /= (scaling_factor * *tau2_xi);
      }
      RandomNormal(&int_one, &n_responses, weighted_sum, diag_n_responses, mu_xi);
      
      for(i = 0; i < n_clusters; i++) {
      for(j = 0; j < n_responses; j++) {
      k = i * n_responses + j;
      e_xi[k] = xi[k] - mu_xi[j];
      }
      }
      LogPriorXi(dim_thetastar, log_dets, e_xi, prec_clusters, tau2_xi, lp_xi);
      
      Free(dummy_e);
      Free(lambda);
      Free(diag_n_responses);
      Free(colsums_clusters);
      Free(weighted_sum);
      }
      
      void UpdateRho(int *dim_thetastar, double *log_dets, double *e_thetastar,
      double *lp_thetastar, double *prec_bf, double *prec_preds,
      double *prec_time, double *prec_clusters, double *prec_responses,
      double *cov_clusters, double *prec_thetastar,
      double *tau_low, double *tau_high, double *lp_xi_low, double *lp_xi_high,
      double *tau2_xi_low, double *tau2_xi_high,
      double *e_xi_low, double *e_xi_high,
      double *rho, double *max_rho, double *sd_rho, int *acc_rho,
      double *shape1, double *shape2, int *correlation_type,
      double *distance_matrix) {
      /**********************************************
      # Metropolis updates rho.
      #
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: a vector of length 5 containing the log determinants of
      #      prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   e_thetastar: dim_thetastar array of residuals of thetastar.
      #   lp_thetastar: a scalar containing the current log prior thetastar value.
      #   prec_bf: the precision of the basis functions.
      #   prec_preds: the precision across predictors.
      #   prec_time: the precision across time.
      #   prec_clusters: the precision across clusters.
      #   prec_responses: the precision across responses.
      #   cov_clusters: the covariance across clusters
      #    (not checked to be the inverse of the precision).
      #   prec_thetastar: the precision matrix of thetastar
      #    (not checked against all the composite precisions).
      #   tau_low: the threshold below which a GEV tail is fit.
      #   tau_high: the threshold above which a GEV tail is fit.
      #   lp_xi_low: the log prior of the xi.low parameter.
      #   lp_xi_high: the log prior of the xi.high parameter.
      #   tau2_xi_low: the precision of the xi.low parameter.
      #   tau2_xi_high: the precision of the xi.high parameter.
      #   e_xi_low: the errors of the xi.low parameter.
      #   e_xi_high: the errors of the xi.high parameter.
      #   rho: the autocorrelation/range parameter.
      #   max_rho: the maximum value for rho.
      #   sd_rho: the candidate standard deviation.
      #   acc_rho: counter for the number accepted.
      #   shape1: first prior shape parameter.
      #   shape2: second prior shape parameter.
      #   correlation_type: an integer of either 0/1/2 indicating if the
      #     cluster correlation structure is independent, time series or spatial.
      #   distance_matrix: a valid distance matrix of dimension
      #    n_clusters x n_clusters.
      # Updates rho, acc_rho, prec_clusters, cov_clusters, prec_thetastar,
      #         log_dets, lp_thetastar, lp_xi_low, lp_xi_high.
      **********************************************/
      int i;
      int int_one = 1;
      int n_clusters = dim_thetastar[3];
      int length_prec_clusters = n_clusters * n_clusters;
      int length_prec_thetastar = dim_thetastar[0];
      for(i = 1; i < 5; i++) {
      length_prec_thetastar *= dim_thetastar[i];
      }
      length_prec_thetastar = length_prec_thetastar * length_prec_thetastar;
      double rho_1, rho_prime, can_rho_1, can_rho_prime, can_rho, ff, can_ff,
      can_lp_xi_low, can_lp_xi_high, rho_post, can_rho_post;
      double can_lp_thetastar = 0.0;
      double lp_rho = 0.0;
      double can_lp_rho = 0.0;
      double double_zero = 0.0;
      double double_one = 1.0;
      double z = 0.0;
      double e = 0.0;
      double dummy_log_det = 0.0;
      
      double *can_log_dets = (double *)Calloc(5, double);
      double *can_prec_clusters = (double *)Calloc(length_prec_clusters, double);
      double *can_cov_clusters = (double *)Calloc(length_prec_clusters, double);
      double *can_prec_thetastar = (double *)Calloc(length_prec_thetastar, double);
      //
      for(i = 0; i < 5; i++) {
      can_log_dets[i] = log_dets[i];
      }
      
      for(i = 0; i < length_prec_clusters; i++) {
      can_prec_clusters[i] = 0.0;
      can_cov_clusters[i] = 0.0;
      }
      for(i = 0; i < length_prec_thetastar; i++) {
      can_prec_thetastar[i] = 0.0;
      }
      
      rho_1 = *rho / *max_rho; /*in unit interval*/
      rho_prime = log(rho_1) - log(1.0 - rho_1);
      RandomNormal(&int_one, &int_one, &double_zero, &double_one, &z);
      can_rho_prime = rho_prime + *sd_rho * z;
      can_rho_1 = exp(can_rho_prime) / (1.0 + exp(can_rho_prime));
      can_rho = can_rho_1 * *max_rho;
      
      //  Rprintf("before make precision.\n");
      //  Rprintf("*rho rho_1, z, sd_rho, can_rho_prime, can_rho_1, can_rho, *max_rho = %f %f %f %f %f %f %f %f\n", *rho, rho_1, z, *sd_rho, can_rho_prime, can_rho_1, can_rho, *max_rho);
      
      if(*correlation_type == 1) {
      MakeAutoregressivePrecision(&n_clusters, &can_rho, distance_matrix,
      can_prec_clusters, &can_log_dets[3]);
      } else {
      MakeSpatialPrecision(&n_clusters, &can_rho, distance_matrix,
      can_prec_clusters, &can_log_dets[3]);
      }
      MakePrecisionThetastar(dim_thetastar, prec_bf, prec_preds, prec_time,
      can_prec_clusters, prec_responses, can_prec_thetastar);
      LogPriorThetastar(dim_thetastar, can_log_dets, e_thetastar,
      can_prec_thetastar, &can_lp_thetastar);
      LogDBeta(&rho_1, shape1, shape2, &lp_rho);
      LogDBeta(&can_rho_1, shape1, shape2, &can_lp_rho);
      
      ff = log(rho_1) - log(1 - rho_1);
      can_ff = log(can_rho_1) - log(1 - can_rho_1);
      
      if(*tau_low == 0.0) {
      *lp_xi_low = 0.0;
      }
      if(*tau_high == 1.0) {
      *lp_xi_high = 0.0;
      }
      can_lp_xi_low = 0.0;
      if(*tau_low > 0) {
      LogPriorXi(dim_thetastar, can_log_dets, e_xi_low, can_prec_clusters,
      tau2_xi_low, &can_lp_xi_low);
      }
      can_lp_xi_high = 0.0;
      if(*tau_high < 1.0) {
      LogPriorXi(dim_thetastar, can_log_dets, e_xi_high, can_prec_clusters,
      tau2_xi_high, &can_lp_xi_high);
      }
      rho_post = *lp_thetastar + *lp_xi_low + *lp_xi_high + ff + lp_rho;
      can_rho_post = can_lp_thetastar + can_lp_xi_low +
      can_lp_xi_high + can_ff + can_lp_rho;
      RandomExponential(&double_one, &e);
      
      if(e > (rho_post  - can_rho_post)) {
      *rho = can_rho;
      set_equal(&length_prec_clusters, can_prec_clusters, prec_clusters);
      InvertSymmetricMatrix(&n_clusters, prec_clusters,
      cov_clusters, &dummy_log_det);
      set_equal(&length_prec_thetastar, can_prec_thetastar, prec_thetastar);
      log_dets[3] = can_log_dets[3];
      *lp_thetastar = can_lp_thetastar;
      *lp_xi_low = can_lp_xi_low;
      *lp_xi_high = can_lp_xi_high;
      *acc_rho += 1;
      }
      Free(can_log_dets);
      Free(can_prec_clusters);
      Free(can_cov_clusters);
      Free(can_prec_thetastar);
      }
      
      void UpdatePriorThetastar (int *dim_thetastar, double *log_dets,
      double *thetastar, double *e_thetastar, double *lp_thetastar,
      double *prec_bf, double *prec_preds, double *prec_time,
      double *prec_clusters, double *prec_responses,
      double *cov_clusters, double *cov_responses,
      double *prec_thetastar, int *indicator_mu, double *mu, double *mu_array,
      double *mu_0, double *prec_mu_0, double *nu_0, double *sigma_0,
      double *tau_low, double *tau_high, double *lp_xi_low, double *lp_xi_high,
      double *tau2_xi_low, double *tau2_xi_high, double *e_xi_low,
      double *e_xi_high, double *rho, double *max_rho, double *sd_rho,
      int *acc_rho, double *shape1, double *shape2,
      int *correlation_type, double *distance_matrix) {
      /************************************************
      # updates the prior for thetastar
      # Args:
      #   dim_thetastar: a vector of length 6 containing the
      #    dimension of thetastar
      #    and the number of timepoints.
      #   log_dets: a vector of length 5 containing the log determinants of
      #      prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   thetastar: dim_thetastar array of regression parameters.
      #   e_thetastar: dim_thetastar array of residuals of thetastar.
      #   lp_thetastar: a scalar containing the current log prior
      #     thetastar value.
      #     This value is assumed correct as an input.
      #   prec_bf: The precision of the basis functions.
      #   prec_preds: The precision across predictors.
      #   prec_time: The precision across time.
      #   prec_clusters: The precision across clusters.
      #   prec_responses: the precision across responses.
      #   cov_clusters: The covariance across clusters.
      #   cov_responses: the covariance across responses.
      #   prec_thetastar: the precision matrix of thetastar
      #     (not checked against all the composite precisions).
      #   indicator_mu: an indicator if mu
      #     is constant for basis functions greater than 1.
      #   mu: an array of means of dimension n_bf x num_preds x dim.time
      #     or 2 x num_preds x dim_time if indicator_mu = 1.
      #   mu_array: an array of mu of dimension dim_thetastar.
      #   mu_0: the prior mean for mu.
      #   prec_mu_0: the prior precision for mu.
      #   nu_0: the prior scale for cov_responses.
      #   sigma_0: the prior location for cov_responses.
      #   tau_low: the threshold below which a GEV tail is fit.
      #   tau_high: the threshold above which a GEV tail is fit.
      #   lp_xi_low: the log prior of the xi_low parameter.
      #   lp_xi_high: the log prior of the xi_high parameter.
      #   tau2_xi_low: the precision of the xi_low parameter.
      #   tau2_xi_high: the precision of the xi_high parameter.
      #   e_xi_low: the errors of the xi_low parameter.
      #   e_xi_high: the errors of the xi_high parameter.
      #   rho: the autocorrelation/range parameter.
      #   max_rho: the maximum value for rho.
      #   sd_rho: the candidate standard deviation.
      #   acc_rho: counter for the number accepted.
      #   shape1: first prior shape parameter.
      #   shape2: second prior shape parameter.
      #   correlation_type: an integer of either 0/1/2 indicating if the
      #     cluster correlation structure is
      #     independent, time series or spatial.
      #   distance_matrix: a valid distance matrix of dimension n_clusters x
      #     n_clusters.
      # Updates:
      #   mu, mu_array, e_thetastar, cov_responses, prec_responses,
      #   prec_thetastar, lp_thetastar, log_dets,
      #   rho, prec_clusters, cov_clusters, prec_thetastar,
      #   lp_xi_low, lp_xi_high.
      ************************************************/
      
      UpdateMu(dim_thetastar, thetastar, e_thetastar, prec_thetastar,
      indicator_mu, mu_0, prec_mu_0,
      mu, mu_array);
      UpdateCovResponses(dim_thetastar, log_dets, e_thetastar,
      prec_bf, prec_preds, prec_time, prec_clusters, prec_responses,
      cov_responses,
      nu_0, sigma_0,
      prec_thetastar, lp_thetastar);
      if(*correlation_type > 0) {
      UpdateRho(dim_thetastar, log_dets, e_thetastar, lp_thetastar,
      prec_bf, prec_preds, prec_time, prec_clusters, prec_responses,
      cov_clusters, prec_thetastar,
      tau_low, tau_high, lp_xi_low, lp_xi_high,
      tau2_xi_low, tau2_xi_high, e_xi_low, e_xi_high,
      rho, max_rho, sd_rho, acc_rho,
      shape1, shape2, correlation_type, distance_matrix);
      }
      }
      
      void StoreParameters (int *i, int *burn, int *num_keep, int *thin,
      int *dim_thetastar,
      double *pred_mns, double *pred_sds,
      double *response_mns, double *response_sds,
      int *indicator_mu, double *mu,
      double *cov_responses, double *rho, double *thetastar,
      double *mu_keep, double *cov_responses_keep,
      double *rho_keep, double *thetastar_keep,
      double *xi_low, double *xi_high,
      double *xi_low_keep, double *xi_high_keep,
      double *mu_xi_low, double *mu_xi_high,
      double *mu_xi_low_keep, double *mu_xi_high_keep,
      double *tau2_xi_low, double *tau2_xi_high,
      double *tau2_xi_low_keep, double *tau2_xi_high_keep) {
      /************************************************
      # stores parameters.
      # Args:
      #   i: current iteration.
      #   burn: number of samples to burn.
      #   num_keep: number of samples to keep.
      #   thin: number of samples to thin.
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   pred_mns: n_preds x 1 array of standardizing locations of the predictors.
      #   pred_sds: n_preds x 1 array of standardizing scales of the predictors.
      #   response_mns: n_responses x 1 array of standardizing locations
      #     of the responses.
      #   response_sds: n_responses x 1 array of standardizing scales
      #     of the responses.
      #   indicator_mu: an indicator if mu
      #     is constant for basis functions greater than 1.
      #   mu: an array of means of dimension n_bf x num_preds x dim_time
      #     or 2 x num_preds x dim_time if indicator_mu = 1.
      #   cov_responses: the covariance across responses.
      #   rho: the autocorrelation/range parameter.
      #   thetastar: dim_thetastar array of regression parameters.
      #   mu_keep: an array for storing mu of dimension
      #     num_keep x n_bf x num_preds x dim_time
      #     or num_keep x 2 x num_preds x dim_time if indicator_mu = 1.
      #   cov_responses_keep: an array for storing cov_responses of dimension
      #     num_keep x n_responses x n_responses.
      #   rho_keep: an array for storing rho of dimension
      #     num_keep x 1.
      #   thetastar_keep: array for storing thetastar of dimension
      #     num_keep x dim_thetastar.
      #   xi_low: array of dimension n_clusters x n_responses
      #     of the lower tail parameter.
      #   xi_high: array of dimension n_clusters x n_responses
      #     of the upper tail parameter.
      #   xi_low_keep: array for storing xi_low of dimension
      #     num_keep x n_clusters x n_responses.
      #   xi_high_keep: array for storing xi_high of dimension
      #     num_keep x n_clusters x n_responses.
      #   mu_xi_low: n_responses x 1 mean of the lower tail parameters.
      #   mu_xi_high: n_responses x 1 mean of the upper tail parameters.
      #   mu_xi_low_keep: array for storing mu_xi_low of dimension
      #     num_keep x n_responses.
      #   mu_xi_high_keep: array for storing mu_xi_high of dimension
      #     num_keep x n_responses.
      #   tau2_xi_low: scalar precision of the lower tail parameters.
      #   tau2_xi_high: scalar precision of the upper tail parameters.
      #   tau2_xi_low_keep: array for storing tau2_xi_low of dimension num_keep x 1.
      #   tau2_xi_high_keep: array for storing tau2_xi_high of
      #     dimension num_keep x 1.
      # Updates:
      #    mu_keep, cov_responses_keep, rho_keep, thetastar_keep,
      #    xi_low_keep, xi_high_keep, mu_xi_low_keep, mu_xi_high_keep,
      #    tau2_xi_low_keep, tau2_xi_high_keep.
      ************************************************/
      int d, g, h, k, k1, k2, m, p;
      int i2 = (*i + 1 - *burn) / *thin - 1;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_bf_preds = n_bf * n_preds;
      int n_responses_2 = n_responses * n_responses;
      
      int length_thetastar = dim_thetastar[0];
      for(k = 1; k < 5; k++) {
      length_thetastar *= dim_thetastar[k];
      }
      
      int dim_mu = n_bf_preds * dim_time;
      if(*indicator_mu == 1) {
      dim_mu = 2 * n_preds * dim_time;
      }
      int dim_xi = n_clusters * n_responses;
      
      double *this_thetastar = (double *)Calloc(length_thetastar, double);
      
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      for(d = 0; d < dim_time; d++) {
      /*m = 0, p = 0 case*/
      k1 = 0 + n_bf * 0;
      k2 = d + dim_time * (g + n_clusters * h);
      k = k1 + n_bf_preds * k2;
      this_thetastar[k] = pred_sds[0] * thetastar[k] +
      response_mns[h] / response_sds[h];
      for(p = 0; p < n_preds; p++) {
      this_thetastar[k] += thetastar[(0 + n_bf * p) +
      n_bf_preds * k2] * pred_mns[p];
      }
      this_thetastar[k] *= response_sds[h];
      /*m > 0, p = 0 case*/
      for(m = 1; m < n_bf; m++) {
      k1 = m + n_bf * 0;
      k = k1 + n_bf_preds * k2;
      this_thetastar[k] = thetastar[k] * pred_sds[0];
      for(p = 0; p < n_preds; p++) {
      this_thetastar[k] += thetastar[m + n_bf * p +
      n_bf_preds * k2] * pred_mns[p];
      }
      this_thetastar[k] *= response_sds[h];
      }
      /*p > 0 case*/
      for(p = 1; p < n_preds; p++) {
      for(m = 0; m < n_bf; m++) {
      k1 = m + n_bf * p;
      k = k1 + n_bf_preds * k2;
      this_thetastar[k] = thetastar[k] * pred_sds[p] * response_sds[h];
      }
      }
      }
      }
      }
      
      for(k = 0; k < dim_mu; k++) {
      mu_keep[i2 + *num_keep * k] = mu[k];
      }
      for(k = 0; k < n_responses_2; k++) {
      cov_responses_keep[i2 + *num_keep * k] = cov_responses[k];
      }
      rho_keep[i2] = *rho;
      for(k = 0; k < length_thetastar; k++) {
      thetastar_keep[i2 + *num_keep * k] = this_thetastar[k];
      }
      for(k = 0; k < dim_xi; k++) {
      xi_low_keep[i2 + *num_keep * k] = xi_low[k];
      xi_high_keep[i2 + *num_keep * k] = xi_high[k];
      }
      for(k = 0; k < n_responses; k++) {
      mu_xi_low_keep[i2 + *num_keep * k] = mu_xi_low[k];
      mu_xi_high_keep[i2 + *num_keep * k] = mu_xi_high[k];
      }
      tau2_xi_low_keep[i2] = *tau2_xi_low;
      tau2_xi_high_keep[i2] = *tau2_xi_high;
      
      Free(this_thetastar);
      }
      
      void CheckPrior (int *mcmc_ints, double *mcmc_doubles,
      double *thetastar,
      double *prec_bf, double *prec_preds, double *prec_time,
      double *mu_0, double *prec_mu_0,
      double *nu_0, double *sigma_0,
      double *distance_matrix,
      double *xi_low, double *xi_high,
      double *mu_keep, double *cov_responses_keep, double *rho_keep,
      double *thetastar_keep,
      double *xi_low_keep, double *xi_high_keep,
      double *mu_xi_low_keep, double *mu_xi_high_keep,
      double *tau2_xi_low_keep, double *tau2_xi_high_keep
      ) {
      /************************************************
      # updates the prior for thetastar
      # Args:
      #   mcmc_ints: a vector of length 10 containing dim_thetastar,
      #     how many samples to burn, keep and thin, indicator_mu,
      #     correlation_type, and cop.
      #   mcmc_doubles: a vector where:
      #     1) the first n_preds elements are the locations of the
      #        untransformed predictors.
      #     2) the next n_preds elements are the scales of the
      #        untransformed predictors.
      #     3) the next n_responses elements are the locations of the
      #        untransformed responses.
      #     4) the next n_responses elements are the scales of the
      #        untransformed responses.
      #     5) the next 3 values are max_rho, shape1, shape2.
      #     6) the next 8 values are tau_low, tau_high, a_xi_low, a_xi_high,
      #        b_xi_low, b_xi_high, kappa_low, kappa_high.
      #     7) the next n_responses values are mu_0_xi_low.
      #     8) the next n_responses values are mu_0_xi_high.
      #   thetastar: dim_thetastar array of regression parameters.
      #   prec_bf: The precision of the basis functions.
      #   prec_preds: The precision across predictors.
      #   prec_time: The precision across time.
      #   mu_0: the prior mean for mu.
      #   prec_mu_0: the prior precision for mu.
      #   nu_0: the prior scale for cov_responses.
      #   sigma_0: the prior location for cov_responses.
      #   distance_matrix: a valid distance matrix of dimension
      #     n_clusters x n_clusters.
      #   xi_low: n_clusters x n_responses array of lower tail parameters.
      #   xi_high: n_clusters x n_responses array of upper tail parameters.
      #   mu_keep: an array for storing mu of dimension
      #     num_keep x n_bf x num_preds x dim_time
      #     or num_keep x 2 x num_preds x dim_time if indicator_mu = 1.
      #   cov_responses_keep: an array for storing cov_responses of dimension
      #     num_keep x n_responses x n_responses.
      #   rho_keep: an array for storing rho of dimension
      #     num_keep x 1.
      #   thetastar_keep: array for storing thetastar of dimension
      #     num_keep x dim_thetastar.
      #   xi_low_keep: array for storing xi_low of dimension
      #     num_keep x n_clusters x n_responses.
      #   xi_high_keep: array for storing xi_high of dimension
      #     num_keep x n_clusters x n_responses.
      #   mu_xi_low_keep: array for storing mu_xi_low of dimension
      #     num_keep x n_responses.
      #   mu_xi_high_keep: array for storing mu_xi_high of dimension
      #     num_keep x n_responses.
      #   tau2_xi_low_keep: array for storing tau2_xi_low of dimension num_keep x 1.
      #   tau2_xi_high_keep: array for storing tau2_xi_high of dimension
      #     num_keep x 1.
      # Updates:
      #    mu_keep, cov_responses_keep, rho_keep, thetastar_keep,
      #    xi_low_keep, xi_high_keep, mu_xi_low_keep, mu_xi_high_keep,
      #    tau2_xi_low_keep, tau2_xi_high_keep.
      ************************************************/
      int d, g, h, i, j, k, l, m, p;
      int *dim_thetastar = (int *)Calloc(6, int);
      for(i = 0; i < 6; i++) {
      dim_thetastar[i] = mcmc_ints[i];
      }
      int length_thetastar = 1;
      for(i = 0; i < 5; i++) {
      length_thetastar *= dim_thetastar[i];
      }
      
      int length_thetastar_2 = length_thetastar * length_thetastar;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_clusters_2 = n_clusters * n_clusters;
      int n_responses_2 = n_responses * n_responses;
      
      int max_dim = n_bf;
      if(max_dim < n_preds) {
      max_dim = n_preds;
      }
      if(max_dim < dim_time) {
      max_dim = dim_time;
      }
      int max_dim_2 = max_dim * max_dim;
      
      int burn = mcmc_ints[6];
      int keep = mcmc_ints[7];
      int thin = mcmc_ints[8];
      int iters = burn + keep;
      int num_keep =  keep / thin;
      
      int indicator_mu = mcmc_ints[9];
      int length_mu = n_bf * n_preds * dim_time;
      if(indicator_mu == 1) {
      length_mu = 2 * n_preds * dim_time;
      }
      
      int correlation_type = mcmc_ints[10];
      int acc_rho = 0;
      k = 2 * (n_preds + n_responses);
      
      double max_rho = mcmc_doubles[k];
      double shape1 = mcmc_doubles[k + 1];
      double shape2 = mcmc_doubles[k + 2];
      double tau_low = mcmc_doubles[k + 3];
      double tau_high = mcmc_doubles[k + 4];
      double a_xi_low = mcmc_doubles[k + 5];
      double a_xi_high = mcmc_doubles[k + 6];
      double b_xi_low = mcmc_doubles[k + 7];
      double b_xi_high = mcmc_doubles[k + 8];
      double kappa_xi_low = mcmc_doubles[k + 9];
      double kappa_xi_high = mcmc_doubles[k + 10];
      
      double rho = 0.5 * max_rho;
      double lp_thetastar = 0.0;
      double sd_rho = 1.0;
      double tau2_xi_low = 1.0;
      double tau2_xi_high = 1.0;
      double lp_xi_low = 0.0;
      double lp_xi_high = 0.0;
      
      double *pred_mns = (double *)Calloc(n_preds, double);
      double *pred_sds = (double *)Calloc(n_preds, double);
      double *response_mns = (double *)Calloc(n_responses, double);
      double *response_sds = (double *)Calloc(n_responses, double);
      double *mu = (double *)Calloc(length_mu, double);
      double *mu_array = (double *)Calloc(length_thetastar, double);
      double *e_thetastar = (double *)Calloc(length_thetastar, double);
      double *cov_clusters = (double *)Calloc(n_clusters_2, double);
      double *prec_clusters = (double *)Calloc(n_clusters_2, double);
      double *cov_responses = (double *)Calloc(n_clusters_2, double);
      double *prec_responses = (double *)Calloc(n_clusters_2, double);
      double *dummy_matrix = (double *)Calloc(max_dim_2, double);
      double *prec_thetastar = (double *)Calloc(length_thetastar_2, double);
      double *mu_xi_low = (double *)Calloc(n_responses, double);
      double *mu_xi_high = (double *)Calloc(n_responses, double);
      double *mu_0_xi_low = (double *)Calloc(n_responses, double);
      double *mu_0_xi_high = (double *)Calloc(n_responses, double);
      double *e_xi_low = (double *)Calloc(n_clusters * n_responses, double);
      double *e_xi_high = (double *)Calloc(n_clusters * n_responses, double);
      double *log_dets = (double *)Calloc(5, double);
      
      for(i = 0; i < n_preds; i++) {
      pred_mns[i] = mcmc_doubles[i];
      pred_sds[i] = mcmc_doubles[i + n_preds];
      }
      k = 2 * n_preds;
      l = 2 * (n_preds + n_responses) + 11;
      for(i = 0; i < n_responses; i++) {
      response_mns[i] = mcmc_doubles[k + i];
      response_sds[i] = mcmc_doubles[k + n_responses + i];
      mu_0_xi_low[i] = mcmc_doubles[l + i];
      mu_0_xi_high[i] = mcmc_doubles[l + n_responses + i];
      }
      for(i = 0; i < length_mu; i++) {
      mu[i] = mu_0[i];
      }
      if(indicator_mu == 1) {
      for(m = 0; m < n_bf; m++) {
      if(m == 0) {
      i = 0;
      } else {
      i = 1;
      }
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = i + 2 * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      } else {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = m + n_bf * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      }
      for(i = 0; i < n_responses_2; i++) {
      cov_responses[i] = 0.0;
      prec_responses[i] = 0.0;
      }
      for(i = 0; i < n_responses; i++) {
      k = i + n_responses * i;
      cov_responses[k] = 1.0;
      prec_responses[k] = 1.0;
      }
      
      for(i = 0; i < 5; i++) {
      log_dets[i] = 0.0;
      }
      for(i = 0; i < max_dim_2; i++) {
      dummy_matrix[i] = 0.0;
      }
      for(i = 0; i < n_clusters_2; i++) {
      cov_clusters[i] = 0.0;
      prec_clusters[i] = 0.0;
      }
      for(i = 0; i < n_clusters; i++) {
      k = i + n_clusters * i;
      cov_clusters[k] = 1.0;
      prec_clusters[k] = 1.0;
      }
      if(correlation_type == 1) {
      MakeAutoregressivePrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      if(correlation_type == 2) {
      MakeSpatialPrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      InvertSymmetricMatrix(&n_bf, prec_bf, dummy_matrix, &log_dets[0]);
      InvertSymmetricMatrix(&n_preds, prec_preds, dummy_matrix, &log_dets[1]);
      InvertSymmetricMatrix(&dim_time, prec_time, dummy_matrix, &log_dets[2]);
      InvertSymmetricMatrix(&n_clusters, prec_clusters, cov_clusters, &log_dets[3]);
      InvertSymmetricMatrix(&n_responses, prec_responses, cov_responses,
      &log_dets[4]);
      
      for(i = 0; i < 5; i++) {
      log_dets[i] *= -1.0;
      }
      for(i = 0; i < length_thetastar_2; i++) {
      prec_thetastar[i] = 0.0;
      }
      MakePrecisionThetastar(dim_thetastar,
      prec_bf, prec_preds, prec_time, prec_clusters, prec_responses,
      prec_thetastar);
      
      LogPriorThetastar(dim_thetastar, log_dets, e_thetastar, prec_thetastar,
      &lp_thetastar);
      
      for(h = 0; h < n_responses; h++) {
      mu_xi_low[h] = mu_0_xi_low[h];
      mu_xi_high[h] = mu_0_xi_high[h];
      }
      for(i = 0; i < n_clusters; i++) {
      for(j = 0; j < n_responses; j++) {
      k = i * n_responses + j;
      e_xi_low[k] = xi_low[k] - mu_xi_low[j];
      e_xi_high[k] = xi_high[k] - mu_xi_high[j];
      }
      }
      
      if(tau_low > 0.0) {
      LogPriorXi(dim_thetastar, log_dets, e_xi_low, prec_clusters, &tau2_xi_low,
      &lp_xi_low);
      }
      if(tau_high < 1.0) {
      LogPriorXi(dim_thetastar, log_dets, e_xi_high, prec_clusters, &tau2_xi_high,
      &lp_xi_high);
      }
      for(i = 0; i < iters; i++) {
      UpdatePriorThetastar (dim_thetastar, log_dets, thetastar, e_thetastar,
      &lp_thetastar, prec_bf, prec_preds, prec_time, prec_clusters,
      prec_responses,
      cov_clusters, cov_responses,
      prec_thetastar,
      &indicator_mu, mu, mu_array,
      mu_0, prec_mu_0,
      nu_0, sigma_0,
      &tau_low, &tau_high, &lp_xi_low, &lp_xi_high,
      &tau2_xi_low, &tau2_xi_high, e_xi_low, e_xi_high,
      &rho, &max_rho, &sd_rho, &acc_rho,
      &shape1, &shape2, &correlation_type, distance_matrix);
      
      if(tau_low > 0.0) {
      UpdatePriorXi (dim_thetastar, log_dets,
      xi_low, mu_xi_low, e_xi_low, mu_0_xi_low,
      prec_clusters, &kappa_xi_low,
      &tau2_xi_low,
      &a_xi_low, &b_xi_low, &lp_xi_low);
      }
      if(tau_high < 1.0) {
      UpdatePriorXi (dim_thetastar, log_dets,
      xi_high, mu_xi_high, e_xi_high, mu_0_xi_high,
      prec_clusters, &kappa_xi_high,
      &tau2_xi_high,
      &a_xi_high, &b_xi_high, &lp_xi_high);
      }
      
      if(i >= burn) {
      if((i + 1 - burn) % thin == 0) {
      StoreParameters (&i, &burn, &num_keep, &thin,
      dim_thetastar,
      pred_mns, pred_sds,
      response_mns, response_sds,
      &indicator_mu,
      mu, cov_responses, &rho, thetastar,
      mu_keep, cov_responses_keep, rho_keep, thetastar_keep,
      xi_low, xi_high,
      xi_low_keep, xi_high_keep,
      mu_xi_low, mu_xi_high,
      mu_xi_low_keep, mu_xi_high_keep,
      &tau2_xi_low, &tau2_xi_high,
      tau2_xi_low_keep, tau2_xi_high_keep);
      }
      }
      }
      Free(pred_mns);
      Free(pred_sds);
      Free(response_mns);
      Free(response_sds);
      Free(mu);
      Free(mu_array);
      Free(e_thetastar);
      Free(cov_clusters);
      Free(prec_clusters);
      Free(cov_responses);
      Free(prec_responses);
      Free(dummy_matrix);
      Free(prec_thetastar);
      Free(mu_xi_low);
      Free(mu_xi_high);
      Free(mu_0_xi_low);
      Free(mu_0_xi_high);
      Free(e_xi_low);
      Free(e_xi_high);
      Free(log_dets);
      }
      
      ///*******data functions ****/
      void ThetastarToTheta(int *dim_thetastar, double *thetastar, double *theta,
      double *timepoints, int *this_cluster, int *this_response) {
      /*************************************************************
      # projects thetastar linearly over time onto theta, and thresholds
      #   if necessary.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   thetastar: dim_thetastar[1:5] array of regression parameters.
      #   theta: array of dimension
      #     n_bf x n_preds x n_timepoints x n_clusters x n_responses
      #     of thresholded regression parameters.
      #     Not checked to come from thetastar.
      #   timepoints: vector of length n_timepoints sorted in increasing order.
      #   this_cluster: cluster index that is updated for thetastar.
      #   this_response: response index that is updated for thetastar.
      # Returns:
      #   theta: n_bf x n_preds x n_timepoints x n_clusters x n_responses array
      #    of thresholded regression coefficients.
      *************************************************************/
      int j, k, m, p;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_timepoints = dim_thetastar[5];
      int this_level_d = n_bf * n_preds * dim_time *
      (*this_cluster + n_clusters * *this_response);
      int this_level_j = n_bf * n_preds * n_timepoints *
      (*this_cluster + n_clusters * *this_response);
      double neg_part;
      
      if(dim_time == 1) {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(j = 0; j < n_timepoints; j++) {
      theta[m + n_bf * (p + n_preds * j) + this_level_j] =
      thetastar[m + n_bf * (p + n_preds * 0) + this_level_d];
      }
      }
      }
      } else {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(j = 0; j < n_timepoints; j++) {
      theta[m + n_bf * (p + n_preds * j) + this_level_j] =
      thetastar[m + n_bf * (p + n_preds * 0) + this_level_d] +
      timepoints[j] *
      thetastar[m + n_bf * (p + n_preds * 1) + this_level_d];
      }
      }
      }
      }
      for(j = 0; j < n_timepoints; j++) {
      for(m = 1; m < n_bf; m++) {
      k = m + n_bf * (0 + n_preds * j) + this_level_j;
      if(theta[k] < 0.0) {
      theta[k] = 0.001;
      }
      neg_part = theta[k];
      for(p = 1; p < n_preds; p++) {
      k = m + n_bf * (p + n_preds * j) + this_level_j;
      if(theta[k] < 0.0) {
      neg_part += theta[k];
      } else {
      neg_part -= theta[k];
      }
      }
      if(neg_part < 0.0) {
      for(p = 1; p < n_preds; p++) {
      theta[m + n_bf * (p + n_preds * j) + this_level_j] = 0.0;
      }
      }
      }
      }
      }
      
      void GenerateResponse(int *dim_thetastar, int *n_y, int *basis_ind,
      double *basis_knots, double *spline_coefs, double *y, int *status,
      double *u, double *x, double *theta, double *shape,
      double *tau_low, double *tau_high, int *pareto_low, int *pareto_high,
      double *xi_low, double *xi_high) {
      /************************************************
      # generates the response.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #  basis_ind: integer indicator where (-1, 0, 1, 2, 3, 4, 5) maps to
      #      ("spline", "Gaussian", "t", "logistic", "ALAP", "Weibull", "gamma").
      #  basis_knots: vector of knots.
      #  y: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     to be updated.
      #  status: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of status indicators where:
      #     -99 indicates that the value should be ignored;
      #     0 indicates that the value is missing and treated as continuous;
      #     1 indicates that the value is observed;
      #     2 indicates that the value is right-censored;
      #     3 indicates that the value is left-censored;
      #     4 indicates that the value is interval-censored.
      #  u: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #    of uniform random variables.
      #  x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters of
      #    predictors.
      #  theta: n_bf x n_preds x n_timepoints x n_clusters x n_repsonses array
      #    of thresholded regression parameters.
      #    theta is not checked to be in the space of valid regression parameters.
      #  shape: n_responses length vector of shape parameters for bases.
      #  tau_low: the quantile level below which a parametric tail is fit.
      #  tau_high: the quantile level above which a parametric tail is fit.
      #  xi_low: array of dimension n.clusters x n.responses of
      #   lower tail parameters.
      #  xi_high: array of dimension n.clusters x n.responses of
      #   upper tail parameters.
      # Updates y.
      ************************************************/
      int c_theta, c_x, c_y, g, h, i, j, k, m, p;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_bp = n_bf * n_preds;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int basis_matrix_length = max_n_y * n_bf;
      int this_x_length = max_n_y * n_preds;
      int int_zero = 0;
      int int_one = 1;
      int spline_knots_length = n_bf + 4;
      int bin = 0;
      
      double mn, s, z;
      double thresh = 0.0;
      double scale = 0.0;
      double double_zero = 0.0;
      double double_one = 1.0;
      /*pointer to the log density function*/
      void (*dens_ptr)(double *, double *, double *, double *, double*);
      /*pointer to the quantile function*/
      void (*q_ptr)(double *, double *, double *, double *, double*);
      dens_ptr = &LogDNorm;
      q_ptr = &QNorm;
      
      double *this_u = (double *)Calloc(max_n_y, double);
      double *this_theta = (double *)Calloc(n_bf * n_preds, double);
      double *basis_matrix = (double *)Calloc(basis_matrix_length, double);
      double *w = (double *)Calloc(basis_matrix_length, double);
      /*quantile function evaluated at the knots*/
      double *this_x = (double *)Calloc(this_x_length, double);
      double *q_knots = (double *)Calloc(n_bf, double);
      double *m_low = (double *)Calloc(n_bf, double);
      double *m_high = (double *)Calloc(n_bf, double);
      double *i_low = (double *)Calloc(n_bf, double);
      double *i_high = (double *)Calloc(n_bf, double);
      double *m_low_p = (double *)Calloc(n_preds, double);
      double *m_high_p = (double *)Calloc(n_preds, double);
      double *i_low_p = (double *)Calloc(n_preds, double);
      double *i_high_p = (double *)Calloc(n_preds, double);
      double *mn_low_p = (double *)Calloc(n_preds, double);
      double *mn_high_p = (double *)Calloc(n_preds, double);
      double *sd_low_p = (double *)Calloc(n_preds, double);
      double *sd_high_p = (double *)Calloc(n_preds, double);
      
      for(i = 0; i < max_n_y; i++) {
      this_u[i] = 0.0;
      }
      for(i = 0; i < n_bp; i++) {
      this_theta[i] = 0.0;
      }
      
      for(i = 0; i < basis_matrix_length; i++) {
      basis_matrix[i] = 0.0;
      w[i] = 0.0;
      }
      for(i = 0; i < this_x_length; i++) {
      this_x[i] = 0.0;
      }
      for(m = 0; m < (n_bf - 1); m++) {
      q_knots[m] = -999.0;
      }
      q_knots[n_bf - 1] = 999.0;
      
      for(p = 0; p < n_preds; p++) {
      m_low_p[p] = 0.0;
      m_high_p[p] = 0.0;
      i_low_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      }
      
      for(m = 0; m < n_bf; m++) {
      m_low[m] = 0.0;
      m_high[m] = 0.0;
      i_low[m] = 0.0;
      i_high[m] = 0.0;
      }
      i_low[0] = 1.0;
      i_high[0] = 1.0;
      
      if(basis_ind[0] == -1) {
      MakeBasis(&basis_ind[0], &n_bf, basis_knots,
      &int_one, tau_low, &shape[0], i_low);
      MakeBasis(&basis_ind[0], &n_bf, basis_knots,
      &int_one, tau_high, &shape[0], i_high);
      for(m = 0; m < (n_bf - 1); m++) {
      CubicMSpline2(tau_low, &spline_knots_length, basis_knots,
      spline_coefs, &m, &bin, &m_low[m + 1]);
      CubicMSpline2(tau_high, &spline_knots_length, basis_knots,
      spline_coefs, &m, &bin, &m_high[m + 1]);
      m_low[m + 1] *= *tau_low;
      m_high[m + 1] *= (1 - *tau_high);
      }
      }
      
      for(j = 0; j < n_timepoints; j++) {
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      switch(basis_ind[h]) {
      case 0:
      dens_ptr = &LogDNorm;
      q_ptr = &QNorm;
      break;
      case 1:
      dens_ptr = &LogDT;
      q_ptr = &QT;
      break;
      case 2:
      dens_ptr = &LogDLogistic;
      q_ptr = &QLogistic;
      break;
      case 3:dens_ptr = &LogDAsymmetricLaplace;
      q_ptr = &QAsymmetricLaplace;
      break;
      case 4:
      dens_ptr = &LogDWeibull;
      q_ptr = &QWeibull;
      break;
      case 5:
      dens_ptr = &LogDGamma;
      q_ptr = &QGamma;
      break;
      }
      if(basis_ind[0] != -1)  {
      if(n_bf > 2) {
      for(m = 1; m < n_bf - 1; m++) {
      q_ptr(&basis_knots[m], &double_zero, &double_one, &shape[h],
      &q_knots[m]);
      }
      }
      }
      
      c_theta = n_bf * n_preds * (j + n_timepoints * (g + n_clusters * h));
      c_x = max_n_y * n_preds * (j + n_timepoints * g);
      c_y = max_n_y * (j + n_timepoints * (g + n_clusters * h));
      
      for(i = 0; i < max_n_y; i++) {
      this_u[i] = u[i + c_y];
      for(p = 0; p < n_preds; p++) {
      this_x[i + max_n_y * p] = x[i + max_n_y * p + c_x];
      }
      }
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      k = m + n_bf * p;
      this_theta[k] = theta[k + c_theta];
      }
      }
      if(*tau_low > 0.0) {
      if(basis_ind[h] == -1) {
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, m_low, m_low_p);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, i_low, i_low_p);
      } else {
      MakeBasis(&basis_ind[h], &n_bf, basis_knots, &int_one, tau_low,
      &shape[h], i_low);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, i_low, i_low_p);
      for(p = 0; p < n_preds; p++) {
      mn_low_p[p] = this_theta[0 + n_bf * p];
      sd_low_p[p] = this_theta[1 + n_bf * p];
      }
      }
      }
      if(*tau_high < 1.0) {
      if(basis_ind[h] == -1) {
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, m_high, m_high_p);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, i_high, i_high_p);
      } else {
      MakeBasis(&basis_ind[h], &n_bf, basis_knots, &int_one, tau_high,
      &shape[h], i_high);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf,
      &int_one, this_theta, i_high, i_high_p);
      for(p = 0; p < n_preds; p++) {
      mn_high_p[p] = this_theta[0 + n_bf * p];
      sd_high_p[p] = this_theta[(n_bf - 1) + n_bf * p];
      }
      if(n_bf > 2) {
      for(m = 1; m < (n_bf - 1); m++) {
      for(p = 0; p < n_preds; p++) {
      mn_high_p[p] += (this_theta[(m + 1) + n_bf * p] -
      this_theta[m + n_bf * p]) * q_knots[m];
      }
      }
      }
      }
      }
      matrix_multiply(&int_zero, &int_one, &max_n_y, &n_preds, &n_bf,
      &n_preds, this_x, this_theta, w);
      MakeBasis(&basis_ind[h], &n_bf, basis_knots, &max_n_y, this_u,
      &shape[h], basis_matrix);
      for(i = 0; i < n_y[g]; i++) {
      if(status[i + c_y] == 0 || status[i + c_y] > 1) {
      y[i + c_y] = 0.0;
      for(m = 0; m < n_bf; m++) {
      k = i + max_n_y * m;
      y[i + c_y] += basis_matrix[k] * w[k];
      }
      }
      if(this_u[i] < *tau_low) {
      if(basis_ind[h] == -1) {
      thresh = 0.0;
      scale = 0.0;
      for(p = 0; p < n_preds; p++) {
      k = i + max_n_y * p;
      thresh += this_x[k] * i_low_p[p];
      scale += this_x[k] * m_low_p[p];
      }
      if(*pareto_low == 0) {
      y[i + c_y] = thresh + scale * (log(this_u[i]) - log(*tau_low));
      } else {
      y[i + c_y] = thresh -
      (scale / exp(xi_low[g + n_clusters * h])) *
      (pow((this_u[i] / *tau_low),
      -exp(xi_low[g + n_clusters * h])) - 1.0);
      }
      } else {
      thresh = 0.0;
      mn = 0.0;
      s = 0.0;
      for(p = 0; p < n_preds; p++) {
      k = i + max_n_y * p;
      thresh += this_x[k] * i_low_p[p];
      mn += this_x[k] * mn_low_p[p];
      s += this_x[k] * sd_low_p[p];
      }
      dens_ptr(&thresh, &mn, &s, &shape[h], &scale);
      scale = *tau_low * exp(-1.0 * scale);
      if(*pareto_low == 0) {
      y[i + c_y] = thresh + scale * (log(this_u[i]) - log(*tau_low));
      } else {
      y[i + c_y] = thresh -
      (scale / exp(xi_low[g + n_clusters * h])) *
      (pow((this_u[i] / *tau_low),
      -exp(xi_low[g + n_clusters * h])) - 1.0);
      }
      }
      }
      if(this_u[i] > *tau_high) {
      if(basis_ind[h] == -1) {
      thresh = 0.0;
      scale = 0.0;
      for(p = 0; p < n_preds; p++) {
      k = i + max_n_y * p;
      thresh += this_x[k] * i_high_p[p];
      scale += this_x[k] * m_high_p[p];
      }
      if(*pareto_high == 0) {
      y[i + c_y] = thresh - scale * (log(1.0 - this_u[i]) -
      log(1.0 - *tau_high));
      } else {
      y[i + c_y] = thresh +
      (scale / exp(xi_high[g + n_clusters * h])) *
      (pow((  (1.0 - this_u[i])  / (1.0 - *tau_high)),
      -exp(xi_high[g + n_clusters * h])) - 1.0);
      }
      } else {
      thresh = 0.0;
      mn = 0.0;
      s = 0.0;
      for(p = 0; p < n_preds; p++) {
      k = i + max_n_y * p;
      thresh += this_x[k] * i_high_p[p];
      mn += this_x[k] * mn_high_p[p];
      s += this_x[k] * sd_high_p[p];
      }
      z = (thresh - mn) / s;
      dens_ptr(&z, &double_zero, &double_one, &shape[h], &scale);
      scale -= log(s);
      scale = (1.0 - *tau_high) * exp(-1.0 * scale);
      if(*pareto_high == 0) {
      y[i + c_y] = thresh -
      scale * (log(1.0 - this_u[i]) - log(1.0 - *tau_high));
      } else {
      y[i + c_y] = thresh +
      (scale / exp(xi_low[g + n_clusters * h])) *
      (pow(((1.0 - this_u[i]) / (1.0 - *tau_high)),
      -exp(xi_high[g + n_clusters * h])) - 1.0);
      }
      }
      }
      }
      }
      }
      }
      
      Free(this_u);
      Free(this_theta);
      Free(basis_matrix);
      Free(w);
      Free(this_x);
      Free(q_knots);
      Free(m_low);
      Free(m_high);
      Free(i_low);
      Free(i_high);
      Free(m_low_p);
      Free(m_high_p);
      Free(i_low_p);
      Free(i_high_p);
      Free(mn_low_p);
      Free(mn_high_p);
      Free(sd_low_p);
      Free(sd_high_p);
      }
      
      void FindCubicRoot(double *a, double *b, double *c, double *d,
      double *y, double *kappa1, double *kappa2, double *u) {
      /*******************************************************
      # find the solution to y = a * u ^ 3 + b * u ^ 2 + c * u + d,
      #  where u in (kappa1, kappa2).
      #  Does not check if the solution exists.
      # Args:
      #   a: cubic polynomial coefficient.
      #   b: quadratic polynomial coefficient.
      #   c: linear polynomial coefficient.
      #   d: constant polynomial coefficient.
      #   y: scalar solution.
      #   kappa1: lower interval limit.
      #   kappa2: upper interval limit.
      #   u: solution.
      # Updates u, where u in (kappa1, kappa2) and
      #   y = a * u ^ 3 + b * u ^ 2 + c * u + d.
      *******************************************************/
      int iter = 0;
      double eps = pow(10, -5);
      double f = *a * pow(*u, 3.0) + *b * pow(*u, 2.0) + *c * *u + *d;
      double f_prime = 3.0 * *a * pow(*u, 2) + 2.0 * *b * *u + *c;
      double u1 = *kappa1;
      double u2 = *kappa2;
      
      while(fabs(*y - f) > eps && iter < 1000) {
      iter ++;
      *u += (*y - f) / f_prime;
      /*if Newton-Raphson fails, use bisection.*/
      if(*u < *kappa1 || *u > *kappa2) {
      *u = 0.5 * (u1 + u2);
      while(fabs(*y - f) > eps && iter < 1000) {
      iter ++;
      f = *a * pow(*u, 3.0) + *b * pow(*u, 2.0) + *c * *u + *d;
      if(f < *y) {
      u1 = *u;
      } else {
      u2 = *u;
      }
      *u = 0.5 * (u1 + u2);
      }
      }
      f = *a * pow(*u, 3.0) + *b * pow(*u, 2.0) + *c * *u + *d;
      f_prime = 3.0 * *a * pow(*u, 2) + 2.0 * *b * *u + *c;
      }
      if(iter == 1000) {
      *u = 0.5 * (*kappa1 + *kappa2);
      Rprintf("max_iter!");
      }
      }
      
      void LowerParetoDensity(double *y, double *tau_low, double *thresh_low,
      double *scale_low, double *xi_low, int *pareto_low,
      double *u, double *log_like) {
      /****************************************************
      # returns the log density and quantile level for the lower tail.
      # Args:
      #   y: response that is below the lower threshold.
      #   tau_low: quantile level below which the lower threshold is fit.
      #   thresh_low: lower threshold.
      #   scale_low: lower scale parameter.
      #   xi_low: lower shape parameter.
      #   pareto_low: indicator if pareto tail is used.
      # Returns:
      #   u: the quantile level.
      #   log_like: the log likelihood.
      ****************************************************/
      double z = (*thresh_low - *y) / *scale_low;
      if(*pareto_low == 0) {
      *u = *tau_low * exp(-z);
      *log_like = log(*tau_low) - log(*scale_low) - z;
      } else {
      *u = *tau_low * (pow(1 +  exp(*xi_low) * z, -1.0 / exp(*xi_low)));
      *log_like = log(*tau_low) - log(*scale_low) -
      (1.0 / exp(*xi_low) + 1.0) * log(1.0 + exp(*xi_low) * z);
      }
      if(*u <= 0) {
      *u = 0.0000001;
      }
      }
      
      void UpperParetoDensity(double *y, double *tau_high, double *thresh_high,
      double *scale_high, double *xi_high, int *pareto_high,
      double *u, double *log_like) {
      /****************************************************
      # returns the log density and quantile level for the upper tail.
      # Args:
      #   y: response that is above the upper threshold.
      #   tau_high: quantile level above which the upper threshold is fit.
      #   thresh_high: upper threshold.
      #   scale_high: upper scale parameter.
      #   xi_high: upper shape parameter.
      #   pareto_high: indicator if pareto tail is used.
      # Updates:
      #   u: the quantile level.
      #   log_like: the log likelihood.
      ****************************************************/
      double z = (*y - *thresh_high) / *scale_high;
      if(*pareto_high == 0) {
      *u = *tau_high + (1.0 - *tau_high) * (1.0 - exp(-z));
      *log_like = log(1.0 - *tau_high) - log(*scale_high) - z;
      } else {
      *u = *tau_high + (1.0 - *tau_high) *
      (1.0 - pow(1.0 +  exp(*xi_high) * z, -1.0 / exp(*xi_high)));
      *log_like = log(1 - *tau_high) - log(*scale_high) -
      (1.0 / exp(*xi_high) + 1) * log(1.0 + exp(*xi_high) * z);
      }
      if(*u >= 1.0) {
      *u = 1.0 - 0.0000001;
      }
      }
      
      void MiddleSplineDensity(int *n_bf, int *n_preds, double *y, double *u,
      double *log_like, int *bin, double *x, double *theta,
      double *break_matrix, double *basis_knots, double *spline_coefs) {
      /************************************************************
      # returns the log density and quantile level
      #   for the middle of the distribution.
      # Args:
      #   n_bf: number of basis functions.
      #   n_preds: number of predictors.
      #   y: response.
      #   u: uniform random variable.
      #   log_like: log likelihood.
      #   bin: cell y is located in.
      #   x: n_preds length vector of covariates.
      #   theta: n_bf x n_preds matrix of regression coefficients.
      #   break_matrix: n_bf x n_preds matrix of the basis functions times
      #     the regression coefficients.
      #   basis_knots: vector of spline knots.
      #   spline_coefs: 4 x n_bf x (n_bf - 1) array of cubic spline coefficients.
      # Updates u, log_like, bin.
      ************************************************************/
      int k, m, p;
      int length_breaks = *n_bf - 2;
      double a, b, c, d;
      
      double *w = (double *)Calloc(*n_bf, double);
      for(m = 0; m < *n_bf; m++) {
      w[m] = 0.0;
      }
      /*first use w for breaks, then use w for weights.*/
      for(m = 0; m < length_breaks; m++) {
      for(p = 0; p < *n_preds; p++) {
      w[m] += break_matrix[m + length_breaks * p] * x[p];
      }
      }
      find_interval(&length_breaks, w, y, bin);
      for(m = 0; m < *n_bf; m++) {
      w[m] = 0.0;
      for(p = 0; p < *n_preds; p++) {
      w[m] += theta[m + *n_bf * p] * x[p];
      }
      }
      a = 0.0;
      b = 0.0;
      c = 0.0;
      d = w[0];
      for(m = 1; m < *n_bf; m++) {
      k = 4 * (m + *n_bf * (*bin + 2));
      a += spline_coefs[0 + k] * w[m];
      b += spline_coefs[1 + k] * w[m];
      c += spline_coefs[2 + k] * w[m];
      d += spline_coefs[3 + k] * w[m];
      }
      if(*u <= basis_knots[*bin + 2] || *u >= basis_knots[*bin + 3]) {
      *u = 0.5 * (basis_knots[*bin + 2] + basis_knots[*bin + 3]);
      }
      FindCubicRoot(&a, &b, &c, &d, y,
      &basis_knots[*bin + 2], &basis_knots[*bin + 3], u);
      *log_like = -1.0 * log(3 * a  * *u * *u + 2 * b * *u + c);
      Free(w);
      }
      
      void MiddleParametricDensity(int *n_bf, int *n_preds, int *basis_ind,
      double *y, double *u, double *log_like, int *bin, double *x, double *theta,
      double *shape, double *break_matrix,
      double *basis_knots, double *spline_coefs) {
      /******************************************************************
      # returns the log density and quantile level for the middle
      #   of the distribution.
      # Args:
      #   n_bf: number of basis functions.
      #   n_preds: number of predictors.
      #   basis_ind: parametric basis indicator.
      #   y: response.
      #   u: uniform random variable.
      #   log_like: log likelihood.
      #   bin: cell y is located in.
      #   x: n_preds length vector of covariates.
      #   theta: n_bf x n_preds matrix of regression coefficients.
      #   shape: shape parameter.
      #   break_matrix: n_bf x n_preds matrix of the basis functions times
      #      the regression coefficients.
      #   basis_knots: sequence of knots for the parametric basis.
      #   spline_coefs: 4 x n_bf x (n_bf - 1) array of cubic spline coefficients.
      # Updates u, log_like, bin.
      ******************************************************************/
      int m, p;
      double loc, scale, z;
      double double_zero = 0.0;
      double double_one = 1.0;
      double *breaks = (double *)Calloc(*n_bf, double);
      /*quantile function evaluated at the knots*/
      double *q_knots = (double *)Calloc(*n_bf, double);
      /*pointer to the log density function*/
      void (*dens_ptr)(double *, double *, double *, double *, double*);
      /*pointer to the distribution function*/
      void (*cdf_ptr)(double *, double *, double *, double *, double*);
      /*pointer to the quantile function*/
      void (*q_ptr)(double *, double *, double *, double *, double*);
      dens_ptr = &LogDNorm;
      cdf_ptr = &PNorm;
      q_ptr = &QNorm;
      
      switch(*basis_ind) {
      case 0:
      dens_ptr = &LogDNorm;
      cdf_ptr = &PNorm;
      q_ptr = &QNorm;
      break;
      case 1:
      dens_ptr = &LogDT;
      cdf_ptr = &PT;
      q_ptr = &QT;
      break;
      case 2:
      dens_ptr = &LogDLogistic;
      cdf_ptr = &PLogistic;
      q_ptr = &QLogistic;
      break;
      case 3:
      dens_ptr = &LogDAsymmetricLaplace;
      cdf_ptr = &PAsymmetricLaplace;
      q_ptr = &QAsymmetricLaplace;
      break;
      case 4:
      dens_ptr = &LogDWeibull;
      cdf_ptr = &PWeibull;
      q_ptr = &QWeibull;
      break;
      case 5:
      dens_ptr = &LogDGamma;
      cdf_ptr = &PGamma;
      q_ptr = &QGamma;
      break;
      }
      for(m = 0; m < *n_bf - 1; m++) {
      q_knots[m] = -999.0;
      }
      q_knots[*n_bf - 1] = 999.0;
      if(*n_bf > 2) {
      for(m = 1; m < *n_bf - 1; m++) {
      q_ptr(&basis_knots[m], &double_zero, &double_one, shape, &q_knots[m]);
      }
      }
      if(*n_bf == 2) {
      *bin = 0;
      } else {
      for(m = 0; m < *n_bf; m++) {
      breaks[m] = 0.0;
      }
      for(m = 0; m < *n_bf; m++) {
      for(p = 0; p < *n_preds; p++) {
      breaks[m] += break_matrix[m + *n_bf * p] * x[p];
      }
      }
      find_interval(n_bf, breaks, y, bin);
      }
      scale = 0.0;
      for(p = 0; p < *n_preds; p++) {
      scale += theta[(*bin + 1) + *n_bf * p] * x[p];
      }
      loc = 0.0;
      if(*bin == 0) {
      for(p = 0; p < *n_preds; p++) {
      loc += theta[*n_bf * p] * x[p];
      }
      } else {
      loc = breaks[*bin] - scale * q_knots[*bin];
      }
      z = (*y - loc) / scale;
      cdf_ptr(&z, &double_zero, &double_one, shape, u);
      dens_ptr(&z, &double_zero, &double_one, shape, log_like);
      *log_like -= log(scale);
      Free(breaks);
      }
      
      void SplineLikelihood(int *n_bf, int *n_preds, int *basis_ind, double *y,
      double *u, double *log_like,
      int *bin, double *x, double *theta, double *shape,
      double *break_matrix, double *basis_knots, double *spline_coefs,
      double *tau_low, double *tau_high,
      double *xi_low, double *xi_high,
      int *pareto_low, int *pareto_high,
      double *mn_low_p, double *mn_high_p,
      double *sd_low_p, double *sd_high_p,
      double *m_low_p, double *m_high_p,
      double *i_low_p, double *i_high_p) {
      /******************************************************************
      # returns the log density and quantile level for one observation
      #    using the spline basis.
      # Args:
      #   n_bf: number of basis functions.
      #   n_preds: number of predictors.
      #   basis_ind: parametric basis indicator.
      #   y: response.
      #   u: uniform random variable.
      #   bin: cell y is located in.
      #   x: n_preds length vector of covariates.
      #   theta: n_bf x n_preds matrix of regression coefficients.
      #   shape: shape parameter.
      #   break_matrix: n_bf x n_preds matrix of the basis functions times
      #     the regression coefficients.
      #   basis_knots: vector of knot coefficients.
      #   spline_coefs: 4 x n_bf x (n_bf - 1) array of cubic spline coefficients.
      #   tau_low: quantile level below which GEV tail is fit.
      #   xi_low: lower tail shape parameter.
      #   pareto_low: indicator if pareto tail is used.
      #   tau_high: quantile level above which the upper threshold is fit.
      #   xi_high: upper tail shape parameter.
      #   pareto_high: indicator if pareto tail is used.
      #   mn_low_p: vector of length p where t(mn_low_p) %*% x is the mean
      #     in the lower tail.
      #   mn_high_p: vector of length p where t(mn_high_p) %*% x is the mean
      #     in the upper tail.
      #     Assumed correct.
      #   sd_low_p: vector of length p where t(sd_low_p) %*% x is the sd
      #     in the lower tail.
      #     Assumed correct.
      #   sd_high_p: vector of length p where t(sd_high_p) %*% x is the sd
      #     in the upper tail.
      #     Assumed correct.
      #   m_low_p: vector of length p where t(m_low_p) %*% x is the spline sparsity
      #     in the lower tail.
      #     Assumed correct.
      #   m_high_p: vector of length p where t(m_high_p) %*% x is the
      #     spline sparsity in the upper tail.
      #     Assumed correct.
      #   i_low_p: vector of length p where t(i_low_p) %*% x is the lower threshold.
      #     Assumed correct.
      #   i_high_p: vector of length p where t(i_high_p) %*% x
      #     is the upper threshold.
      #     Assumed correct.
      # Updates: u, log_like, bin.
      *****************************************************************************/
      int p;
      double pareto_scale, thresh_low, thresh_high;
      
      thresh_low = 0.0;
      thresh_high = 0.0;
      for(p = 0; p < *n_preds; p++) {
      thresh_low += i_low_p[p] * x[p];
      thresh_high += i_high_p[p] * x[p];
      }
      if(*y <= thresh_low) {
      pareto_scale = 0.0;
      for(p = 0; p < *n_preds; p++) {
      pareto_scale += m_low_p[p] * x[p];
      }
      LowerParetoDensity(y, tau_low, &thresh_low,
      &pareto_scale, xi_low, pareto_low, u, log_like);
      } else if(*y >= thresh_high)  {
      pareto_scale = 0.0;
      for(p = 0; p < *n_preds; p++) {
      pareto_scale += m_high_p[p] * x[p];
      }
      UpperParetoDensity(y, tau_high, &thresh_high, &pareto_scale, xi_high,
      pareto_high, u, log_like);
      } else {
      MiddleSplineDensity(n_bf, n_preds, y, u, log_like, bin, x, theta,
      break_matrix, basis_knots, spline_coefs);
      }
      }
      
      void ParametricLikelihood(int *n_bf, int *n_preds, int *basis_ind, double *y,
      double *u, double *log_like, int *bin, double *x, double *theta,
      double *shape, double *break_matrix, double *basis_knots,
      double *spline_coefs, double *tau_low, double *tau_high,
      double *xi_low, double *xi_high,
      int *pareto_low, int *pareto_high,
      double *mn_low_p, double *mn_high_p,
      double *sd_low_p, double *sd_high_p,
      double *m_low_p, double *m_high_p,
      double *i_low_p, double *i_high_p) {
      /******************************************************************
      # returns the log density and quantile level for one observation
      #    using the parametric basis.
      # Args:
      #   n_bf: number of basis functions.
      #   n_preds: number of predictors.
      #   basis_ind: parametric basis indicator.
      #   y: response.
      #   u: uniform random variable.
      #   bin: cell y is located in.
      #   x: n_preds length vector of covariates.
      #   theta: n_bf x n_preds matrix of regression coefficients.
      #   shape: shape parameter.
      #   break_matrix: n_bf x n_preds matrix of the basis functions times
      #     the regression coefficients.
      #   spline_coefs: ignored.
      #     Not checked to be correct.
      #   basis_knots: vector of knot coefficients.
      #   tau_low: quantile level below which GEV tail is fit.
      #   xi_low: lower shape parameter.
      #   pareto_low: indicator if pareto tail is used.
      #   tau_high: quantile level above which the upper threshold is fit.
      #   xi_high: upper shape parameter.
      #   pareto_high: indicator if pareto tail is used.
      #   mn_low_p: vector of length p where t(mn_low_p) %*% x is the mean
      #     in the lower tail.
      #   mn_high_p: vector of length p where t(mn_high_p) %*% x is the mean
      #     in the upper tail.
      #     Assumed correct.
      #   sd_low_p: vector of length p where t(sd_low_p) %*% x is the sd
      #     in the lower tail.
      #     Assumed correct.
      #   sd_high_p: vector of length p where t(sd_high_p) %*% x is the sd
      #     in the upper tail.
      #     Assumed correct.
      #   m_low_p: vector of length p where t(m_low_p) %*% x is the spline sparsity
      #     in the lower tail.
      #     Assumed correct.
      #   m_high_p: vector of length p where t(m_high_p) %*% x
      #     is the spline sparsity in the upper tail.
      #     Assumed correct.
      #   i_low_p: vector of length p where t(i_low_p) %*% x is the lower threshold.
      #     Assumed correct.
      #   i_high_p: vector of length p where t(i_high_p) %*% x
      #     is the upper threshold.
      #     Assumed correct.
      # Updates: u, log_like, bin.
      *****************************************************************************/
      int p;
      double loc, scale, pareto_scale, thresh_low, thresh_high;
      /*pointer to the log density function*/
      void (*dens_ptr)(double *, double *, double *, double *, double*);
      dens_ptr = &LogDNorm;
      switch(*basis_ind) {
      case 0:
      dens_ptr = &LogDNorm;
      break;
      case 1:
      dens_ptr = &LogDT;
      break;
      case 2:
      dens_ptr = &LogDLogistic;
      break;
      case 3:
      dens_ptr = &LogDAsymmetricLaplace;
      break;
      case 4:
      dens_ptr = &LogDWeibull;
      break;
      case 5:
      dens_ptr = &LogDGamma;
      break;
      }
      
      thresh_low = 0.0;
      thresh_high = 0.0;
      for(p = 0; p < *n_preds; p++) {
      thresh_low += i_low_p[p] * x[p];
      thresh_high += i_high_p[p] * x[p];
      }
      if(*tau_low > 0.0 && *y <= thresh_low) {
      loc = 0.0;
      scale = 0.0;
      for(p = 0; p < *n_preds; p++) {
      loc += mn_low_p[p] * x[p];
      scale += sd_low_p[p] * x[p];
      }
      dens_ptr(&thresh_low, &loc, &scale, shape, &pareto_scale);
      pareto_scale = exp(-pareto_scale);
      pareto_scale *= *tau_low;
      LowerParetoDensity(y, tau_low, &thresh_low, &pareto_scale, xi_low,
      pareto_low, u, log_like);
      } else if(*tau_high < 1.0 && *y >= thresh_high)  {
      loc = 0.0;
      scale = 0.0;
      for(p = 0; p < *n_preds; p++) {
      loc += mn_high_p[p] * x[p];
      scale += sd_high_p[p] * x[p];
      }
      dens_ptr(&thresh_high, &loc, &scale, shape, &pareto_scale);
      pareto_scale = exp(-pareto_scale);
      pareto_scale *= (1.0 - *tau_high);
      UpperParetoDensity(y, tau_high, &thresh_high, &pareto_scale, xi_high,
      pareto_high, u, log_like);
      } else {
      MiddleParametricDensity(n_bf, n_preds, basis_ind, y, u, log_like,
      bin, x, theta, shape,
      break_matrix, basis_knots, spline_coefs);
      }
      }
      
      void Likelihood(int *dim_thetastar, int *n_y, int *basis_ind,
      double *basis_knots, double *y, double *y_low, double *y_high,
      int *status, double *x, double *u, double *u_low, double *u_high,
      int *bin, int *bin_low, int *bin_high, double *log_like, double *ll_sum,
      double *theta, double *shape, double *tau_low, double *tau_high,
      int *pareto_low, int *pareto_high, double *xi_low, double *xi_high,
      int *this_cluster, int *this_response) {
      /**************************************************************
      # returns the log density and quantile level for all observations
      #  for a given cluster by response combination.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   n_y: vector of length n_clusters containing the number of observations
      #    in each cluster.
      #  basis_ind: integer indicator where (-1, 0, 1, 2, 3, 4, 5) maps to
      #      ("spline", "Gaussian", "t", "logistic", "ALAP", "Weibull", "gamma").
      #  basis_knots: vector of knots.
      #   y: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #    of responses.
      #   y_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of left-censoring points.
      #   y_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of right-censoring points.
      #   status: array of dimension
      #     max(n.y) x n.timepoints x n.clusters x n_responses
      #     of status indicators where:
      #       -99 indicates that the value should be ignored;
      #       0 indicates that the value is missing and treated as continuous;
      #       1 indicates that the value is observed;
      #       2 indicates that the value is right-censored;
      #       3 indicates that the value is left-censored;
      #       4 indicates that the value is interval-censored.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   u: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables.
      #   u_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_low.
      #   u_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_high.
      #   bin: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u is in.
      #   bin_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u_low is in.
      #   bin_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u_high is in.
      #   log_like: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of loglikelihood (or censored loglikelihood) values.
      #   ll_sum: scalar summed loglikelihood (or censored loglikelihood) values.
      #   theta: n_bf x n_preds x n_timepoints x n_clusters x n_repsonses array
      #    of thresholded regression parameters.
      #   shape: n_responses length vector of shape parameters for bases.
      #   tau_low: the quantile level below which a parametric tail is fit.
      #   tau_high: the quantile level above which a parametric tail is fit.
      #   pareto_low: indicator if pareto tail is fit below tau_low.
      #   pareto_high: indicator if pareto tail is fit above tau_high.
      #   xi_low: array of dimension n_clusters x n_responses
      #     of lower tail parameters.
      #   xi_high: array of dimension n_clusters x n_responses
      #     of upper tail parameters.
      #   this_cluster: cluster indicator.
      #   this_response: response indicator.
      #   c_indicator: an indicator if the calculation should be performed in C.
      # Updates u, u_low, u_high, bin, bin_low, bin_high, loglike, ll_sum.
      **************************************************************/
      int i, j, k, m, p;
      int c_x, this_observation, this_status;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_timepoints = dim_thetastar[5];
      int n_bf_preds = n_bf * n_preds;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int c_y = max_n_y * n_timepoints *
      (*this_cluster + n_clusters * *this_response);
      int int_one = 1;
      int length_spline_coefs  = 1;
      if(*basis_ind == -1) {
      length_spline_coefs = 4 * n_bf * (n_bf - 1);
      }
      int length_interior_knots = 1;
      int nrow_this_basis = n_bf;
      if(*basis_ind == -1) {
      length_interior_knots = n_bf - 4;
      nrow_this_basis = n_bf - 2;
      }
      if(*basis_ind != -1) {
      if(n_bf > 2) {
      length_interior_knots = n_bf - 2;
      }
      }
      int length_interior_knot_basis = length_interior_knots * n_bf;
      int length_this_basis = nrow_this_basis * n_bf;
      int length_break_matrix = nrow_this_basis * n_preds;
      int spline_knots_length = n_bf + 4;
      int dummy_bin = 2;
      int int_zero = 0;
      
      double double_zero = 0.0;
      double double_one = 1.0;
      
      double *spline_coefs = (double *)Calloc(length_spline_coefs, double);
      /*quantile function evaluated at the knots*/
      double *q_knots = (double *)Calloc(n_bf, double);
      double *this_basis = (double *)Calloc(length_this_basis, double);
      double *interior_knots = (double *)Calloc(length_interior_knots, double);
      double *this_theta = (double *)Calloc(n_bf * n_preds, double);
      double *break_matrix = (double *)Calloc(length_break_matrix, double);
      double *this_x = (double *)Calloc(n_preds, double);
      double *interior_knot_basis =
      (double *)Calloc(length_interior_knot_basis, double);
      
      double *m_low = (double *)Calloc(n_bf, double);
      double *m_high = (double *)Calloc(n_bf, double);
      double *i_low = (double *)Calloc(n_bf, double);
      double *i_high = (double *)Calloc(n_bf, double);
      double *m_low_p = (double *)Calloc(n_preds, double);
      double *m_high_p = (double *)Calloc(n_preds, double);
      double *i_low_p = (double *)Calloc(n_preds, double);
      double *i_high_p = (double *)Calloc(n_preds, double);
      double *mn_low_p = (double *)Calloc(n_preds, double);
      double *mn_high_p = (double *)Calloc(n_preds, double);
      double *sd_low_p = (double *)Calloc(n_preds, double);
      double *sd_high_p = (double *)Calloc(n_preds, double);
      /*pointer to the quantile function*/
      void (*q_ptr)(double *, double *, double *, double *, double*);
      /*pointer to the parametric or spline likelihood*/
      void (*like_ptr)(int *, int *, int *, double *, double *, double *,
      int *, double *, double *, double *,
      double *, double *, double *,
      double *, double *,
      double *, double *,
      int *, int *,
      double *, double *,
      double *, double *,
      double *, double *,
      double *, double *);
      
      spline_coefs[0] = 0.0;
      if(*basis_ind == -1) {
      for(i = 1; i < length_spline_coefs; i++) {
      spline_coefs[i] = 0.0;
      }
      MakeCubicSplineCoefficients (&n_bf, basis_knots, spline_coefs);
      }
      
      
      for(m = 0; m < (n_bf - 1); m++) {
      q_knots[m] = -999.0;
      }
      q_knots[n_bf - 1] = 999.0;
      
      q_ptr = &QNorm;
      switch(*basis_ind) {
      case 0:
      q_ptr = &QNorm;
      break;
      case 1:
      q_ptr = &QT;
      break;
      case 2:
      q_ptr = &QLogistic;
      break;
      case 3:
      q_ptr = &QAsymmetricLaplace;
      break;
      case 4:
      q_ptr = &QWeibull;
      break;
      case 5:
      q_ptr = &QGamma;
      break;
      }
      if(*basis_ind != -1)  {
      if(n_bf > 2) {
      for(m = 1; m < n_bf - 1; m++) {
      q_ptr(&basis_knots[m], &double_zero, &double_one, shape, &q_knots[m]);
      }
      }
      }
      like_ptr = &SplineLikelihood;
      if(*basis_ind != -1) {
      like_ptr = &ParametricLikelihood;
      }
      for(i = 0; i < length_this_basis; i++) {
      this_basis[i] = 0.0;
      }
      
      for(i = 0; i < length_interior_knots; i++) {
      interior_knots[i] = 0.0;
      }
      for(i = 0; i < n_bf_preds; i++) {
      this_theta[i] = 0.0;
      }
      for(i = 0; i < length_break_matrix; i++) {
      break_matrix[i] = 0.0;
      }
      
      
      if(n_bf > 2) {
      for(i = 0; i < length_interior_knot_basis; i++) {
      interior_knot_basis[i] = 0.0;
      }
      if(*basis_ind == -1) {
      for(m = 0; m < length_interior_knots; m++) {
      interior_knots[m] = basis_knots[m + 3];
      }
      } else {
      for(m = 0; m < length_interior_knots; m++) {
      interior_knots[m] = basis_knots[m + 1];
      }
      }
      MakeBasis(basis_ind, &n_bf, basis_knots, &length_interior_knots,
      interior_knots, shape, interior_knot_basis);
      /*copy rows of basis matrix down*/
      for(i = 1; i < (nrow_this_basis - 1); i++) {
      for(m = 0; m < n_bf; m++) {
      this_basis[i + nrow_this_basis * m] =
      interior_knot_basis[(i - 1) + length_interior_knots * m];
      }
      }
      }
      
      /*fill in top and bottom rows of this_basis*/
      for(m = 0; m < n_bf; m++) {
      this_basis[0 + nrow_this_basis * m] = -999999.0;
      this_basis[(nrow_this_basis - 1) + nrow_this_basis * m] = 999999.0;
      }
      for(p = 0; p < n_preds; p++) {
      m_low_p[p] = 0.0;
      m_high_p[p] = 0.0;
      i_low_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      i_high_p[p] = 0.0;
      }
      for(m = 0; m < n_bf; m++) {
      m_low[m] = 0.0;
      m_high[m] = 0.0;
      i_low[m] = 0.0;
      i_high[m] = 0.0;
      }
      i_low[0] = 1.0;
      i_high[0] = 1.0;
      
      if(*tau_low > 0.0) {
      MakeBasis(basis_ind, &n_bf, basis_knots, &int_one, tau_low, shape, i_low);
      for(m = 0; m < n_bf; m++) {
      this_basis[nrow_this_basis * m] = i_low[m];
      }
      if(*basis_ind == -1) {
      for(m = 0; m < (n_bf - 1); m++) {
      CubicMSpline2(tau_low, &spline_knots_length, basis_knots, spline_coefs,
      &m, &dummy_bin, &m_low[m + 1]);
      m_low[m + 1] *= *tau_low;
      }
      }
      }
      
      if(*tau_high < 1.0) {
      MakeBasis(basis_ind, &n_bf, basis_knots, &int_one, tau_high, shape, i_high);
      for(m = 0; m < n_bf; m++) {
      this_basis[(nrow_this_basis - 1) + nrow_this_basis * m] = i_high[m];
      }
      if(*basis_ind == -1) {
      for(m = 0; m < (n_bf - 1); m++) {
      CubicMSpline2(tau_high, &spline_knots_length, basis_knots, spline_coefs,
      &m, &dummy_bin, &m_high[m + 1]);
      m_high[m + 1] *= (1 - *tau_high);
      }
      }
      }
      *ll_sum = 0.0;
      for(j = 0; j < n_timepoints; j++) {
      c_x = max_n_y * n_preds * (j + n_timepoints * *this_cluster);
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      k = m + n_bf * p;
      this_theta[k] = theta[k + n_bf * n_preds * n_timepoints *
      (*this_cluster + n_clusters * *this_response)];
      }
      }
      matrix_multiply(&int_zero, &int_zero, &nrow_this_basis, &n_bf, &n_bf,
      &n_preds, this_basis, this_theta, break_matrix);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf, &int_one,
      this_theta, m_low, m_low_p);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf, &int_one,
      this_theta, m_high, m_high_p);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf, &int_one,
      this_theta, i_low, i_low_p);
      matrix_multiply(&int_one, &int_zero, &n_bf, &n_preds, &n_bf, &int_one,
      this_theta, i_high, i_high_p);
      
      for(p = 0; p < n_preds; p++) {
      mn_low_p[p] = this_theta[0 + n_bf * p];
      sd_low_p[p] = this_theta[1 + n_bf * p];
      mn_high_p[p] = this_theta[0 + n_bf * p];
      sd_high_p[p] = this_theta[(n_bf - 1) + n_bf * p];
      }
      if(n_bf > 2) {
      for(m = 1; m < (n_bf - 1); m++) {
      for(p = 0; p < n_preds; p++) {
      mn_high_p[p] += (this_theta[(m + 1) + n_bf * p] -
      this_theta[m + n_bf * p]) * q_knots[m];
      }
      }
      }
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(p = 0; p < n_preds; p++) {
      this_x[p] = x[i + max_n_y * p + c_x];
      }
      this_observation = i + max_n_y * j + c_y;
      this_status = status[this_observation];
      if(this_status == -99) {
      log_like[this_observation] = 0.0;
      }
      if(this_status == 0 || this_status == 1) {
      like_ptr(&n_bf, &n_preds, basis_ind, &y[this_observation],
      &u[this_observation], &log_like[this_observation],
      &bin[this_observation], this_x, this_theta, shape,
      break_matrix, basis_knots, spline_coefs,
      tau_low, tau_high,
      xi_low, xi_high,
      pareto_low, pareto_high,
      mn_low_p, mn_high_p,
      sd_low_p, sd_high_p,
      m_low_p, m_high_p,
      i_low_p, i_high_p);
      }
      if(this_status == 2) {
      like_ptr(&n_bf, &n_preds, basis_ind, &y_low[this_observation],
      &u_low[this_observation], &log_like[this_observation],
      &bin_low[this_observation], this_x, this_theta, shape,
      break_matrix, basis_knots, spline_coefs,
      tau_low, tau_high,
      xi_low, xi_high,
      pareto_low, pareto_high,
      mn_low_p, mn_high_p,
      sd_low_p, sd_high_p,
      m_low_p, m_high_p,
      i_low_p, i_high_p);
      }
      if(this_status == 3) {
      like_ptr(&n_bf, &n_preds, basis_ind, &y_high[this_observation],
      &u_high[this_observation], &log_like[this_observation],
      &bin_high[this_observation], this_x, this_theta, shape,
      break_matrix, basis_knots, spline_coefs,
      tau_low, tau_high,
      xi_low, xi_high,
      pareto_low, pareto_high,
      mn_low_p, mn_high_p,
      sd_low_p, sd_high_p,
      m_low_p, m_high_p,
      i_low_p, i_high_p);
      }
      if(this_status == 4) {
      like_ptr(&n_bf, &n_preds, basis_ind, &y_low[this_observation],
      &u_low[this_observation], &log_like[this_observation],
      &bin_low[this_observation], this_x, this_theta, shape,
      break_matrix, basis_knots, spline_coefs,
      tau_low, tau_high,
      xi_low, xi_high,
      pareto_low, pareto_high,
      mn_low_p, mn_high_p,
      sd_low_p, sd_high_p,
      m_low_p, m_high_p,
      i_low_p, i_high_p);
      
      like_ptr(&n_bf, &n_preds, basis_ind, &y_high[this_observation],
      &u_high[this_observation], &log_like[this_observation],
      &bin_high[this_observation], this_x, this_theta, shape,
      break_matrix, basis_knots, spline_coefs,
      tau_low, tau_high,
      xi_low, xi_high,
      pareto_low, pareto_high,
      mn_low_p, mn_high_p,
      sd_low_p, sd_high_p,
      m_low_p, m_high_p,
      i_low_p, i_high_p);
      }
      if(this_status > 1) {
      log_like[this_observation] =
      log(u_high[this_observation] - u_low[this_observation]);
      RandomUniform(&u_low[this_observation],
      &u_high[this_observation], &u[this_observation]);
      }
      *ll_sum += log_like[this_observation];
      }
      }
      Free(spline_coefs);
      Free(q_knots);
      Free(this_basis);
      Free(interior_knots);
      Free(this_theta);
      Free(break_matrix);
      Free(this_x);
      Free(interior_knot_basis);
      Free(m_low);
      Free(m_high);
      Free(i_low);
      Free(i_high);
      Free(m_low_p);
      Free(m_high_p);
      Free(i_low_p);
      Free(i_high_p);
      Free(mn_low_p);
      Free(mn_high_p);
      Free(sd_low_p);
      Free(sd_high_p);
      }
      
      void MakeZ(int *dim_thetastar, int *n_preds_q, double *x, double *z) {
      /**************************************************************
      # Turns x into a regression matrix for the copula.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   x: n_timepoints x n_preds_q matrix of predictors.
      #   z: (n_timepoints * n_responses) x (n_preds_q x n_responses)
      #     matrix of predictors.
      # Updates z.
      **************************************************************/
      int h, j, k, l;
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int nrow_z = n_responses * n_timepoints;
      
      for(j = 0; j < n_timepoints; j++) {
      for(h = 0; h < n_responses; h++) {
      k = j * n_responses + h;
      for(l = 0; l < *n_preds_q; l++) {
      z[k + nrow_z * (h * *n_preds_q + l)] = x[j + n_timepoints * l];
      }
      }
      }
      }
      
      void MakeD2(int *dim_thetastar, int *n_preds_q, int *n_y, double *x,
      double *delta, double *lambda, double *alpha,
      double *distance_matrix_j, double *d2, int *this_cluster) {
      /**************************************************************
      # Creates the copula weights to ensure uniform random variables.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   d2: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      #   this_cluster: cluster indicator.
      # Updates d2.
      **************************************************************/
      int i, j, k, q;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_responses_2 = n_responses * n_responses;
      int n_timepoints_2 = n_timepoints * n_timepoints;
      int n_timepoints_responses = n_timepoints * n_responses;
      int n_col_z = n_responses * *n_preds_q;
      int length_cov_eta = n_timepoints_responses * n_timepoints_responses;
      int length_this_x = n_timepoints * *n_preds_q;
      int length_z = n_timepoints_responses * n_col_z;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      double this_z;
      double *ar_matrix = (double *)Calloc(n_timepoints_2, double);
      double *this_lambda = (double *)Calloc(n_responses_2, double);
      double *cov_eta = (double *)Calloc(length_cov_eta, double);
      double *this_x = (double *)Calloc(length_this_x, double);
      double *z = (double *)Calloc(length_z, double);
      double *v = (double *)Calloc(length_cov_eta, double);
      
      for(i = 0; i < n_timepoints_2; i++) {
      ar_matrix[i] = pow(alpha[*this_cluster], distance_matrix_j[i]);
      }
      for(i = 0; i < n_responses_2; i++) {
      this_lambda[i] = lambda[*this_cluster + n_clusters * i];
      }
      for(i = 0; i < length_cov_eta; i++) {
      cov_eta[i] = 0.0;
      v[i] = 0.0;
      }
      Kronecker(&n_timepoints, &n_timepoints, &n_responses, &n_responses,
      ar_matrix, this_lambda, cov_eta);
      for(i = 0; i < length_this_x; i++) {
      this_x[i] = 0.0;
      }
      for(i = 0; i < length_z; i++) {
      z[i] = 0.0;
      }
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(q = 0; q < *n_preds_q; q++) {
      this_x[j + n_timepoints * q] = x[i + max_n_y * q +
      max_n_y * n_preds * (j + n_timepoints * *this_cluster)];
      }
      }
      MakeZ(dim_thetastar, n_preds_q, this_x, z);
      for(j = 0; j < n_timepoints_responses; j++) {
      k = i + max_n_y * *this_cluster + max_n_y * n_clusters * j;
      d2[k] = 0.0;
      for(q = 0; q < n_col_z; q++) {
      this_z = z[j + n_timepoints_responses * q];
      d2[k] += this_z * this_z * delta[*this_cluster + n_clusters * q];
      }
      d2[k] += cov_eta[j + n_timepoints_responses * j] + 1.0;
      d2[k] = pow(d2[k], -0.5);
      }
      }
      Free(ar_matrix);
      Free(this_lambda);
      Free(cov_eta);
      Free(this_x);
      Free(z);
      Free(v);
      }
      
      void GenerateU(int *dim_thetastar, int *n_preds_q, int *n_y,
      int *status, double *x,
      double *u, double *u2, double *w, double *w2,
      double *gamma, double *eta, double *d2) {
      /**********************************************************
      # generates uniform random variables from the copula.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   status: array of dimension
      #     max_n_y x n_timepoints x n_clusters x n_responses
      #     of status indicators where:
      #       -99 indicates that the value should be ignored;
      #       0 indicates that the value is missing and treated as continuous;
      #       1 indicates that the value is observed;
      #       2 indicates that the value is right-censored;
      #       3 indicates that the value is left-censored;
      #       4 indicates that the value is interval-censored.
      #   x: array of dimension max_n_y x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   u: array of dimension max_n_y x n_timepoints x n_clusters x n_responses
      #     of uniform random variables.
      #   u2: array of dimension max_n_y x n_clusters x (n_timepoints * n_responses)
      #     of uniform random variables.  Not checked to be equal to u
      #     where appropriate.
      #   w: array of dimension max_n_y x n_timepoints x n_clusters x n_responses
      #     of qnorm(u).  Not checked to be equal to u where appropriate.
      #   w2: array of dimension max_n_y x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).  Not checked to be equal to u2 where appropriate.
      #  gamma: array of dimension max_n_y x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #  eta: array of dimension max_n_y x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #  d2: array of dimension max_n_y x n_timepoints x n_clusters x n_responses
      #    of scaling factors.
      # Updates u, u2, w, w2.
      **********************************************************/
      int g, h, i, j, k, k2, l, q;
      int int_zero = 0;
      int int_one = 1;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_timepoints_responses = n_timepoints * n_responses;
      int n_col_z = *n_preds_q * n_responses;
      int length_z = n_timepoints_responses * n_col_z;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_this_x = n_timepoints * *n_preds_q;
      
      double double_zero = 0.0;
      double double_one = 1.0;
      double *z = (double *)Calloc(length_z, double);
      double *this_w = (double *)Calloc(n_timepoints_responses, double);
      double *this_gamma = (double *)Calloc(n_col_z, double);
      double *this_x = (double *)Calloc(length_this_x, double);
      
      for(i = 0; i < length_z; i++) {
      z[i] = 0.0;
      }
      for(i = 0; i < n_timepoints_responses; i++) {
      this_w[i] = 0.0;
      }
      for(i = 0; i < n_col_z; i++) {
      this_gamma[i] = 0.0;
      }
      for(i = 0; i < length_this_x; i++) {
      this_x[i] = 0.0;
      }
      
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_col_z; j++) {
      this_gamma[j] = gamma[i + max_n_y * (g + n_clusters * j)];
      }
      for(j = 0; j < n_timepoints; j++) {
      for(q = 0; q < *n_preds_q; q++) {
      this_x[j + n_timepoints * q] =
      x[i + max_n_y * q + max_n_y * n_preds * (j + n_timepoints * g)];
      }
      }
      MakeZ(dim_thetastar, n_preds_q, this_x, z);
      matrix_multiply(&int_zero, &int_zero, &n_timepoints_responses,
      &n_col_z, &n_col_z, &int_one, z, this_gamma, this_w);
      for(j = 0; j < n_timepoints_responses; j++) {
      k = i + max_n_y * (g + n_clusters * j);
      this_w[j] += eta[k] + rnorm(0, 1);
      this_w[j] *= d2[k];
      }
      for(j = 0; j < n_timepoints; j++) {
      for(h = 0; h < n_responses; h++) {
      k = i + max_n_y * j + max_n_y * n_timepoints * (g + n_clusters * h);
      if(status[k] == -99 || status[k] == 0) {
      l = j * n_responses + h;
      k2 = i + max_n_y * g + max_n_y * n_clusters * l;
      w2[k2] = this_w[l];
      w[k] = this_w[l];
      PNorm(&this_w[l], &double_zero, &double_one, &double_one, &u2[k2]);
      u[k] = u2[k2];
      }
      }
      }
      }
      }
      PutRNGstate();
      Free(z);
      Free(this_w);
      Free(this_gamma);
      Free(this_x);
      }
      
      void CopyUniform(int *dim_thetastar, int *n_y,
      double *u, double *u2, double *w, double *w2,
      int *this_cluster, int *to_copula) {
      /**********************************************************
      # generates uniform random variables from the copula.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   u: array of dimension max_n_y x n_timepoints x n_clusters x n_responses
      #     of uniform random variables.
      #   u2: array of dimension max_n_y x n_clusters x (n_timepoints * n_responses)
      #     of uniform random variables.  Not checked to be equal to u
      #     where appropriate.
      #   w: array of dimension max_n_y x n_timepoints x n_clusters x n_responses
      #     of qnorm(u).  Not checked to be equal to u where appropriate.
      #   w2: array of dimension max_n_y x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).  Not checked to be equal to u2 where appropriate.
      #   this_cluster: indicator for which cluster to copy in.
      #   to_copula: binary indicator
      #     if the data variables are copied to the copula, or vice versa.
      # Updates u, u2, w, w2.
      **********************************************************/
      int i, j, k, l1, l2;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      
      if(*to_copula == 1) {
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(k = 0; k < n_responses; k++) {
      l1 = i + max_n_y * j +
      max_n_y * n_timepoints * (*this_cluster + n_clusters * k);
      l2 = i + max_n_y * *this_cluster +
      max_n_y * n_clusters * (j * n_responses + k);
      u2[l2] = u[l1];
      w2[l2] = w[l1];
      }
      }
      }
      } else {
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(k = 0; k < n_responses; k++) {
      l1 = i + max_n_y * j +
      max_n_y * n_timepoints * (*this_cluster + n_clusters * k);
      l2 = i + max_n_y * *this_cluster +
      max_n_y * n_clusters * (j * n_responses + k);
      u[l1] = u2[l2];
      w[l1] = w2[l2];
      }
      }
      }
      }
      }
      
      void UpdateGamma(int *dim_thetastar, int *n_preds_q, int *n_y, double *x,
      double *w2, double *gamma, double *delta, double *eta,
      double *d2) {
      
      /**********************************************************
      # updates random regression coefficients in the copula.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   w2: array of dimension max(n_y) x n_clusters x
      #     (n_timepoints * n_responses)
      #     of qnorm(u2).  Not checked to be equal to u2 where appropriate.
      #  gamma: array of dimension max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #  delta: n_clusters x n_responses x n_preds_q array of variances.
      #  eta: array of dimension max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   d2: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      # Updates gamma.
      **********************************************************/
      int g, i, j, k, q;
      int int_zero = 0;
      int int_one = 1;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int n_row_z = n_timepoints * n_responses;
      int n_col_z = *n_preds_q * n_responses;
      int n_col_z_2 = n_col_z * n_col_z;
      int length_z = n_row_z * n_col_z;
      int length_this_x = n_timepoints * *n_preds_q;
      
      double dummy_double = 0.0;
      double *z = (double *)Calloc(length_z, double);
      double *prec_gamma = (double *)Calloc(n_col_z_2, double);
      double *cov_gamma = (double *)Calloc(n_col_z_2, double);
      double *omega = (double *)Calloc(n_row_z, double);
      double *omega2 = (double *)Calloc(n_col_z, double);
      double *this_gamma = (double *)Calloc(n_col_z, double);
      double *this_x = (double *)Calloc(length_this_x, double);
      
      for(i = 0; i < length_z; i++) {
      z[i] = 0.0;
      }
      for(i = 0; i < n_col_z_2; i++) {
      prec_gamma[i] = 0.0;
      cov_gamma[i] = 0.0;
      }
      for(i = 0; i < n_row_z; i++) {
      omega[i] = 0.0;
      }
      for(i = 0; i < n_col_z; i++) {
      omega2[i] = 0.0;
      this_gamma[i] = 0.0;
      }
      for(i = 0; i < length_this_x; i++) {
      this_x[i] = 0.0;
      }
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(q = 0; q < *n_preds_q; q++) {
      this_x[j + n_timepoints * q] =
      x[i + max_n_y * q + max_n_y * n_preds * (j + n_timepoints * g)];
      }
      }
      MakeZ(dim_thetastar, n_preds_q, this_x, z);
      matrix_multiply(&int_one, &int_zero, &n_row_z, &n_col_z,
      &n_row_z, &n_col_z, z, z, prec_gamma);
      for(q = 0; q < n_col_z; q++) {
      prec_gamma[q + n_col_z * q] += 1.0 / delta[g + n_clusters * q];
      }
      InvertSymmetricMatrix(&n_col_z, prec_gamma, cov_gamma, &dummy_double);
      for(q = 0; q < n_row_z; q++) {
      k = i + max_n_y * (g + n_clusters * q);
      omega[q] = w2[k]/ d2[k] - eta[k] ;
      }
      matrix_multiply(&int_one, &int_zero, &n_row_z, &n_col_z, &n_row_z,
      &int_one, z, omega, omega2);
      matrix_multiply(&int_zero, &int_zero, &n_col_z, &n_col_z, &n_col_z,
      &int_one, cov_gamma, omega2, omega);
      RandomNormal(&int_one, &n_col_z, omega, cov_gamma, this_gamma);
      for(q = 0; q < n_col_z; q++) {
      gamma[i + max_n_y * (g + n_clusters * q)] = this_gamma[q];
      }
      }
      }
      PutRNGstate();
      Free(z);
      Free(prec_gamma);
      Free(cov_gamma);
      Free(omega);
      Free(omega2);
      Free(this_x);
      }
      
      void UpdateEta(int *dim_thetastar, int *n_preds_q, int *n_y, double *x,
      double *w2, double *gamma, double *eta,
      double *lambda, double *alpha, double *distance_matrix_j, double *d2) {
      /**********************************************************
      # updates random regression coefficients in the copula.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   w2: array of dimension
      #     max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).  Not checked to be equal to u2 where appropriate.
      #   gamma: array of dimension
      #    max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #   eta: array of dimension max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses of scaling factors.
      # Updates eta.
      **********************************************************/
      int g, i, j, k, q;
      int int_zero = 0;
      int int_one = 1;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int n_row_z = n_timepoints * n_responses;
      int n_col_z = *n_preds_q * n_responses;
      int n_row_z_2 = n_row_z * n_row_z;
      int length_z = n_row_z * n_col_z;
      int length_this_x = n_timepoints * *n_preds_q;
      int length_ar_matrix = n_timepoints * n_timepoints;
      int length_this_lambda = n_responses * n_responses;
      
      double dummy_double = 0.0;
      double *z = (double *)Calloc(length_z, double);
      double *prec_eta = (double *)Calloc(n_row_z_2, double);
      double *cov_eta = (double *)Calloc(n_row_z_2, double);
      double *this_eta = (double *)Calloc(n_row_z, double);
      double *omega = (double *)Calloc(n_row_z, double);
      double *this_x = (double *)Calloc(length_this_x, double);
      double *ar_matrix = (double *)Calloc(length_ar_matrix, double);
      double *this_gamma = (double *)Calloc(n_col_z, double);
      double *this_lambda = (double *)Calloc(length_this_lambda, double);
      
      for(i = 0; i < length_z; i++) {
      z[i] = 0.0;
      }
      for(i = 0; i < n_row_z_2; i++) {
      prec_eta[i] = 0.0;
      cov_eta[i] = 0.0;
      }
      for(i = 0; i < n_row_z; i++) {
      this_eta[i] = 0.0;
      omega[i] = 0.0;
      }
      for(i = 0; i < length_this_x; i++) {
      this_x[i] = 0.0;
      }
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < length_ar_matrix; i++) {
      ar_matrix[i] = pow(alpha[g], distance_matrix_j[i]);
      }
      for(i = 0; i < length_this_lambda; i++) {
      this_lambda[i] = lambda[g + n_clusters * i];
      }
      Kronecker(&n_timepoints, &n_timepoints, &n_responses, &n_responses,
      ar_matrix, this_lambda, cov_eta);
      InvertSymmetricMatrix(&n_row_z, cov_eta, prec_eta, &dummy_double);
      for(q = 0; q < n_row_z; q++) {
      prec_eta[q + n_row_z * q] += 1.0;
      }
      InvertSymmetricMatrix(&n_row_z, prec_eta, cov_eta, &dummy_double);
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(q = 0; q < *n_preds_q; q++) {
      this_x[j + n_timepoints * q] =
      x[i + max_n_y * q + max_n_y * n_preds * (j + n_timepoints * g)];
      }
      }
      MakeZ(dim_thetastar, n_preds_q, this_x, z);
      for(q = 0; q < n_col_z; q++) {
      this_gamma[q] = gamma[i + max_n_y * (g + n_clusters * q)];
      //        if(i == 0) {
      //          Rprintf("this_gamma[q] = %f\n", this_gamma[q]);
      //        }
      }
      matrix_multiply(&int_zero, &int_zero, &n_row_z, &n_col_z, &n_col_z,
      &int_one, z, this_gamma, this_eta);
      
      //        if(i == 0) {
      //          Rprintf("this_eta0 = %f\n", this_eta[n_row_z - 1]);
      //        }
      
      for(q = 0; q < n_row_z; q++) {
      k = i + max_n_y * (g + n_clusters * q);
      this_eta[q] *= -1.0;
      this_eta[q] += w2[k] / d2[k];
      }
      //        if(i == 0) {
      //          Rprintf("this_eta1 = %f\n", this_eta[n_row_z - 1]);
      //        }
      
      
      
      matrix_multiply(&int_zero, &int_zero, &n_row_z, &n_row_z, &n_row_z,
      &int_one, cov_eta, this_eta, omega);
      
      if(i == 0) {
      
      //          double zzz1 = 1.0;
      //          double zzz2 = 1.0;
      //
      //          for(q = 0; q < max_n_y * n_clusters * *n_preds_q * n_timepoints; q++) {
      //              if(d2[q] < zzz1) {
      //                zzz1 = d2[q];
      //              }
      //              if(d2[q] > zzz2) {
      //                zzz2 = d2[q];
      //              }
      //
      //          }
      //          Rprintf("min d max d this_eta[q] w[q] = %f %f %f \n", zzz1, zzz2, this_eta[n_row_z - 1], w[]);
      
      
      
      
      //        for(q = 0; q < n_timepoints * n_responses * n_timepoints * n_responses; q++) {
      //          Rprintf("cov_eta[q] = %f\n", cov_eta[q]);
      //        }
      
      
      }
      
      
      
      RandomNormal(&int_one, &n_row_z, omega, cov_eta, this_eta);
      
      if(i == 0) {
      //        for(q = 0; q < n_timepoints * n_responses; q++) {
      //          Rprintf("omega[q] = %f\n", omega[q]);
      //        }
      
      //        for(q = 0; q < n_timepoints * n_responses * n_timepoints * n_responses; q++) {
      //          Rprintf("cov_eta[q] = %f\n", cov_eta[q]);
      //        }
      
      
      }
      
      for(q = 0; q < n_row_z; q++) {
      eta[i + max_n_y * (g + n_clusters * q)] = this_eta[q];
      }
      }
      }
      PutRNGstate();
      Free(z);
      Free(prec_eta);
      Free(cov_eta);
      Free(this_eta);
      Free(omega);
      Free(this_x);
      Free(ar_matrix);
      Free(this_gamma);
      Free(this_lambda);
      }
      
      void EtaLogLikelihood(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *eta,
      double *lambda, double *alpha, double *distance_matrix_j,
      int *this_cluster, double *eta_ll) {
      /**********************************************************
      # Prior for eta.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #    in each cluster.
      #   eta: array of dimension max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   this_cluster: cluster indicator.
      #   eta_ll: likelihood for eta.
      # Updates eta_ll.
      **********************************************************/
      int i, q;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_eta = n_timepoints * n_responses;
      int length_eta_2 = length_eta * length_eta;
      int length_ar_matrix = n_timepoints * n_timepoints;
      int length_this_lambda = n_responses * n_responses;
      
      double log_det = 0.0;
      double qf = 0.0;
      double this_qf = 0.0;
      double *prec_eta = (double *)Calloc(length_eta_2, double);
      double *cov_eta = (double *)Calloc(length_eta_2, double);
      double *this_eta = (double *)Calloc(length_eta, double);
      double *ar_matrix = (double *)Calloc(length_ar_matrix, double);
      double *this_lambda = (double *)Calloc(length_this_lambda, double);
      
      for(i = 0; i < length_eta_2; i++) {
      prec_eta[i] = 0.0;
      cov_eta[i] = 0.0;
      }
      for(i = 0; i < length_eta; i++) {
      this_eta[i] = 0.0;
      }
      for(i = 0; i < length_ar_matrix; i++) {
      ar_matrix[i] = pow(alpha[*this_cluster], distance_matrix_j[i]);
      }
      for(i = 0; i < length_this_lambda; i++) {
      this_lambda[i] = lambda[*this_cluster + n_clusters * i];
      }
      Kronecker(&n_timepoints, &n_timepoints, &n_responses, &n_responses,
      ar_matrix, this_lambda, cov_eta);
      InvertSymmetricMatrix(&length_eta, cov_eta, prec_eta, &log_det);
      
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(q = 0; q < length_eta; q++) {
      this_eta[q] = eta[i + max_n_y * (*this_cluster + n_clusters * q)];
      }
      QuadraticForm(&length_eta, this_eta, prec_eta, &this_qf);
      qf += this_qf;
      }
      *eta_ll = 0.5 * (n_y[*this_cluster] * log_det - qf);
      Free(prec_eta);
      Free(cov_eta);
      Free(this_eta);
      Free(ar_matrix);
      Free(this_lambda);
      }
      
      void GammaLogLikelihood(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *gamma, double *delta,
      int *this_cluster, double *gamma_ll) {
      /**********************************************************
      # Prior for gamma.
      # Args:
      #   dim_thetastar: a vector of length 6 containing
      #     the dimension of thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters
      #     containing the number of observations in each cluster.
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   this_cluster: cluster indicator.
      # Returns the prior for gamma.
      # Updates eta_ll.
      **********************************************************/
      int i, q;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_gamma = n_responses * *n_preds_q;
      
      double qf = 0.0;
      double log_det = 0.0;
      double *prec_delta = (double *)Calloc(length_gamma, double);
      double *this_gamma = (double *)Calloc(length_gamma, double);
      
      for(i = 0; i < length_gamma; i++) {
      prec_delta[i] = 1.0 / delta[*this_cluster + n_clusters * i];
      log_det += log(prec_delta[i]);
      }
      
      for(i = 0; i < n_y[*this_cluster]; i++) {
      for(q = 0; q < length_gamma; q++) {
      this_gamma[q] = gamma[i + max_n_y * (*this_cluster + n_clusters * q)];
      qf += prec_delta[q] * this_gamma[q] * this_gamma[q];
      }
      *gamma_ll = 0.5 * (n_y[*this_cluster] * log_det - qf);
      }
      Free(prec_delta);
      Free(this_gamma);
      }
      
      void CopulaLogLikelihood(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *x, double *w2,
      double *gamma, double *eta,
      double *d2, double *log_like_copula, double *ll_sum_copula,
      int *this_cluster) {
      /****************************************************
      # copula likelihood function.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of
      #     thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number
      #     of observations in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses).
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses of
      #     scaling factors.
      #   log_like_copula: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses of
      #     scaling factors.
      #   ll_sum_copula: scalar log-likelihood sum.
      #   this_cluster: cluster indicator.
      # Updates log_like_copula, ll_sum_copula.
      ****************************************************/
      int i, j, k, q;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_timepoints_responses = n_timepoints * n_responses;
      int n_col_z = *n_preds_q * n_responses;
      int length_z = n_timepoints_responses * n_col_z;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_this_x = n_timepoints * *n_preds_q;
      int int_zero = 0;
      int int_one = 1;
      
      double qf;
      double marginal_qf;
      double log_det;
      double *z = (double *)Calloc(length_z, double);
      double *this_mean = (double *)Calloc(n_timepoints_responses, double);
      double *this_e = (double *)Calloc(n_timepoints_responses, double);
      double *this_gamma = (double *)Calloc(n_col_z, double);
      double *this_x = (double *)Calloc(length_this_x, double);
      
      for(i = 0; i < length_z; i++) {
      z[i] = 0.0;
      }
      for(i = 0; i < n_timepoints_responses; i++) {
      this_mean[i] = 0.0;
      this_e[i] = 0.0;
      }
      *ll_sum_copula = 0.0;
      for(i = 0; i < n_y[*this_cluster]; i++) {
      qf = 0.0;
      marginal_qf = 0.0;
      log_det = 0.0;
      for(j = 0; j < n_timepoints; j++) {
      for(q = 0; q < *n_preds_q; q++) {
      this_x[j + n_timepoints * q] =
      x[i + max_n_y * q +
      max_n_y * n_preds * (j + n_timepoints * *this_cluster)];
      }
      }
      MakeZ(dim_thetastar, n_preds_q, this_x, z);
      for(j = 0; j < n_col_z; j++) {
      this_gamma[j] = gamma[i + max_n_y * (*this_cluster + n_clusters * j)];
      }
      matrix_multiply(&int_zero, &int_zero, &n_timepoints_responses,
      &n_col_z, &n_col_z, &int_one, z, this_gamma, this_mean);
      for(j = 0; j < n_timepoints_responses; j++) {
      k = i + max_n_y * (*this_cluster + n_clusters * j);
      this_mean[j] += eta[k];
      this_mean[j] *= d2[k];
      this_e[j] = w2[k] - this_mean[j];
      qf += this_e[j] * this_e[j] / (d2[k] * d2[k]);
      marginal_qf += w2[k] * w2[k];
      log_det += log(d2[k]);
      }
      log_det *= -2.0;
      log_like_copula[i + max_n_y * *this_cluster] =
      0.5 * (log_det - qf + marginal_qf);
      *ll_sum_copula += log_like_copula[i + max_n_y * *this_cluster];
      }
      Free(z);
      Free(this_mean);
      Free(this_e);
      Free(this_gamma);
      Free(this_x);
      }
      
      void UpdateDelta(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *x, double *w2,
      double *gamma, double *delta, double *eta, double *lambda, double *alpha,
      double *distance_matrix_j,
      double *zeta_a, double *zeta_b,
      double *sd_delta, int *acc_delta,
      double *d2, double *log_like_copula, double *ll_sum_copula) {
      /********************************************************************
      # Updates delta.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension
      #     of thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of
      #     observations in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #      of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds.q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   zeta_a: array of dimension
      #     n_clusters x (n_preds_q * n_responses)
      #     of prior shape parameters.
      #   zeta_b: array of dimension
      #     n_clusters x (n_preds_q * n_responses)
      #     of prior scale parameters.
      #   sd_delta: n_clusters x n_responses x n_preds_q array of
      #     candidate standard deviations.
      #   acc_delta: n_clusters x n_responses x n_preds_q array of
      #     acceptance counts.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      #   log_like_copula: array of dimension
      #      max(n_y) x n_clusters of copula likelihoods.
      #     of scaling factors.
      #   ll_sum_copula: vector of length n_clusters of log-likelihood sums.
      # Updates delta.
      ********************************************************************/
      int i, g, j, k, l, m, q;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int n_timepoints_responses = n_timepoints * n_responses;
      int n_responses_preds_q = n_responses * *n_preds_q;
      int length_delta = n_clusters * n_responses_preds_q;
      int length_d2 = max_n_y * n_timepoints * n_clusters * n_responses;
      int length_llc = max_n_y * n_clusters;
      int int_one = 1;
      
      double log_delta, can_log_delta, lp_delta, can_lp_delta, ff, can_ff,
      delta_post, can_delta_post;
      double e = 0.0;
      double z = 0.0;
      double ll_delta = 0.0;
      double can_ll_delta = 0.0;
      double double_zero = 0.0;
      double double_one = 1.0;
      
      double *can_delta = (double *)Calloc(length_delta, double);
      double *can_d2 = (double *)Calloc(length_d2, double);
      double *can_log_like_copula = (double *)Calloc(length_llc, double);
      double *can_ll_sum_copula = (double *)Calloc(n_clusters, double);
      
      for(i = 0; i < length_delta; i++) {
      can_delta[i] = delta[i];
      }
      for(i = 0; i < length_d2; i++) {
      can_d2[i] = d2[i];
      }
      for(i = 0; i < length_llc; i++) {
      can_log_like_copula[i] = log_like_copula[i];
      }
      for(i = 0; i < n_clusters; i++) {
      can_ll_sum_copula[i] = ll_sum_copula[i];
      }
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      for(q = 0; q < n_responses_preds_q; q++) {
      k = g + n_clusters * q;
      log_delta = log(delta[k]);
      RandomNormal(&int_one, &int_one, &double_zero, &double_one, &z);
      can_log_delta = log_delta + sd_delta[k] * z;
      can_delta[k] = exp(can_log_delta);
      lp_delta = (zeta_a[k] - 1.0) * log(delta[k]) -
      delta[k] / zeta_b[k];
      can_lp_delta = (zeta_a[k] - 1.0) * log(can_delta[k]) -
      can_delta[k] / zeta_b[k];
      ff = 1.0 / delta[k];
      can_ff = 1.0 / can_delta[k];
      GammaLogLikelihood(dim_thetastar, n_preds_q, n_y, gamma, delta,
      &g, &ll_delta);
      GammaLogLikelihood(dim_thetastar, n_preds_q, n_y, gamma, can_delta,
      &g, &can_ll_delta);
      MakeD2(dim_thetastar, n_preds_q, n_y, x, can_delta, lambda, alpha,
      distance_matrix_j, can_d2, &g);
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x, w2, gamma, eta,
      can_d2, can_log_like_copula, &can_ll_sum_copula[g], &g);
      delta_post = lp_delta + ff + ll_delta + ll_sum_copula[g];
      can_delta_post = can_lp_delta + can_ff +
      can_ll_delta + can_ll_sum_copula[g];
      RandomExponential(&double_one, &e);
      if(e > (delta_post - can_delta_post)) {
      delta[k] = can_delta[k];
      acc_delta[k] ++;
      ll_sum_copula[g] = can_ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      log_like_copula[l] = can_log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      d2[m] = can_d2[m];
      }
      }
      } else {
      can_delta[k] = delta[k];
      can_ll_sum_copula[g] = ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      can_log_like_copula[l] = log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      can_d2[m] = d2[m];
      }
      }
      }
      }
      }
      PutRNGstate();
      Free(can_delta);
      Free(can_d2);
      Free(can_log_like_copula);
      Free(can_ll_sum_copula);
      }
      
      void UpdateLambda(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *x, double *w2,
      double *gamma, double *delta, double *eta, double *lambda, double *alpha,
      double *distance_matrix_j,
      double *lambda_0, double *nu_0,
      double *tuning_nu, int *acc_lambda,
      double *d2, double *log_like_copula, double *ll_sum_copula) {
      /********************************************************************
      # Updates lambda.
      # Args:
      #   dim_thetastar: a vector of length 6 containing
      #     the dimension of thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing
      #     the number of observations in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).  Not checked to be equal to u2 where appropriate.
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds.q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   lambda_0: n_clusters x n_responses x n_responses array
      #     of prior variances.
      #   nu_0: n_clusters vector of prior degrees of freedom.
      #   tuning_nu: n_clusters vector of tuning parameters.
      #   acc_lambda: n_clusters vector of acceptance counts.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      #   log_like_copula: array of dimension
      #      max(n_y) x n_clusters of copula likelihoods.
      #   ll_sum_copula: vector of length n_clusters of log-likelihood sums.
      # Updates lambda, acc_lambda, copula_log_like, copula_ll_sum, d2.
      *******************************************************************/
      int h, i, g, j, k, l, m, q;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_responses_2 = n_responses * n_responses;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int n_timepoints_responses = n_timepoints * n_responses;
      int length_lambda = n_clusters * n_responses_2;
      int length_d2 = max_n_y * n_timepoints * n_clusters * n_responses;
      int length_llc = max_n_y * n_clusters;
      
      double lp_lambda, can_lp_lambda,
      lambda_post, can_lambda_post, log_prob_cand_given_current,
      log_prob_current_given_cand;
      double e = 0.0;
      double ll_lambda = 0.0;
      double can_ll_lambda = 0.0;
      double double_one = 1.0;
      double n_responses_double = (double) (n_responses);
      
      double *sai = (double *)Calloc(n_responses_2, double);
      double *this_lambda = (double *)Calloc(n_responses_2, double);
      double *can_this_lambda = (double *)Calloc(n_responses_2, double);
      double *this_prec_lambda = (double *)Calloc(n_responses_2, double);
      double *can_this_prec_lambda = (double *)Calloc(n_responses_2, double);
      double *this_lambda_0 = (double *)Calloc(n_responses_2, double);
      double *dummy_matrix = (double *)Calloc(n_responses_2, double);
      double *can_lambda = (double *)Calloc(length_lambda, double);
      double *can_d2 = (double *)Calloc(length_d2, double);
      double *can_log_like_copula = (double *)Calloc(length_llc, double);
      double *can_ll_sum_copula = (double *)Calloc(n_clusters, double);
      
      for(i = 0; i < n_responses_2; i++) {
      sai[i] = 0.0;
      this_lambda[i] = 0.0;
      can_this_lambda[i] = 0.0;
      this_prec_lambda[i] = 0.0;
      can_this_prec_lambda[i] = 0.0;
      this_lambda_0[i] = 0.0;
      dummy_matrix[i] = 0.0;
      }
      for(i = 0; i < length_lambda; i++) {
      can_lambda[i] = lambda[i];
      }
      for(i = 0; i < length_d2; i++) {
      can_d2[i] = d2[i];
      }
      for(i = 0; i < length_llc; i++) {
      can_log_like_copula[i] = log_like_copula[i];
      }
      for(i = 0; i < n_clusters; i++) {
      can_ll_sum_copula[i] = ll_sum_copula[i];
      }
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      k = g + n_clusters * q;
      for(h = 0; h < n_responses_2; h++) {
      sai[h] = lambda[g + n_clusters * h] *
      (tuning_nu[g] - n_responses_double - 1.0);
      }
      RandomInverseWishart(&n_responses, &tuning_nu[g], sai, can_this_lambda);
      for(h = 0; h < n_responses_2; h++) {
      k = g + n_clusters * h;
      this_lambda[h] = lambda[k];
      can_lambda[k] = can_this_lambda[h];
      this_lambda_0[h] = lambda_0[k];
      }
      LogDInverseWishart(
      &n_responses, &nu_0[g], this_lambda_0, this_lambda, &lp_lambda);
      LogDInverseWishart(
      &n_responses, &nu_0[g], this_lambda_0, can_this_lambda, &can_lp_lambda);
      LogDInverseWishart(
      &n_responses, &tuning_nu[g],
      sai, can_this_lambda, &log_prob_cand_given_current);
      for(h = 0; h < n_responses_2; h++) {
      sai[h] = can_lambda[g + n_clusters * h] *
      (tuning_nu[g] - n_responses_double - 1.0);
      }
      LogDInverseWishart(
      &n_responses, &tuning_nu[g],
      sai, this_lambda, &log_prob_current_given_cand);
      EtaLogLikelihood(dim_thetastar, n_preds_q, n_y, eta, lambda, alpha,
      distance_matrix_j, &g, &ll_lambda);
      EtaLogLikelihood(dim_thetastar, n_preds_q, n_y, eta, can_lambda, alpha,
      distance_matrix_j, &g, &can_ll_lambda);
      MakeD2(dim_thetastar, n_preds_q, n_y, x,
      delta, can_lambda, alpha, distance_matrix_j, can_d2, &g);
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x, w2, gamma, eta,
      can_d2, can_log_like_copula, &can_ll_sum_copula[g], &g);
      lambda_post = lp_lambda + log_prob_cand_given_current +
      ll_lambda + ll_sum_copula[g];
      can_lambda_post = can_lp_lambda + log_prob_current_given_cand +
      can_ll_lambda + can_ll_sum_copula[g];
      RandomExponential(&double_one, &e);
      if(e > (lambda_post - can_lambda_post)) {
      for(h = 0; h < n_responses_2; h++) {
      k = g + n_clusters * h;
      lambda[k] = can_lambda[k];
      }
      acc_lambda[g] ++;
      ll_sum_copula[g] = can_ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      log_like_copula[l] = can_log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      d2[m] = can_d2[m];
      }
      }
      } else {
      for(h = 0; h < n_responses_2; h++) {
      k = g + n_clusters * h;
      can_lambda[k] = lambda[k];
      }
      can_ll_sum_copula[g] = ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      can_log_like_copula[l] = log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      can_d2[m] = d2[m];
      }
      }
      }
      }
      PutRNGstate();
      Free(sai);
      Free(this_lambda);
      Free(can_this_lambda);
      Free(this_prec_lambda);
      Free(can_this_prec_lambda);
      Free(this_lambda_0);
      Free(dummy_matrix);
      Free(can_lambda);
      Free(can_d2);
      Free(can_log_like_copula);
      Free(can_ll_sum_copula);
      }
      
      void UpdateAlpha(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *x, double *w2,
      double *gamma, double *delta, double *eta, double *lambda, double *alpha,
      double *distance_matrix_j,
      double *alpha_a, double *alpha_b,
      double *sd_alpha, int *acc_alpha,
      double *d2, double *log_like_copula, double *ll_sum_copula) {
      /********************************************************************
      # Updates alpha.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension
      #     of thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of
      #     observations in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #      of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds.q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   alpha_a: vector of length n_clusters
      #     of shape1 hyperparameters for alpha.
      #   alpha_b: vector of length n_clusters
      #     of shape2 hyperparameters for alpha.
      #   sd_alpha: n_clusters vector of tuning parameters.
      #   acc_alpha: n_clusters vector of acceptance counts.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      #   log_like_copula: array of dimension
      #      max(n_y) x n_clusters of copula likelihoods.
      #     of scaling factors.
      #   ll_sum_copula: vector of length n_clusters of log-likelihood sums.
      # Updates delta.
      ********************************************************************/
      int i, g, j, l, m;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int n_timepoints_responses = n_timepoints * n_responses;
      int length_d2 = max_n_y * n_timepoints * n_clusters * n_responses;
      int length_llc = max_n_y * n_clusters;
      int int_one = 1;
      
      double alpha_prime, can_alpha_prime, lp_alpha, can_lp_alpha, ff, can_ff,
      ll_alpha, can_ll_alpha, alpha_post, can_alpha_post;
      double e = 0.0;
      double z = 0.0;
      double double_zero = 0.0;
      double double_one = 1.0;
      
      double *can_alpha = (double *)Calloc(n_clusters, double);
      double *can_d2 = (double *)Calloc(length_d2, double);
      double *can_log_like_copula = (double *)Calloc(length_llc, double);
      double *can_ll_sum_copula = (double *)Calloc(n_clusters, double);
      
      for(i = 0; i < n_clusters; i++) {
      can_alpha[i] = alpha[i];
      }
      for(i = 0; i < length_d2; i++) {
      can_d2[i] = d2[i];
      }
      for(i = 0; i < length_llc; i++) {
      can_log_like_copula[i] = log_like_copula[i];
      }
      for(i = 0; i < n_clusters; i++) {
      can_ll_sum_copula[i] = ll_sum_copula[i];
      }
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      alpha_prime = log(alpha[g]) - log(1.0 - alpha[g]);
      RandomNormal(&int_one, &int_one, &double_zero, &double_one, &z);
      can_alpha_prime = alpha_prime + z * sd_alpha[g];
      can_alpha[g] = exp(can_alpha_prime) / (1.0 + exp(can_alpha_prime));
      lp_alpha = (alpha_a[g] - 1.0) * log(alpha[g]) +
      (alpha_b[g] - 1.0) * log(1.0 - alpha[g]);
      can_lp_alpha = (alpha_a[g] - 1.0) * log(can_alpha[g]) +
      (alpha_b[g] - 1.0) * log(1.0 - can_alpha[g]);
      ff = log(alpha[g]) - log(1.0 - alpha[g]);
      can_ff = log(can_alpha[g]) - log(1.0 - can_alpha[g]);
      
      EtaLogLikelihood(dim_thetastar, n_preds_q, n_y, eta, lambda, alpha,
      distance_matrix_j, &g, &ll_alpha);
      EtaLogLikelihood(dim_thetastar, n_preds_q, n_y, eta, lambda, can_alpha,
      distance_matrix_j, &g, &can_ll_alpha);
      MakeD2(dim_thetastar, n_preds_q, n_y, x,
      delta, lambda, can_alpha, distance_matrix_j, can_d2, &g);
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x, w2, gamma, eta,
      can_d2, can_log_like_copula, &can_ll_sum_copula[g], &g);
      alpha_post = lp_alpha + ff + ll_alpha + ll_sum_copula[g];
      can_alpha_post = can_lp_alpha + can_ff +
      can_ll_alpha + can_ll_sum_copula[g];
      RandomExponential(&double_one, &e);
      if(e > (alpha_post - can_alpha_post)) {
      alpha[g] = can_alpha[g];
      acc_alpha[g] ++;
      ll_sum_copula[g] = can_ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      log_like_copula[l] = can_log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      d2[m] = can_d2[m];
      }
      }
      } else {
      can_alpha[g] = alpha[g];
      can_ll_sum_copula[g] = ll_sum_copula[g];
      for(i = 0; i < max_n_y; i++) {
      l = i + max_n_y * g;
      can_log_like_copula[l] = log_like_copula[l];
      for(j = 0; j < n_timepoints_responses; j++) {
      m = i + max_n_y * (g + n_clusters * j);
      can_d2[m] = d2[m];
      }
      }
      }
      }
      PutRNGstate();
      Free(can_alpha);
      Free(can_d2);
      Free(can_log_like_copula);
      Free(can_ll_sum_copula);
      }
      
      void StoreCopulaParameters (
      int *i, int *burn, int *num_keep, int *thin,
      int *dim_thetastar, int *n_preds_q, int *n_y,
      double *gamma, double *delta, double *eta, double *lambda, double *alpha,
      double *gamma_keep, double *delta_keep,
      double *eta_keep,
      double *lambda_keep, double *alpha_keep) {
      /************************************************
      # stores parameters.
      # Args:
      #   i: current iteration.
      #   burn: number of samples to burn.
      #   num_keep: number of samples to keep.
      #   thin: number of samples to thin.
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of observations
      #     in each cluster.
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x (n_responses x n_preds_q) array of variances.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_timepoints x n_responses)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   gamma_keep: an array for storing gamma means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   delta_keep: an array for storing delta of dimension
      #     num_keep x n_clusters x (n_responses x n_preds_q).
      #   eta_keep: an array for storing eta means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   lambda_keep: an array for storing lambda of dimension
      #     num_keep x n_clusters x n_responses x n_responses.
      #   alpha_keep: an array for storing alpha of dimension
      #     num_keep x n_clusters.
      # Updates:
      #    gamma_keep, delta_keep, eta_keep, lambda_keep, alpha_keep.
      ************************************************/
      int k, l;
      int i1 = (*i + 1 - *burn) / *thin - 1;
      int i2 = i1 + 1;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_responses_2 = n_responses * n_responses;
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_gamma = max_n_y * n_clusters * (n_responses * *n_preds_q);
      int length_delta = n_clusters * (n_responses * *n_preds_q);
      int length_eta = max_n_y * n_clusters * (n_timepoints * n_responses);
      int length_lambda = n_clusters * n_responses_2;
      
      if(i2 > 1) {
      for(k = 0; k < length_gamma; k++) {
      l = 1 + 2 * k;
      gamma_keep[l] = (
      (i2 - 2) * gamma_keep[l] / (i2 - 1) +
      pow(gamma[k] - gamma_keep[l - 1], 2) / i2
      );
      }
      for(k = 0; k < length_eta; k++) {
      l = 1 + 2 * k;
      eta_keep[l] = (
      (i2 - 2) * eta_keep[l] / (i2 - 1) +
      pow(eta[k] - eta_keep[l - 1], 2) / i2
      );
      }
      }
      for(k = 0; k < length_gamma; k++) {
      l = 0 + 2 * k;
      gamma_keep[l] = (gamma_keep[l] * (i2 - 1) + gamma[k]) / i2;
      }
      for(k = 0; k < length_eta; k++) {
      l = 0 + 2 * k;
      eta_keep[l] = (eta_keep[l] * (i2 - 1) + eta[k]) / i2;
      }
      for(k = 0; k < length_delta; k++) {
      delta_keep[i1 + *num_keep * k] = delta[k];
      }
      for(k = 0; k < length_lambda; k++) {
      lambda_keep[i1 + *num_keep * k] = lambda[k];
      }
      for(k = 0; k < n_clusters; k++) {
      alpha_keep[i1 + *num_keep * k] = alpha[k];
      }
      }
      
      void RunCopula(int *dim_thetastar, int *n_preds_q, int *n_y,
      double *x, double *w2,
      double *gamma, double *delta, double *eta, double *lambda, double *alpha,
      int *cop,
      double *distance_matrix_j,
      double *zeta_a, double *zeta_b,
      double *sd_delta, int *acc_delta,
      double *lambda_0, double *nu_0,
      double *tuning_nu, int *acc_lambda,
      double *alpha_a, double *alpha_b,
      double *sd_alpha, int *acc_alpha,
      double *d2, double *log_like_copula, double *ll_sum_copula
      ) {
      /********************************************************************
      # Updates the copula prior.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension
      #     of thetastar and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   n_y: vector of length n_clusters containing the number of
      #     observations in each cluster.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #      of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds.q * n_responses)
      #     of random regression parameters.
      #   delta: n_clusters x n_responses x n_preds_q array of variances.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      #   lambda: n_clusters x n_responses x n_responses array of variances.
      #   alpha: n_clusters x 1 array of correlation parameters.
      #   cop: integer indicator for copula type.
      #     0: independent.
      #     1: random ints and slopes.
      #     2: structured dependence across time or responses.
      #     3: types 1 and 2.
      #   distance_matrix_j: n_timepoints x n_timepoints matrix of distances.
      #   zeta_a: array of dimension
      #     n_clusters x (n_preds_q * n_responses)
      #     of prior shape parameters.
      #   zeta_b: array of dimension
      #     n_clusters x (n_preds_q * n_responses)
      #     of prior scale parameters.
      #   sd_delta: n_clusters x n_responses x n_preds_q array of
      #     candidate standard deviations.
      #   acc_delta: n_clusters x n_responses x n_preds_q array of
      #     acceptance counts.
      #   lambda_0: n_clusters x n_responses x n_responses array
      #     of prior variances.
      #   nu_0: n_clusters vector of prior degrees of freedom.
      #   tuning_nu: n_clusters vector of tuning parameters.
      #   acc_lambda: n_clusters vector of acceptance counts.
      #   alpha_a: vector of length n_clusters
      #     of shape1 hyperparameters for alpha.
      #   alpha_b: vector of length n_clusters
      #     of shape2 hyperparameters for alpha.
      #   sd_alpha: n_clusters vector of tuning parameters.
      #   acc_alpha: n_clusters vector of acceptance counts.
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of scaling factors.
      #   log_like_copula: array of dimension
      #      max(n_y) x n_clusters of copula likelihoods.
      #     of scaling factors.
      #   ll_sum_copula: vector of length n_clusters of log-likelihood sums.
      # Updates gamma, delta, eta, lambda, alpha,
      #   acc_delta, acc_lambda, acc_alpha, log_like_copula, ll_sum_copula, d2.
      ********************************************************************/
      int g;
      int n_clusters = dim_thetastar[3];
      int n_timepoints = dim_thetastar[5];
      
      if(*cop == 1 || *cop == 3) {
      UpdateGamma(dim_thetastar, n_preds_q, n_y, x,
      w2, gamma, delta, eta, d2);
      for(g = 0; g < n_clusters; g++) {
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x, w2, gamma, eta,
      d2, log_like_copula, &ll_sum_copula[g], &g);
      }
      UpdateDelta(dim_thetastar, n_preds_q, n_y,
      x, w2,
      gamma, delta, eta, lambda, alpha,
      distance_matrix_j,
      zeta_a, zeta_b,
      sd_delta, acc_delta,
      d2, log_like_copula, ll_sum_copula);
      }
      if(*cop > 1) {
      UpdateEta(dim_thetastar, n_preds_q, n_y, x,
      w2, gamma, eta,
      lambda, alpha, distance_matrix_j, d2);
      for(g = 0; g < n_clusters; g++) {
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x, w2, gamma, eta,
      d2, log_like_copula, &ll_sum_copula[g], &g);
      }
      UpdateLambda(dim_thetastar, n_preds_q, n_y,
      x, w2,
      gamma, delta, eta, lambda, alpha,
      distance_matrix_j,
      lambda_0, nu_0,
      tuning_nu, acc_lambda,
      d2, log_like_copula, ll_sum_copula);
      if(n_timepoints > 1) {
      UpdateAlpha(dim_thetastar, n_preds_q, n_y,
      x, w2,
      gamma, delta, eta, lambda, alpha,
      distance_matrix_j,
      alpha_a, alpha_b,
      sd_alpha, acc_alpha,
      d2, log_like_copula, ll_sum_copula);
      }
      }
      }
      
      void UpdateTuningParms(int *dim_thetastar, int *n_preds_q,
      double *tuning_thetastar, int *acc_thetastar,
      double *sd_rho, int *acc_rho,
      double *tuning_xi_low, int *acc_xi_low,
      double *tuning_xi_high, int *acc_xi_high,
      double *tuning_shape, int *acc_shape,
      double *sd_delta, int *acc_delta,
      double *tuning_nu, int *acc_lambda,
      double *sd_alpha, int *acc_alpha
      ) {
      /*********************************************************
      # Updates the candidate standard deviations and acceptance counts.
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #     and the number of timepoints.
      #   n_preds_q: number of predictors in the copula.
      #   tuning_thetastar: array of dimension dim_thetastar[1:5] of
      #     tuning parameters for thetastar.
      #   acc_thetastar: array of dimension dim_thetastar[1:5] of
      #     acceptance counts for thetastar.
      #   sd_rho: the candidate standard deviation.
      #   acc_rho: counter for the number accepted.
      #   tuning_xi_low: tuning paramters for xi_low.
      #   acc_xi_low: acceptance counts for xi_low.
      #   tuning_xi_high: tuning paramters for xi_high.
      #   acc_xi_high: acceptance counts for xi_high.
      #   tuning_shape: n_responses vector of tuning parameters.
      #   acc_shape: n_responses vector of acceptance counts.
      #   sd_delta: n_clusters x n_responses x n_preds_q array of
      #     candidate standard deviations.
      #   acc_delta: n_clusters x n_responses x n_preds_q array of
      #     acceptance counts.
      #   tuning_nu: n_clusters vector of tuning parameters.
      #   acc_lambda: n_clusters vector of acceptance counts.
      #   sd_alpha: n_clusters vector of tuning parameters.
      #   acc_alpha: n_clusters vector of acceptance counts.
      # Updates: tuning_thetastar, acc_thetastar, sd_rho, acc_rho,
      #   tuning_xi_low, acc_xi_low, tuning_xi_high, acc_xi_high,
      #   tuning_shape, acc_shape, sd_delta, acc_delta,
      #   tuning_nu, acc_lambda, sd_alpha, acc_alpha.
      *********************************************************/
      int i;
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int length_thetastar = dim_thetastar[0];
      for(i = 1; i < 5; i++) {
      length_thetastar *= dim_thetastar[i];
      }
      int length_tail = n_clusters * n_responses;
      int length_delta = n_clusters * n_responses * *n_preds_q;
      double n_responses_double = (double) (n_responses);
      
      for(i = 0; i < length_thetastar; i++) {
      if(acc_thetastar[i] > 50) {
      tuning_thetastar[i] *= 1.2;
      } else if(acc_thetastar[i] < 30) {
      tuning_thetastar[i] *= 0.8;
      }
      acc_thetastar[i] = 0;
      }
      if(*acc_rho > 50) {
      *sd_rho *= 1.2;
      } else if(*acc_rho < 30) {
      *sd_rho *= 0.8;
      }
      *acc_rho = 0;
      for(i = 0; i < length_tail; i++) {
      if(acc_xi_low[i] > 50) {
      tuning_xi_low[i] *= 1.2;
      } else if (acc_xi_low[i] < 30) {
      tuning_xi_low[i] *= 0.8;
      }
      if(acc_xi_high[i] > 50) {
      tuning_xi_high[i] *= 1.2;
      } else if (acc_xi_high[i] < 30) {
      tuning_xi_high[i] *= 0.8;
      }
      acc_xi_low[i] = 0;
      acc_xi_high[i] = 0;
      }
      for(i = 0; i < n_responses; i++) {
      if(acc_shape[i] > 50) {
      tuning_shape[i] *= 1.2;
      } else if (acc_shape[i] < 30) {
      tuning_shape[i] *= 0.8;
      }
      acc_shape[i] = 0;
      }
      for(i = 0; i < length_delta; i++) {
      if(acc_delta[i] > 50) {
      sd_delta[i] *= 1.2;
      } else if(acc_delta[i] < 30) {
      sd_delta[i] *= 0.8;
      }
      acc_delta[i] = 0;
      }
      for(i = 0; i < n_clusters; i++) {
      if(acc_lambda[i] > 50) {
      tuning_nu[i] *= 0.8;
      } else if(acc_lambda[i] < 30) {
      tuning_nu[i] *= 1.2;
      }
      if(tuning_nu[i] <= (n_responses_double + 1.0)) {
      tuning_nu[i] = n_responses_double + 2.0;
      }
      acc_lambda[i] = 0;
      if(acc_alpha[i] > 50) {
      sd_alpha[i] *= 1.2;
      } else if (acc_alpha[i] < 30) {
      sd_alpha[i] *= 0.8;
      }
      acc_alpha[i] = 0;
      }
      }
      
      void CheckCopula (int *mcmc_ints, double *x, double *w2, double *zeta,
      double *lambda_0, double *alpha_array, double *timepoints,
      double *gamma_keep, double *delta_keep,
      double *eta_keep,
      double *lambda_keep, double *alpha_keep
      ) {
      /************************************************
      #   mcmc_ints: a vector containing dim_thetastar, n_preds_q, n_y,
      #     burn, keep, thin, indicator_mu, correlation_type, cop.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #      of predictors.
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).
      #   zeta: 2 x n_clusters x (n_preds_q * n_responses) array of prior shape
      #     and scale parameters for delta.
      #   lambda_0: array of dimension n_clusters x (n_responses * n_responses)
      #     of prior variances.
      #   nu_0: n_clusters length vector of prior df.
      #   alpha_array: 3 x n_clusters array of
      #    shape1 and shape2 (for alpha) and nu_0 (for lambda) hyperparameters.
      #   timepoints: vector of timepoints of length n_timepoints.
      #   gamma_keep: an array for storing gamma means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   delta_keep: an array for storing delta of dimension
      #     num_keep x n_clusters x (n_responses x n_preds_q).
      #   eta_keep: an array for storing eta means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   lambda_keep: an array for storing lambda of dimension
      #     num_keep x n_clusters x n_responses x n_responses.
      #   alpha_keep: an array for storing alpha of dimension
      #     num_keep x n_clusters.
      ************************************************/
      int g, i, k;
      int *dim_thetastar = (int *)Calloc(6, int);
      int length_thetastar = 1;
      for(i = 0; i < 6; i++) {
      dim_thetastar[i] = mcmc_ints[i];
      if(i < 5) {
      length_thetastar *= mcmc_ints[i];
      }
      }
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_responses_2 = n_responses * n_responses;
      int n_timepoints_2 = n_timepoints * n_timepoints;
      int n_preds_q = mcmc_ints[6];
      int length_xi = n_clusters * n_responses;
      int *n_y = (int *)Calloc(n_clusters, int);
      for(g = 0; g < n_clusters; g++) {
      n_y[g] = mcmc_ints[7 + g];
      }
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int burn = mcmc_ints[n_clusters + 7];
      int keep = mcmc_ints[n_clusters + 8];
      int thin = mcmc_ints[n_clusters + 9];
      //  int indicator_mu = mcmc_ints[n_clusters + 10];
      //  int correlation_type = mcmc_ints[n_clusters + 11];
      int cop = mcmc_ints[n_clusters + 12];
      
      int iters = burn + keep;
      int num_keep =  keep / thin;
      int length_gamma = max_n_y * n_clusters * n_preds_q * n_responses;
      int length_delta = n_clusters * (n_preds_q * n_responses);
      int length_eta = max_n_y * n_clusters * n_timepoints * n_responses;
      int length_lambda = n_clusters * n_responses_2;
      int length_log_like_copula = max_n_y * n_clusters;
      
      double delta_constant = 0.0;
      double lambda_constant = 0.0;
      double alpha_constant = 0.0;
      if(cop == 1 || cop == 3) {
      delta_constant = 1.0;
      }
      if(cop > 1) {
      lambda_constant = 1.0;
      }
      if(n_timepoints > 1 && cop > 1) {
      alpha_constant = 0.5;
      }
      double *zeta_a = (double *) Calloc(length_delta, double);
      double *zeta_b = (double *) Calloc(length_delta, double);
      double *alpha_a = (double *) Calloc(n_clusters, double);
      double *alpha_b = (double *) Calloc(n_clusters, double);
      double *nu_0 = (double *) Calloc(n_clusters, double);
      
      int *acc_thetastar = (int *) Calloc(length_thetastar, int);
      int acc_rho = 0;
      int *acc_xi_low = (int *) Calloc(length_xi, int);
      int *acc_xi_high = (int *) Calloc(length_xi, int);
      int *acc_shape = (int *) Calloc(n_responses, int);
      int *acc_delta = (int *) Calloc(length_delta, int);
      int *acc_lambda = (int *) Calloc(n_responses, int);
      int *acc_alpha = (int *) Calloc(n_clusters, int);
      
      double *tuning_thetastar = (double *) Calloc(length_thetastar, double);
      double sd_rho = 1.0;
      double *tuning_xi_low = (double *) Calloc(length_xi, double);
      double *tuning_xi_high = (double *) Calloc(length_xi, double);
      double *tuning_shape = (double *) Calloc(n_responses, double);
      double *sd_delta = (double *) Calloc(length_delta, double);
      double *tuning_nu = (double *) Calloc(n_responses, double);
      double *sd_alpha = (double *) Calloc(n_clusters, double);
      
      double *gamma = (double *) Calloc(length_gamma, double);
      double *delta = (double *) Calloc(length_delta, double);
      double *eta = (double *) Calloc(length_eta, double);
      double *lambda = (double *) Calloc(length_lambda, double);
      double *alpha = (double *) Calloc(n_clusters, double);
      double *distance_matrix_j = (double *) Calloc(n_timepoints_2, double);
      double *d2 = (double *) Calloc(length_eta, double);
      double *log_like_copula = (double *) Calloc(length_log_like_copula, double);
      double *ll_sum_copula = (double *) Calloc(n_clusters, double);
      
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < (n_preds_q * n_responses); i++) {
      k = g + n_clusters * i;
      zeta_a[k] = zeta[0 + 2 * k];
      zeta_b[k] = zeta[1 + 2 * k];
      }
      }
      for(g = 0; g < n_clusters; g++) {
      alpha_a[g] = alpha_array[0 + 3 * g];
      alpha_b[g] = alpha_array[1 + 3 * g];
      nu_0[g] = alpha_array[2 + 3 * g];
      }
      for(i = 0; i < length_thetastar; i++) {
      acc_thetastar[i] = 0;
      tuning_thetastar[i] = 0.1;
      }
      for(i = 0; i < length_xi; i++) {
      acc_xi_low[i] = 0;
      acc_xi_high[i] = 0;
      tuning_xi_low[i] = 0.1;
      tuning_xi_high[i] = 0.1;
      }
      for(i = 0; i < n_responses; i++) {
      acc_shape[i] = 0;
      tuning_shape[i] = 0.1;
      }
      for(i = 0; i < length_delta; i++) {
      acc_delta[i] = 0;
      sd_delta[i] = 0.1;
      }
      for(i = 0; i < n_clusters; i++) {
      acc_lambda[i] = 0;
      tuning_nu[i] = (double) (n_responses) - 1.0 + 10.0;
      }
      for(i = 0; i < n_clusters; i++) {
      acc_alpha[i] = 0;
      sd_alpha[i] = 0.1;
      }
      
      for(i = 0; i < length_gamma; i++) {
      gamma[i] = 0.0;
      }
      for(i = 0; i < length_delta; i++) {
      delta[i] = delta_constant;
      }
      for(i = 0; i < length_eta; i++) {
      eta[i] = 0.0;
      }
      for(i = 0; i < length_lambda; i++) {
      lambda[i] = 0.0;
      }
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < n_responses; i++) {
      lambda[g + n_clusters * (i + n_responses * i)] = lambda_constant;
      }
      }
      for(i = 0; i < n_clusters; i++) {
      alpha[i] = alpha_constant;
      }
      for(i = 0; i < n_timepoints_2; i++) {
      distance_matrix_j[i] = 0.0;
      }
      MakeAbsoluteDistance(&n_timepoints, timepoints, distance_matrix_j);
      
      for(i = 0; i < length_eta; i++) {
      d2[i] = 0.0;
      }
      for(i = 0; i < length_log_like_copula; i++) {
      log_like_copula[i] = 0.0;
      }
      for(i = 0; i < n_clusters; i++) {
      ll_sum_copula[i] = 0.0;
      }
      
      for(g = 0; g < n_clusters; g++) {
      MakeD2(dim_thetastar, &n_preds_q, n_y, x,
      delta, lambda, alpha, distance_matrix_j, d2, &g);
      CopulaLogLikelihood(dim_thetastar, &n_preds_q, n_y, x, w2,
      gamma, eta, d2, log_like_copula, &ll_sum_copula[g], &g);
      }
      
      
      
      for(i = 0; i < iters; i++) {
      if( (i + 1) % 100 == 0) {
      Rprintf("iteration %d\n", i);
      }
      RunCopula(dim_thetastar, &n_preds_q, n_y,
      x, w2, gamma, delta, eta, lambda, alpha,
      &cop, distance_matrix_j,
      zeta_a, zeta_b,
      sd_delta, acc_delta,
      lambda_0, nu_0,
      tuning_nu, acc_lambda,
      alpha_a, alpha_b,
      sd_alpha, acc_alpha,
      d2, log_like_copula, ll_sum_copula);
      if ( ((i + 1) % 100 == 0) && (i < burn) ) {
      UpdateTuningParms(dim_thetastar, &n_preds_q,
      tuning_thetastar, acc_thetastar,
      &sd_rho, &acc_rho,
      tuning_xi_low, acc_xi_low,
      tuning_xi_high, acc_xi_high,
      tuning_shape, acc_shape,
      sd_delta, acc_delta,
      tuning_nu, acc_lambda,
      sd_alpha, acc_alpha);
      }
      
      if(i >= burn) {
      if((i + 1 - burn) % thin == 0) {
      StoreCopulaParameters(&i, &burn, &num_keep, &thin,
      dim_thetastar,
      &n_preds_q, n_y,
      gamma, delta, eta, lambda, alpha,
      gamma_keep, delta_keep, eta_keep, lambda_keep, alpha_keep);
      }
      }
      }
      
      for(g = 0; g < n_clusters; g++) {
      Rprintf("tuning_nu[g] = %f\n", tuning_nu[g]);
      }
      Free(zeta_a);
      Free(zeta_b);
      Free(alpha_a);
      Free(alpha_b);
      Free(nu_0);
      Free(acc_thetastar);
      
      Free(acc_xi_low);
      Free(acc_xi_high);
      Free(acc_shape);
      Free(acc_delta);
      Free(acc_lambda);
      Free(acc_alpha);
      Free(tuning_thetastar);
      Free(tuning_xi_low);
      Free(tuning_xi_high);
      Free(tuning_shape);
      Free(sd_delta);
      Free(tuning_nu);
      Free(sd_alpha);
      
      Free(gamma);
      Free(delta);
      Free(eta);
      Free(lambda);
      Free(alpha);
      Free(distance_matrix_j);
      Free(d2);
      Free(sd_alpha);
      Free(log_like_copula);
      Free(ll_sum_copula);
      
      }
      
      void UpdateThetastar(int *dim_thetastar, int *n_y, int *basis_ind,
      double *basis_knots, double *y, double *y_low, double *y_high,
      int *status, double *x,
      double *u, double *u_low, double *u_high,
      double *can_u, double *can_u_low, double *can_u_high,
      int *bin, int *bin_low, int *bin_high,
      double *log_like, double *ll_sum,
      double *can_log_like, double *can_ll_sum,
      double *lp_thetastar, double *log_dets, double *mu_array,
      double *e_thetastar, double *prec_thetastar,
      double *thetastar, double *theta,
      double *can_thetastar, double *can_theta,
      double *shape,
      double *tuning_thetastar, int *acc_thetastar, int *att_thetastar,
      double *tau_low, double *tau_high,
      int *pareto_low, int *pareto_high,  #   nu_0: n_clusters length vector of prior df.
      
      double *xi_low, double *xi_high,
      int *cop, int *n_preds_q, double *timepoints,
      double *u2, double *w, double *w2, double *d2,
      double *log_like_copula, double *ll_sum_copula,
      double *can_u2, double *can_w, double *can_w2,
      double *can_log_like_copula, double *can_ll_sum_copula,
      double *gamma, double *eta) {
      /**************************************************************
      # updates thetastar.
      # Args:
      #   dim_thetastar: a vector of length 6 containing the dimension of thetastar
      #    and the number of timepoints.
      #   n_y: vector of length n_clusters
      #      containing the number of observations in each cluster.
      #  basis_ind: integer indicator where (-1, 0, 1, 2, 3, 4, 5) maps to
      #      ("spline", "Gaussian", "t", "logistic", "ALAP", "Weibull", "gamma").
      #  basis_knots: vector of knots.
      #   y: array of dimension
      #      max(n_y) x n_timepoints x n_clusters x n_responses of responses.
      #   y_low: array of dimension
      #      max(n_y) x n_timepoints x n_clusters x n_responses
      #      of left-censoring points.
      #   y_high: array of dimension
      #      max(n_y) x n_timepoints x n_clusters x n_responses
      #      of right-censoring points.
      #   status: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of status indicators where:
      #     -99 indicates that the value should be ignored;
      #     0 indicates that the value is missing and treated as continuous;
      #     1 indicates that the value is observed;
      #     2 indicates that the value is right-censored;
      #     3 indicates that the value is left-censored;
      #     4 indicates that the value is interval-censored.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   u: array of dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables.
      #   u_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_low.
      #   u_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_high.
      #   can_u: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables.  Should equal u.
      #   can_u_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_low.  Should equal u_low.
      #   can_u_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of uniform random variables for y_high. Should equal u_high.
      #   bin: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u is in.
      #   bin_low: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u_low is in.
      #   bin_high: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of bins indicating which interval u_high is in.
      #   log_like: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of loglikelihood (or censored loglikelihood) values.
      #   can_log_like: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of loglikelihood (or censored loglikelihood) values.
      #     Should equal log_like.
      #   ll_sum: scalar summed loglikelihood (or censored loglikelihood) values.
      #   can_ll_sum: scalar summed loglikelihood
      #     (or censored loglikelihood) values.  Should equal ll_sum.
      #   lp_thetastar: log prior of thetastar.
      #   log_dets: a vector of length 5 containing the log determinants of
      #      prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   mu_array: dim.thetastar array of means for thetastar.
      #   e_thetastar: dim.thetastar array of residuals of thetastar.
      #   prec_thetastar: the precision matrix of thetastar
      #     (not checked against all the composite precisions).
      #   thetastar: dim.thetastar[1:5] array of regression parameters.
      #   can_thetastar: dim.thetastar[1:5] array of regression parameters.
      #     Should equal thetastar.
      #   theta: n_bf x n_preds x n_timepoints x n_clusters x n_repsonses array
      #    of thresholded regression parameters.
      #    theta is not checked to be in the space of valid regression parameters.
      #   can_theta: n_bf x n_preds x n_timepoints x n_clusters x n_repsonses array
      #    of thresholded regression parameters.
      #    can_theta is not checked to be in the space
      #    of valid regression parameters.  Should equal theta.
      #   shape: n_responses length vector of shape parameters for bases.
      #   tuning_thetastar: array of dimension dim_thetastar[1:5] of
      #     tuning parameters for thetastar.
      #   acc_thetastar: array of dimension dim_thetastar[1:5] of
      #     acceptance counts for thetastar.
      #   att_thetastar: array of dimension dim_thetastar[1:5] of
      #     attempts counts for thetastar.
      #   tau_low: the quantile level below which a parametric tail is fit.
      #   tau_high: the quantile level above which a parametric tail is fit.
      #   pareto_low: indicator if pareto tail is fit below tau_low.
      #   pareto_high: indicator if pareto tail is fit above tau_high.
      #   xi_low: array of dimension
      #     n_clusters x n_responses of lower tail parameters.
      #   xi.high: array of dimension
      #     n_clusters x n_responses of upper tail parameters.
      #   log_dets: a vector of length 5 containing the log determinants of
      #      prec_bf, prec_preds, prec_time, prec_clusters, and prec_responses.
      #   mu_array: dim_thetastar array of means for thetastar.
      #   e_thetastar: dim_thetastar array of residuals of thetastar.
      #   prec_thetastar: the precision matrix of thetastar
      #     (not checked against all the composite precisions).
      #   cop: integer indicator for copula type.
      #       0: independent.
      #       1: random ints and slopes.
      #       2: structured dependence across time or responses.
      #       3: types 1 and 2.
      #   n_preds_q: number of predictors in the copula.
      #   timepoints: vector of longitudinal timepoints.
      #   u2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of uniform random variables.  Used only for the copula.
      #   w: array of
      #     dimension max(n_y) x n_timepoints x n_clusters x n_responses
      #     of qnorm(u).
      #   w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(u2).
      #   d2: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses of scaling factors.
      #   log_like_copula: array of dimension
      #      max(n_y) x n_clusters of copula likelihoods.
      #   ll_sum_copula: vector of length n_clusters of log-likelihood sums.
      #   can_u2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of uniform random variables.  Used only for the copula.
      #   can_w: array of
      #     dimension max(n_y) x n_timepoints x n_clusters  x n_responses
      #     of qnorm(can_u).
      #   can_w2: array of
      #     dimension max(n_y) x n_clusters x (n_timepoints * n_responses)
      #     of qnorm(can_u2).
      #   can_log_like_copula: array of dimension
      #      max(n_y) x n.clusters of copula candidate likelihoods.
      #   can_ll_sum_copula: vector of length n_clusters
      #      of candidate log-likelihood sums.
      #   gamma: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q * n_responses)
      #     of random regression parameters.
      #   eta: array of dimension
      #     max(n_y) x n_clusters x (n_preds_q x n_timepoints)
      #     of random effects.
      # Updates att_thetastar, acc_thetastar, lp_thetastar, e_thetastar,
      #   thetastar, theta, can_thetastar, can_theta,
      #   u, u_low, u_high, can_u, can_u_low, can_u_high,
      #   log_like, ll_sum, can_log_like, can_ll_sum,
      #   bin, bin_low, bin_high,
      #   w, w2, can_w, can_w2, u2, can_u2,
      #   log_like_copula, ll_sum_copula, can_log_like_copula, can_ll_sum_copula.
      **************************************************************/
      int d, g, h, i, j, k, l, m, p;
      int tail_index, time_index;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_timepoints_responses = n_timepoints * n_responses;
      int dmp = n_bf * n_preds * dim_time;
      int max_n_y = 0;
      int int_one = 1;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      
      int length_thetastar = 1;
      for(i = 0; i < 5; i++) {
      length_thetastar *= dim_thetastar[i];
      }
      
      double ll_total, can_ll_total;
      double can_lp_thetastar = 0.0;
      double e = 0.0;
      double z = 0.0;
      double double_zero = 0.0;
      double double_one = 1.0;
      double *can_e_thetastar = (double *) Calloc(length_thetastar, double);
      for(i = 0; i < length_thetastar; i++) {
      can_e_thetastar[i] = e_thetastar[i];
      }
      
      GetRNGstate();
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      tail_index = g + n_clusters * h;
      for(d = 0; d < dim_time; d++) {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      k = m + n_bf * p + n_bf * n_preds * d + dmp * (g + n_clusters * h);
      RandomNormal(&int_one, &int_one,
      &double_zero, &double_one, &z);
      att_thetastar[k] ++;
      can_thetastar[k] = thetastar[k] + tuning_thetastar[k] * z;
      ThetastarToTheta(dim_thetastar,
      can_thetastar, can_theta, timepoints, &g, &h);
      can_e_thetastar[k] = can_thetastar[k] - mu_array[k];
      LogPriorThetastar(dim_thetastar, log_dets,
      can_e_thetastar, prec_thetastar, &can_lp_thetastar);
      ll_total = 0.0;
      can_ll_total = 0.0;
      for(j = 0; j < n_timepoints; j++) {
      time_index = j + n_timepoints * tail_index;
      Likelihood(dim_thetastar, n_y, &basis_ind[h], basis_knots,
      y, y_low, y_high, status, x, can_u, can_u_low, can_u_high,
      bin, bin_low, bin_high, can_log_like, &can_ll_sum[time_index],
      can_theta, &shape[h], tau_low, tau_high,
      pareto_low, pareto_high,
      &xi_low[tail_index], &xi_high[tail_index], &g, &h);
      ll_total += ll_sum[time_index];
      can_ll_total += can_ll_sum[time_index];
      }
      
      Rprintf("ll_total, lp_thetastar = %f %f \n", ll_total, *lp_thetastar);
      Rprintf("can_ll_total, can_lp_thetastar = %f %f \n", can_ll_total, can_lp_thetastar);
      
      if(*cop > 0) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      l = i + max_n_y * j + max_n_y * n_timepoints * tail_index;
      QNorm(&can_u[l],
      &double_zero, & double_one, &double_zero, &can_w[l]);
      }
      }
      CopyUniform(dim_thetastar, n_y,
      can_u, can_u2, can_w, can_w2, &g, &int_one);
      CopulaLogLikelihood(dim_thetastar, n_preds_q, n_y, x,
      can_w2, gamma, eta, d2,
      can_log_like_copula, can_ll_sum_copula, &g);
      ll_total += ll_sum_copula[g];
      can_ll_total += can_ll_sum_copula[g];
      }
      RandomExponential(&double_one, &e);
      if(e > ll_total + *lp_thetastar - can_ll_total - can_lp_thetastar) {
      acc_thetastar[k] ++;
      e_thetastar[k] = can_e_thetastar[k];
      thetastar[k] = can_thetastar[k];
      theta[k] = can_theta[k];
      for(j = 0; j < n_timepoints; j++) {
      l = j + n_timepoints * tail_index;
      ll_sum[l] = can_ll_sum[l];
      }
      if(*cop > 0) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      l = i + max_n_y * j + max_n_y * n_timepoints * tail_index;
      w[l] = can_w[l];
      }
      for(j = 0; j < n_timepoints_responses; j++) {
      l = i + max_n_y * g + max_n_y * n_clusters * j;
      u2[l] = can_u2[l];
      w2[l] = can_w2[l];
      }
      l = i + max_n_y * g;
      log_like_copula[i] = can_log_like_copula[i];
      }
      ll_sum_copula[g] = can_ll_sum_copula[g];
      }
      } else {
      can_e_thetastar[k] = e_thetastar[k];
      can_thetastar[k] = thetastar[k];
      can_theta[k] = theta[k];
      for(j = 0; j < n_timepoints; j++) {
      l = j + n_timepoints * tail_index;
      can_ll_sum[l] = ll_sum[l];
      }
      if(*cop > 0) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      l = i + max_n_y * j + max_n_y * n_timepoints * tail_index;
      can_w[l] = w[l];
      }
      for(j = 0; j < n_timepoints_responses; j++) {
      l = i + max_n_y * g + max_n_y * n_clusters * j;
      can_u2[l] = u2[l];
      can_w2[l] = w2[l];
      }
      l = i + max_n_y * g;
      can_log_like_copula[i] = log_like_copula[i];
      }
      can_ll_sum_copula[g] = ll_sum_copula[g];
      }
      }
      }
      }
      }
      }
      }
      PutRNGstate();
      Free(can_e_thetastar);
      }
      
      void MakeArray(int *dim_thetastar, int *n_y, int *basis_ind,
      double *y_mat, double *y_low_mat, double *y_high_mat,
      double *x_mat, int *status_mat,
      double *y, double *y_low, double *y_high,
      double *x, int *status) {
      /******************************************************
      # Turns matrix response into array.
      # Let dim_y = c(max(n_y) x n_timepoints x n_clusters x n_responses).
      # Let dim_x = c(max(n_y) x n_preds x n_timepoints x n_responses).
      # Args:
      #   y_mat: (max_n_y * n_timepoints) x n_responses matrix
      #     of lower endpoints for censored values.
      #   y_low_mat: (max_n_y * n_timepoints) x n_responses matrix
      #     of lower endpoints for censored values.
      #   y_high_mat: (max_n_y * n_timepoints) x n_responses matrix
      #     of upper endpoints for censored values.
      #   x_mat: (max_n_y * n_timepoints) x n_preds matrix
      #     of predictors including intercept.
      #   status_mat: (max_n_y * n_timepoints) x n_responses matrix of statuses.
      #   y: array of y_mat of dimension dim_y.
      #   y_low: array of y_low_mat of dimension dim_y.
      #   y_high: array of y_high_mat of dimension dim_y.
      #   x: array of x_mat of dimension dim_x.
      #   status: array of status of dimension dim_y.
      # Updates y, y_low, y_high, x, status.
      ******************************************************/
      int g, h, i, j, k, l, p;
      int n_preds = dim_thetastar[1];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int big_n_timepoints = max_n_y * n_timepoints;
      int big_n_preds = max_n_y * n_preds;
      int total_n = n_y[0];
      if(n_clusters > 1) {
      for(g = 1; g < n_clusters; g++) {
      total_n += n_y[g];
      }
      }
      total_n *= n_timepoints;
      int this_row = 0;
      for(g = 0; g < n_clusters; g++) {
      for(i = 0; i < n_y[g]; i++) {
      for(j = 0; j < n_timepoints; j++) {
      for(h = 0; h < n_responses; h++) {
      k = i + max_n_y * j + big_n_timepoints * (g + n_clusters * h);
      l = this_row + total_n * h;
      y[k] = y_mat[l];
      y_low[k] = y_low_mat[l];
      y_high[k] = y_high_mat[l];
      status[k] = status_mat[l];
      }
      for(p = 0; p < n_preds; p++) {
      x[i + max_n_y * p + big_n_preds *(j + n_timepoints * g)] =
      x_mat[this_row + total_n * p];
      }
      this_row ++;
      }
      }
      }
      }
      
      void RunMCMC(int *mcmc_ints, double *mcmc_doubles, int *n_y,
      double *prec_bf, double *prec_preds, double *prec_time,
      double *mu_0, double *prec_mu_0,
      double *nu_0, double *sigma_0,
      double *distance_matrix,
      int *basis_ind, int *basis_knots,
      double *y_array, int *status,
      double *x,
      double *timepoints,
      double *zeta,
      double *lambda_0, double *nu_0_copula, double *alpha_array,
      double *mu_keep, double *cov_responses_keep, double *rho_keep,
      double *thetastar_keep, double *theta_keep,
      double *gamma_keep, double *delta_keep,
      double *eta_keep,
      double *lambda_keep, double *alpha_keep
      double *xi_low_keep, double *xi_high_keep,
      double *mu_xi_low_keep, double *mu_xi_high_keep,
      double *tau2_xi_low_keep, double *tau2_xi_high_keep
      ) {
      /*********************************************************
      #   mcmc_ints: a vector containing dim_thetastar, n_preds_q,
      #     how many samples to burn, keep and thin, indicator_mu,
      #     correlation_type, and cop.
      #   mcmc_doubles: a vector where:
      #     1) the first n_preds elements are the locations of the
      #        untransformed predictors.
      #     2) the next n_preds elements are the scales of the
      #        untransformed predictors.
      #     3) the next n_responses elements are the locations of the
      #        untransformed responses.
      #     4) the next n_responses elements are the scales of the
      #        untransformed responses.
      #     5) the next 3 values are max_rho, shape1, shape2.
      #     6) the next 10 values are tau_low, tau_high, a_xi_low, a_xi_high,
      #        b_xi_low, b_xi_high, kappa_xi_low, kappa_xi_high,
      #        pareto_probs_low, pareto_probs_high.
      #     7) the next n_responses values are mu_0_xi_low.
      #     8) the next n_responses values are mu_0_xi_high.
      #   n_y: a vector of length n_clusters.
      #   prec_bf: The precision of the basis functions.
      #   prec_preds: The precision across predictors.
      #   prec_time: The precision across time.
      #   mu_0: the prior mean for mu.
      #   prec_mu_0: the prior precision for mu.
      #   nu_0: the prior scale for cov_responses.
      #   sigma_0: the prior location for cov_responses.
      #   distance_matrix: a valid distance matrix of dimension
      #     n_clusters x n_clusters.
      #   basis_ind: integer indicator where (-1, 0, 1, 2, 3, 4, 5) maps to
      #      ("spline", "Gaussian", "t", "logistic", "ALAP", "Weibull", "gamma").
      #   basis_knots: vector of knots.
      #   y_array : array of dimension
      #     7 x max(n_y) x n_timepoints x n_clusters x n_responses of y, y_low,
      #     y_high, y_mean, y_var, u_mean, u_var.
      #   status: array of dimension
      #     max(n_y) x n_timepoints x n_clusters x n_responses
      #     of status indicators where:
      #     -99 indicates that the value should be ignored;
      #     0 indicates that the value is missing and treated as continuous;
      #     1 indicates that the value is observed;
      #     2 indicates that the value is right-censored;
      #     3 indicates that the value is left-censored;
      #     4 indicates that the value is interval-censored.
      #   x: array of dimension max(n_y) x n_preds x n_timepoints x n_clusters
      #     of predictors.
      #   timepoints: vector of timepoints of length n_timepoints.
      #   zeta: 2 x n_clusters x (n_preds_q * n_responses) array of prior shape
      #     and scale parameters for delta.
      #   lambda_0: array of dimension n_clusters x (n_responses * n_responses)
      #     of prior variances.
      #   nu_0_copula: n_clusters length vector of prior df.
      #   alpha_array: 3 x n_clusters array of
      #    shape1 and shape2 (for alpha) and nu_0 (for lambda) hyperparameters.
      #   mu_keep: an array for storing mu of dimension
      #     num_keep x n_bf x num_preds x dim_time
      #     or num_keep x 2 x num_preds x dim_time if indicator_mu = 1.
      #   cov_responses_keep: an array for storing cov_responses of dimension
      #     num_keep x n_responses x n_responses.
      #   rho_keep: an array for storing rho of dimension
      #     num_keep x 1.
      #   thetastar_keep: array for storing thetastar of dimension
      #     num_keep x dim_thetastar.
      #   theta_keep: array for storing thetastar of dimension
      #     num_keep x dim_theta.
      #   gamma_keep: an array for storing gamma means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   delta_keep: an array for storing delta of dimension
      #     num_keep x n_clusters x (n_responses x n_preds_q).
      #   eta_keep: an array for storing eta means and variances of dimension
      #     2 x max(n_y) x n_clusters x (n_preds_q * n_responses).
      #   lambda_keep: an array for storing lambda of dimension
      #     num_keep x n_clusters x n_responses x n_responses.
      #   alpha_keep: an array for storing alpha of dimension
      #     num_keep x n_clusters.
      #   xi_low_keep: array for storing xi_low of dimension
      #     num_keep x n_clusters x n_responses.
      #   xi_high_keep: array for storing xi_high of dimension
      #     num_keep x n_clusters x n_responses.
      #   mu_xi_low_keep: array for storing mu_xi_low of dimension
      #     num_keep x n_responses.
      #   mu_xi_high_keep: array for storing mu_xi_high of dimension
      #     num_keep x n_responses.
      #   tau2_xi_low_keep: array for storing tau2_xi_low of dimension num_keep x 1.
      #   tau2_xi_high_keep: array for storing tau2_xi_high of dimension
      #     num_keep x 1.
      # Updates:
      #    mu_keep, cov_responses_keep, rho_keep, thetastar_keep, theta_keep,
      #    xi_low_keep, xi_high_keep, mu_xi_low_keep, mu_xi_high_keep,
      #    tau2_xi_low_keep, tau2_xi_high_keep.
      ************************************************/
      int g, h, i, k, p;
      int *dim_thetastar = (int *)Calloc(6, int);
      int length_thetastar = 1;
      for(i = 0; i < 6; i++) {
      dim_thetastar[i] = mcmc_ints[i];
      if(i < 5) {
      length_thetastar *= mcmc_ints[i];
      }
      }
      int length_thetastar_2 = length_thetastar * length_thetastar;
      int n_bf = dim_thetastar[0];
      int n_preds = dim_thetastar[1];
      int dim_time = dim_thetastar[2];
      int n_clusters = dim_thetastar[3];
      int n_responses = dim_thetastar[4];
      int n_timepoints = dim_thetastar[5];
      int n_preds_q = mcmc_ints[6];
      int n_clusters_2 = n_clusters * n_clusters;
      int n_responses_2 = n_responses * n_responses;
      int length_mu = n_bf * n_preds * dim_time;
      int length_theta = n_bf * n_preds * n_timepoints * n_clusters * n_responses;
      int length_xi = n_clusters * n_responses;
      
      
      for(g = 0; g < n_clusters; g++) {
      n_y[g] = mcmc_ints[7 + g];
      }
      k = 7 + n_clusters;
      int burn = mcmc_ints[k];
      int keep = mcmc_ints[k + 1];
      int thin = mcmc_ints[k + 2];
      int indicator_mu = mcmc_ints[k + 3];
      int correlation_type = mcmc_ints[k + 4];
      int cop = mcmc_ints[k + 5];
      
      for(p = 0; p < n_preds; p++) {
      pred_mns[p] = mcmc_doubles[p];
      pred_mns[p] = mcmc_doubles[n_preds + p];
      }
      k = 2 * n_preds;
      for(h = 0; h < n_responses; h++) {
      response_mns[h] = mcmc_doubles[k + h];
      response_sds[h] = mcmc_doubles[k + n_responses + h];
      }
      k = 2 * (n_preds + n_responses);
      
      double max_rho = mcmc_doubles[k];
      double shape1 = mcmc_doubles[k + 1];
      double shape2 = mcmc_doubles[k + 2];
      double tau_low = mcmc_doubles[k + 3];
      double tau_high = mcmc_doubles[k + 4];
      double a_xi_low = mcmc_doubles[k + 5];
      double a_xi_high = mcmc_doubles[k + 6];
      double b_xi_low = mcmc_doubles[k + 7];
      double b_xi_high = mcmc_doubles[k + 8];
      double kappa_xi_low = mcmc_doubles[k + 9];
      double kappa_xi_high = mcmc_doubles[k + 10];
      for(h = 0; h < n_responses; h++) {
      mu_0_xi_low[h] = mcmc_doubles[k + 11 + h];
      mu_0_xi_high[h] = mcmc_doubles[k + 11 + n_responses + h];
      }
      double pareto_prob_low = mcmc_doubles[k + 11 + 2 * n_responses];
      double pareto_prob_high = mcmc_doubles[k + 11 + 2 * n_responses + 1];
      
      int length_ll_sum = n_timepoints * n_clusters * n_responses;
      int length_thetastar_2 = length_thetastar * length_thetastar;
      int n_clusters_2 = n_clusters * n_clusters;
      int n_responses_2 = n_responses * n_responses;
      
      int max_n_y = 0;
      find_integer_max(&n_clusters, n_y, &max_n_y);
      int length_y = max_n_y * n_timepoints * n_clusters * n_responses;
      
      if(indicator_mu == 1) {
      length_mu = 2 * n_preds * dim_time;
      }
      double *log_dets        = (double *)Calloc(5, double);
      double *thetastar       = (double *)Calloc(length_thetastar, double);
      double *e_thetastar     = (double *)Calloc(length_thetastar, double);
      double *mu_array        = (double *)Calloc(length_thetastar, double);
      double *prec_thetastar  = (double *)Calloc(length_thetastar_2, double);
      double *mu              = (double *)Calloc(length_mu, double);
      double *mu_array        = (double *)Calloc(length_thetastar, double);
      
      double *cov_clusters    = (double *)Calloc(n_clusters_2, double);
      double *prec_clusters   = (double *)Calloc(n_clusters_2, double);
      double *cov_responses   = (double *)Calloc(n_responses_2, double);
      double *prec_responses  = (double *)Calloc(n_responses_2, double);
      
      double *xi_low          = (double *)Calloc(length_xi, double);
      double *xi_high         = (double *)Calloc(length_xi, double);
      double *e_xi_low        = (double *)Calloc(length_xi, double);
      double *e_xi_high       = (double *)Calloc(length_xi, double);
      double *mu_0_xi_low     = (double *)Calloc(n_respones, double);
      double *mu_0_xi_high    = (double *)Calloc(n_respones, double);
      double *mu_xi_low       = (double *)Calloc(n_respones, double);
      double *mu_xi_high      = (double *)Calloc(n_respones, double);
      
      int *n_y = (int *)Calloc(n_clusters, int);
      double *pred_mns = (double *) Calloc(n_preds, double);
      double *pred_sds = (double *) Calloc(n_preds, double);
      double *response_mns = (double *) Calloc(n_responses, double);
      double *response_sds = (double *) Calloc(n_responses, double);
      
      double *theta           = (double *)Calloc(length_theta, double);
      double *y               = (double *)Calloc(length_y, double);
      double *y_low           = (double *)Calloc(length_y, double);
      double *y_high          = (double *)Calloc(length_y, double);
      double *y_mean          = (double *)Calloc(length_y, double);
      double *y_var           = (double *)Calloc(length_y, double);
      double *u_mean          = (double *)Calloc(length_y, double);
      double *u_var           = (double *)Calloc(length_y, double);
      int *bin                = (int *)Calloc(length_y, int);
      int *bin_low            = (int *)Calloc(length_y, int);
      int *bin_high           = (int *)Calloc(length_y, int);
      double *u               = (double *)Calloc(length_y, double);
      double *u_low           = (double *)Calloc(length_y, double);
      double *u_high          = (double *)Calloc(length_y, double);
      double *can_u           = (double *)Calloc(length_y, double);
      double *can_u_low       = (double *)Calloc(length_y, double);
      double *can_u_high      = (double *)Calloc(length_y, double);
      double *log_like        = (double *)Calloc(length_y, double);
      double *can_log_like    = (double *)Calloc(length_y, double);
      double *ll_sum          = (double *)Calloc(length_ll_sum, double);
      double *can_ll_sum      = (double *)Calloc(length_ll_sum, double);
      
      
      /***** initialize variables for the prior ****/
      for(i = 0; i < length_thetastar; i++) {
      thetastar[i] = 0.0;
      }
      for(m = 1; m < n_bf; m++) {
      thetastar[m] = 1.0;
      }
      for(i = 0; i < length_mu; i++) {
      mu[i] = mu_0[i];
      }
      if(indicator_mu == 1) {
      for(m = 0; m < n_bf; m++) {
      if(m == 0) {
      i = 0;
      } else {
      i = 1;
      }
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = i + 2 * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      } else {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = m + n_bf * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      }
      
      for(i = 0; i < 5; i++) {
      log_dets[i] = 0.0;
      }
      
      
      double rho = 0.5 * max_rho;
      double sd_rho = 1.0;
      int acc_rho = 0;
      for(i = 0; i < n_clusters_2; i++) {
      cov_clusters[i] = 0.0;
      prec_clusters[i] = 0.0;
      }
      for(i = 0; i < n_clusters; i++) {
      k = i + n_clusters * i;
      cov_clusters[k] = 1.0;
      prec_clusters[k] = 1.0;
      }
      if(correlation_type == 1) {
      MakeAutoregressivePrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      if(correlation_type == 2) {
      MakeSpatialPrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      for(i = 0; i < n_responses_2; i++) {
      cov_responses[i] = 0.0;
      prec_responses[i] = 0.0;
      }
      for(i = 0; i < n_responses; i++) {
      k = i + n_responses * i;
      cov_responses[k] = 1.0;
      prec_responses[k] = 1.0;
      }
      InvertSymmetricMatrix(&n_bf, prec_bf, dummy_matrix, &log_dets[0]);
      InvertSymmetricMatrix(&n_preds, prec_preds, dummy_matrix, &log_dets[1]);
      InvertSymmetricMatrix(&dim_time, prec_time, dummy_matrix, &log_dets[2]);
      InvertSymmetricMatrix(&n_clusters, prec_clusters, cov_clusters, &log_dets[3]);
      InvertSymmetricMatrix(&n_responses, prec_responses, cov_responses,
      &log_dets[4]);
      
      for(i = 0; i < 5; i++) {
      log_dets[i] *= -1.0;
      }
      
      for(i = 0; i < length_thetastar_2; i++) {
      prec_thetastar[i] = 0.0;
      }
      MakePrecisionThetastar(dim_thetastar, prec_bf, prec_preds, prec_time,
      prec_clusters, prec_responses, prec_thetastar);
      
      double lp_thetastar = 0.0;
      LogPriorThetastar(dim_thetastar, log_dets,
      e_thetastar, prec_thetastar, &lp_thetastar);
      for(i = 0; i < length_xi; i++) {
      xi_low[i] = -0.4;
      xi_high[i] = -0.4;
      }
      
      /***** initialize variables for the tail ****/
      double tau2_xi_low = 1.0;
      double tau2_xi_high = 1.0;
      double lp_xi_low = 0.0;
      double lp_xi_high = 0.0;
      
      for(h = 0; h < n_respones; h++) {
      mu_xi_low[h] = mu_0_xi_low[h];
      mu_xi_high[h] = mu_0_xi_high[h];
      for(g = 0; g < n_clusters; g++) {
      k = g + n_clusters * h;
      e_xi_low[k] = xi_low[k] - mu_xi_low[h];
      e_xi_high[k] = xi_high[k] - mu_xi_high[h];
      }
      }
      /*indicator if pareto tail can switch between exponentail and pareto*/
      int update_pareto_low = 1;
      int update_pareto_high = 1;
      if (pareto_prob_low == 0.0 || pareto_prob_low == 1.0 || tau_low == 0.0) {
      update_pareto_low = 0;
      }
      if (pareto_prob_high == 0.0 || pareto_prob_high == 1.0 || tau_high == 1.0) {
      update_pareto_high = 0;
      }
      /*indicator for pareto tail (1) vs exponential tail (0)*/
      int pareto_low = 1;
      int pareto_high = 1;
      if(pareto_prob_low == 0.0) {
      pareto_low = 0;
      }
      if(pareto_prob_high == 0.0) {
      pareto_high = 0;
      }
      if(tau_low > 0.0 && pareto_low == 1) {
      LogPriorXi(dim_thetastar, log_dets,
      e_xi_low, prec_clusters, &tau2_xi_low, &lp_xi_low);
      }
      if(tau_high < 1.0 && pareto_high == 1) {
      LogPriorXi(dim_thetastar, log_dets,
      e_xi_high, prec_clusters, &tau2_xi_high, &lp_xi_high);
      }
      /***initialize variables for Copula***/
      
      
      
      
      /***initialize variables for UpdateThetastar***/
      for(i = 0; i < length_theta; i++) {
      theta[i] = 0.0;
      }
      for(g = 0; g < n_clusters; g++) {
      ThetastarToTheta(dim_thetastar, thetastar, theta, timepoints, &g);
      }
      
      
      for(i = 0; i < length_y; i++) {
      y[i] = y_array[0 + 7 * i];
      y_low[i] = y_array[1 + 7 * i];
      y_high[i] = y_array[2 + 7 * i];
      y_mean[i] = 0.0;
      y_var[i] = 0.0;
      u_mean[i] = 0.0;
      u_var[i] = 0.0;
      u[i] = 0.0;
      u_low[i] = 0.0;
      u_high[i] = 0.0;
      bin[i] = 0;
      bin_low[i] = 0;
      bin_high[i] = 0;
      log_like[i] = 0.0;
      }
      
      for(g = 0; g < length_ll_sum; g++) {
      ll_sum[g] = 0.0;
      }
      for(i = 0; i < 5; i++) {
      log_dets[i] = 0.0;
      }
      
      for(j = 0; j < n_timepoints; j++) {
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      tail_index = g + n_clusters * h;
      time_index = j + n_timepoints * tail_index;
      Likelihood(dim_thetastar, n_y, &basis_ind[h], basis_knots,
      y, y_low, y_high, status, x, u, u_low, u_high,
      bin, bin_low, bin_high, log_like, &ll_sum[time_index],
      theta, &shape[h], tau_low, tau_high,
      pareto_low, pareto_high,
      &xi_low[tail_index], &xi_high[tail_index], &g, &h);
      }
      }
      }
      
      
      
      }
      
      
      
      
      
      
      
      double *u, double *u_low, double *u_high,
      double *can_u, double *can_u_low, double *can_u_high,
      int *bin, int *bin_low, int *bin_high,
      double *log_like, double *ll_sum,
      double *can_log_like, double *can_ll_sum,
      double *lp_thetastar, double *log_dets, double *mu_array,
      double *e_thetastar, double *prec_thetastar,
      double *thetastar, double *theta,
      double *can_thetastar, double *can_theta,
      double *shape,
      double *tuning_thetastar, int *acc_thetastar, int *att_thetastar,
      double *tau_low, double *tau_high,
      int *pareto_low, int *pareto_high,
      double *xi_low, double *xi_high,
      int *cop, int *n_preds_q, double *timepoints,
      double *u2, double *w, double *w2, double *d2,
      double *log_like_copula, double *ll_sum_copula,
      double *can_u2, double *can_w, double *can_w2,
      double *can_log_like_copula, double *can_ll_sum_copula,
      double *gamma, double *eta)
      
      
      
      
      int iters = burn + keep;
      int num_keep =  keep / thin;
      int acc_rho = 0;
      
      
      
      
      
      
      
      
      
      Free(dim_thetastar);
      Free(n_y);
      Free(pred_mns);
      Free(pred_sds);
      Free(response_mns);
      Free(response_sds);
      
      
      /********* helpful below here **************/
      
      int max_dim = n_bf;
      if(max_dim < n_preds) {
      max_dim = n_preds;
      }
      if(max_dim < dim_time) {
      max_dim = dim_time;
      }
      int max_dim_2 = max_dim * max_dim;
      
      
      
      
      double rho = 0.5 * max_rho;
      double lp_thetastar = 0.0;
      double sd_rho = 1.0;
      double tau2_xi_low = 1.0;
      double tau2_xi_high = 1.0;
      double lp_xi_low = 0.0;
      double lp_xi_high = 0.0;
      
      double *pred_mns = (double *)Calloc(n_preds, double);
      double *pred_sds = (double *)Calloc(n_preds, double);
      double *response_mns = (double *)Calloc(n_responses, double);
      double *response_sds = (double *)Calloc(n_responses, double);
      double *mu = (double *)Calloc(length_mu, double);
      double *mu_array = (double *)Calloc(length_thetastar, double);
      double *e_thetastar = (double *)Calloc(length_thetastar, double);
      double *cov_clusters = (double *)Calloc(n_clusters_2, double);
      double *prec_clusters = (double *)Calloc(n_clusters_2, double);
      double *cov_responses = (double *)Calloc(n_clusters_2, double);
      double *prec_responses = (double *)Calloc(n_clusters_2, double);
      double *dummy_matrix = (double *)Calloc(max_dim_2, double);
      double *prec_thetastar = (double *)Calloc(length_thetastar_2, double);
      double *mu_xi_low = (double *)Calloc(n_responses, double);
      double *mu_xi_high = (double *)Calloc(n_responses, double);
      double *mu_0_xi_low = (double *)Calloc(n_responses, double);
      double *mu_0_xi_high = (double *)Calloc(n_responses, double);
      double *e_xi_low = (double *)Calloc(n_clusters * n_responses, double);
      double *e_xi_high = (double *)Calloc(n_clusters * n_responses, double);
      double *log_dets = (double *)Calloc(5, double);
      
      for(i = 0; i < n_preds; i++) {
      pred_mns[i] = mcmc_doubles[i];
      pred_sds[i] = mcmc_doubles[i + n_preds];
      }
      k = 2 * n_preds;
      l = 2 * (n_preds + n_responses) + 11;
      for(i = 0; i < n_responses; i++) {
      response_mns[i] = mcmc_doubles[k + i];
      response_sds[i] = mcmc_doubles[k + n_responses + i];
      mu_0_xi_low[i] = mcmc_doubles[l + i];
      mu_0_xi_high[i] = mcmc_doubles[l + n_responses + i];
      }
      for(i = 0; i < length_mu; i++) {
      mu[i] = mu_0[i];
      }
      if(indicator_mu == 1) {
      for(m = 0; m < n_bf; m++) {
      if(m == 0) {
      i = 0;
      } else {
      i = 1;
      }
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = i + 2 * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      } else {
      for(m = 0; m < n_bf; m++) {
      for(p = 0; p < n_preds; p++) {
      for(d = 0; d < dim_time; d++) {
      k = m + n_bf * (p + n_preds * d);
      for(g = 0; g < n_clusters; g++) {
      for(h = 0; h < n_responses; h++) {
      l = m + n_bf * p + n_bf * n_preds * d +
      n_bf * n_preds * dim_time * (g + n_clusters * h);
      mu_array[l] = mu[k];
      e_thetastar[l] = thetastar[l] - mu_array[l];
      }
      }
      }
      }
      }
      }
      for(i = 0; i < n_responses_2; i++) {
      cov_responses[i] = 0.0;
      prec_responses[i] = 0.0;
      }
      for(i = 0; i < n_responses; i++) {
      k = i + n_responses * i;
      cov_responses[k] = 1.0;
      prec_responses[k] = 1.0;
      }
      
      for(i = 0; i < 5; i++) {
      log_dets[i] = 0.0;
      }
      for(i = 0; i < max_dim_2; i++) {
      dummy_matrix[i] = 0.0;
      }
      for(i = 0; i < n_clusters_2; i++) {
      cov_clusters[i] = 0.0;
      prec_clusters[i] = 0.0;
      }
      for(i = 0; i < n_clusters; i++) {
      k = i + n_clusters * i;
      cov_clusters[k] = 1.0;
      prec_clusters[k] = 1.0;
      }
      if(correlation_type == 1) {
      MakeAutoregressivePrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      if(correlation_type == 2) {
      MakeSpatialPrecision(&n_clusters, &rho, distance_matrix,
      prec_clusters, &log_dets[3]);
      }
      InvertSymmetricMatrix(&n_bf, prec_bf, dummy_matrix, &log_dets[0]);
      InvertSymmetricMatrix(&n_preds, prec_preds, dummy_matrix, &log_dets[1]);
      InvertSymmetricMatrix(&dim_time, prec_time, dummy_matrix, &log_dets[2]);
      InvertSymmetricMatrix(&n_clusters, prec_clusters, cov_clusters, &log_dets[3]);
      InvertSymmetricMatrix(&n_responses, prec_responses, cov_responses,
      &log_dets[4]);
      
      for(i = 0; i < 5; i++) {
      log_dets[i] *= -1.0;
      }
      for(i = 0; i < length_thetastar_2; i++) {
      prec_thetastar[i] = 0.0;
      }
      MakePrecisionThetastar(dim_thetastar,
      prec_bf, prec_preds, prec_time, prec_clusters, prec_responses,
      prec_thetastar);
      
      LogPriorThetastar(dim_thetastar, log_dets, e_thetastar, prec_thetastar,
      &lp_thetastar);
      
      for(h = 0; h < n_responses; h++) {
      mu_xi_low[h] = mu_0_xi_low[h];
      mu_xi_high[h] = mu_0_xi_high[h];
      }
      for(i = 0; i < n_clusters; i++) {
      for(j = 0; j < n_responses; j++) {
      k = i * n_responses + j;
      e_xi_low[k] = xi_low[k] - mu_xi_low[j];
      e_xi_high[k] = xi_high[k] - mu_xi_high[j];
      }
      }
      
      if(tau_low > 0.0) {
      LogPriorXi(dim_thetastar, log_dets, e_xi_low, prec_clusters, &tau2_xi_low,
      &lp_xi_low);
      }
      if(tau_high < 1.0) {
      LogPriorXi(dim_thetastar, log_dets, e_xi_high, prec_clusters, &tau2_xi_high,
      &lp_xi_high);
      }
      
      
      
      
      
      }
      
      
      
      void make_b (
      void q_ptr(double *, double *, double *, double *, double*),
      int *n_bf, double *basis_knots,
      double *shape, double *basis_matrix) {
      
      int i, j, k;
      double mn = 0; /*for identifiability mean is 0 and scale is 1*/
      double scale = 1;
      double previous, current;
      int kappa_length = *n_bf + 4;
      int dummy_bin = 3;
      if (*shape < 0) {
      for (i = 0; i < kappa_length; i++) {
      basis_matrix[i] = 1.0;
      }
      for (i = 0; i < kappa_length; i++) {
      find_interval(&kappa_length, basis_knots, &basis_knots[i], &dummy_bin);
      for (j = 1; j < *n_bf; j++) {
      k = j - 1;
      CubicISpline(&basis_knots[i], &kappa_length, basis_knots, &k, &dummy_bin, &basis_matrix[i + kappa_length * j]);
      }
      }
      } else {
      /*first column*/
      for (i = 1; i < *n_bf; i++) {
      q_ptr(&basis_knots[1], &mn, &scale, shape, &basis_matrix[i]);
      }
      /*second column to (n_bf - 1)th column*/
      for (j = 1; j < (*n_bf - 2); j++) {
      q_ptr(&basis_knots[j], &mn, &scale, shape, &previous);
      q_ptr(&basis_knots[j + 1], &mn, &scale, shape, &current);
      for (i = (j + 1); i < *n_bf; i++) {
      basis_matrix[i + *n_bf * j] = current - previous;
      }
      }
      }
      }
      
      
      ///**** make sure that for interval censored y_high < y_low ***/
      void clike_spline(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      void q_ptr(double *, double *, double *, double *, double*),
      void p_ptr(double *, double *, double *, double *, double*),
      void log_d_ptr(double *, double *, double *, double *, double*),
      double *AAAA,
      double *B,
      int *bin, int *bin_low, int *bin_high,
      double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      double *w, double *M_low, double *M_high, double *I_low, double *I_high,
      int *g, int *j, int *h
      ) {
      double z, a, b, c, d, f, f_prime;
      double dummy_low = 0.0;
      double dummy_high = 0.0;
      
      double eps = pow(10, -5);
      int i, k, m, p, this_i, iter;
      int this_model = *g + *G * *j + *G * *J * *h;
      int this_obs = *G * *J * *H;
      int this_x_a = *g + *G * *j;
      int this_x_b = *G * *J;
      
      int M_1 = *M - 1;
      int kappa_length = *M + 4;
      int dummy_bin = 0;
      
      double *BM = (double *)Calloc((kappa_length * *P), double);
      double *breaks = (double *)Calloc(kappa_length, double);
      
      double *I_low_p = (double *)Calloc(*P, double);
      double *m_low_p = (double *)Calloc(*P, double);
      double *I_high_p = (double *)Calloc(*P, double);
      double *m_high_p = (double *)Calloc(*P, double);
      
      double sig_low, sig_high, low_thresh, high_thresh, U1, U2;
      
      M_low[0] = 0.0;
      M_high[0] = 0.0;
      I_low[0] = 1.0;
      I_high[0] = 1.0;
      
      for (m = 0; m < M_1; m++) {
      CubicMSpline(tau_low, &kappa_length, kappa,   &m,  &dummy_bin, &M_low[m + 1]);
      M_low[m + 1] *= *tau_low;
      CubicMSpline(tau_high, &kappa_length, kappa,   &m, &dummy_bin, &M_high[m + 1]);
      M_high[m + 1] *= (1 - *tau_high);
      CubicISpline(tau_low,  &kappa_length, kappa, &m, &dummy_bin, &I_low[m + 1]);
      CubicISpline(tau_high,  &kappa_length, kappa, &m, &dummy_bin, &I_high[m + 1]);
      }
      
      for (p = 0; p < *P; p++) {
      I_low_p[p] = 0.0;
      m_low_p[p] = 0.0;
      I_high_p[p] = 0.0;
      m_high_p[p] = 0.0;
      for (m = 0; m < *M; m++) {
      k = this_model + this_obs * (m + *M * p);
      m_low_p[p] += theta[k] * M_low[m];
      I_low_p[p] += theta[k] * I_low[m];
      m_high_p[p]   += theta[k] * M_high[m];
      I_high_p[p]     += theta[k] * I_high[m];
      }
      }
      for (i = 0; i < kappa_length; i++) {
      for (p = 0; p < *P; p++) {
      BM[i + kappa_length * p] = 0.0;
      for (m = 0; m < *M; m++) {
      BM[i + kappa_length * p] += B[i + kappa_length * m] * theta[this_model + this_obs * (m + *M * p)];
      }
      dummy_low += BM[i + kappa_length * p];
      }
      }
      *ll_sum = 0;
      
      //*** status 1 is continuous *****//
      //*** status 2 is left censored (censored below) *****//
      //*** status 3 is right censored (censored above) *****//
      //*** status 4 is interval censored *****//
      for (i = 0; i < *n_gjh; i++) {
      this_i = this_model + this_obs * i;
      if (status[this_i] == -99) {
      log_like[this_i] = 0.0;
      goto LL_SUM;
      }
      low_thresh = I_low_p[0] * X[this_x_a + this_x_b * i];
      high_thresh = I_high_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      low_thresh += X[this_x_a + this_x_b * (i + *N * p)] * I_low_p[p]; high_thresh += X[this_x_a + this_x_b * (i + *N * p)] * I_high_p[p];
      }
      if (status[this_i] == 0 || status[this_i] == 1) {
      if (low_thresh > y[this_i]) {
      sig_low = m_low_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_low += X[this_x_a + this_x_b * (i + *N * p)] * m_low_p[p];
      }
      z = (low_thresh - y[this_i]) / sig_low;
      if (*Pareto_low == 0) {
      U[this_i]         =  *tau_low * exp(-z);
      log_like[this_i]    = log(*tau_low) - log(sig_low) - z;
      } else{
      U[this_i]         = *tau_low * (pow(1 +  exp(*xi_low) * z, -1 / exp(*xi_low)));
      log_like[this_i]    = log(*tau_low) - log(sig_low) - (1 / exp(*xi_low) + 1) * log(1 + exp(*xi_low) * z);
      }
      if (U[this_i] <= 0) {
      U[this_i] = 0.0000001;
      }
      goto LL_SUM;
      }
      if (high_thresh < y[this_i]) {
      sig_high = m_high_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_high += X[this_x_a + this_x_b * (i + *N * p)] * m_high_p[p];
      }
      z = (y[this_i] - high_thresh) / sig_high;
      if (*Pareto_high == 0) {
      U[this_i]         =   *tau_high + (1 - *tau_high) * (1 - exp(-z));
      log_like[this_i]    =   log(1 - *tau_high) - log(sig_high) - z;
      }
      else {
      U[this_i]         = *tau_high + (1 - *tau_high) * (1 - pow(1 +  exp(*xi_high) * z,-1 / exp(*xi_high)));
      log_like[this_i]    = log(1 - *tau_high) - log(sig_high) - (1 / exp(*xi_high) + 1) * log(1 + exp(*xi_high) * z);
      }
      if (U[this_i] >= 1.0) {
      U[this_i] = 0.9999999;
      }
      goto LL_SUM;
      }
      for (m = 0; m < *M; m++) {
      w[m] = 0.0;
      dummy_low = 0.0;
      dummy_high = 0.0;
      for (p = 0; p < *P; p++) {
      z = X[this_x_a + this_x_b * (i + *N * p)];
      w[m]        += theta[this_model + this_obs * (m + *M * p)] * z;
      dummy_low   += BM[ bin[this_i]      + kappa_length * p] * z;
      dummy_high  += BM[(bin[this_i] + 1) + kappa_length * p] * z;
      }
      }
      if (dummy_low > y[this_i] || dummy_high < y[this_i]) {
      for (k = 0; k < kappa_length; k++) {
      breaks[k] = 0.0;
      for (p = 0; p < *P; p++) {
      breaks[k] += BM[k + kappa_length * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(&kappa_length, breaks, &y[this_i], &bin[this_i]);
      }
      a = 0.0;
      b = 0.0;
      c = 0.0;
      d = w[0];
      for (m = 1; m < *M; m++) {
      k = 4 * (m + *M * bin[this_i]);
      a += AAAA[0 + k] * w[m];
      b += AAAA[1 + k] * w[m];
      c += AAAA[2 + k] * w[m];
      d += AAAA[3 + k] * w[m];
      }
      if (kappa[bin[this_i]] > U[this_i] || kappa[bin[this_i] + 1] < U[this_i]) {
      U[this_i] =  (kappa[bin[this_i]] +   kappa[bin[this_i] + 1]) / 2;
      }
      f           =     a * pow(U[this_i], 3) +     b * pow(U[this_i], 2) + c * U[this_i] + d;
      f_prime     = 3 * a * pow(U[this_i], 2) + 2 * b *     U[this_i]     + c;
      iter  = 0;
      while(fabs(y[this_i] - f) > eps && iter < 1000) {
      iter ++;
      U[this_i] = U[this_i] + (y[this_i] - f) / f_prime;
      if (U[this_i] < kappa[bin[this_i]] || U[this_i] > kappa[bin[this_i] + 1]) {
      U1 = kappa[bin[this_i]];
      U2 = kappa[bin[this_i] + 1];
      U[this_i] = (U1 + U2) / 2;
      while(fabs(y[this_i] - f) > eps && iter < 1000) {
      iter ++;
      f = a * pow(U[this_i], 3) + b * pow(U[this_i], 2) + c * U[this_i] + d;
      if (f < y[this_i]) {
      U1 = U[this_i];
      }else{
      U2 = U[this_i];
      }
      U[this_i] = (U1 + U2) / 2;
      }
      }
      f           =     a * pow(U[this_i], 3) +     b * pow(U[this_i], 2) + c * U[this_i] + d;
      f_prime     = 3 * a * pow(U[this_i], 2) + 2 * b *     U[this_i]     + c;
      }
      if (iter == 1000) {U[this_i] = (kappa[bin[this_i]] + kappa[bin[this_i] + 1]) / 2; Rprintf("max_iter!");}
      log_like[this_i] = -log(3 * a * pow(U[this_i], 2) + 2 * b * U[this_i]   + c);
      goto LL_SUM;
      }
      if (status[this_i] == 2 || status[this_i] == 4) {
      if (low_thresh > y_high[this_i]) {
      sig_low = m_low_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_low += X[this_x_a + this_x_b * (i + *N * p)] * m_low_p[p];
      }
      z = (low_thresh - y_high[this_i]) / sig_low;
      if (*Pareto_low == 0) {
      U_high[this_i]         =  *tau_low * exp(-z);
      } else {
      U_high[this_i]         = *tau_low * (pow(1 +  exp(*xi_low) * z,-1 / exp(*xi_low)));
      }
      if (U_high[this_i] <= 0) {
      U_high[this_i] = 0.0000001;
      } /*yucky fix*/
      goto LAST_STATUS;
      }
      if (high_thresh < y_high[this_i]) {
      sig_high = m_high_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_high += X[this_x_a + this_x_b * (i + *N * p)] * m_high_p[p];
      }
      z = (y_high[this_i] - high_thresh) / sig_high;
      if (*Pareto_high == 0) {
      U_high[this_i]         =   *tau_high + (1 - *tau_high) * (1 - exp(-z));
      }
      else{
      U_high[this_i]         = *tau_high + (1 - *tau_high) * (1 - pow(1 +  exp(*xi_high) * z,-1 / exp(*xi_high)));
      }
      if (U_high[this_i] >= 1.0) {
      U_high[this_i] = 0.9999999;
      }
      goto LAST_STATUS;
      }
      for (m = 0; m < *M; m++) {
      w[m] = 0.0;
      dummy_low = 0.0;
      dummy_high = 0.0;
      for (p = 0; p < *P; p++) {
      z = X[this_x_a + this_x_b * (i + *N * p)];
      w[m]        += theta[this_model + this_obs * (m + *M * p)] * z;
      dummy_low   += BM[ bin_high[this_i]      + kappa_length * p] * z;
      dummy_high  += BM[(bin_high[this_i] + 1) + kappa_length * p] * z;
      }
      }
      if (dummy_low > y_high[this_i] || dummy_high < y_high[this_i]) {
      for (k = 0; k < kappa_length; k++) {
      breaks[k] = 0.0;
      for (p = 0; p < *P; p++) {
      breaks[k] += BM[k + kappa_length * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(&kappa_length, breaks, &y_high[this_i], &bin_high[this_i]);
      }
      a = 0.0;
      b = 0.0;
      c = 0.0;
      d = w[0];
      for (m = 1; m < *M; m++) {
      k = 4 * (m + *M * bin_high[this_i]);
      a += AAAA[0 + k] * w[m];
      b += AAAA[1 + k] * w[m];
      c += AAAA[2 + k] * w[m];
      d += AAAA[3 + k] * w[m];
      }
      if (kappa[bin_high[this_i]] > U_high[this_i] || kappa[bin_high[this_i] + 1] < U_high[this_i]) {
      U_high[this_i] =  (kappa[bin_high[this_i]] +   kappa[bin_high[this_i] + 1]) / 2;
      }
      f           =     a * pow(U_high[this_i], 3) +     b * pow(U_high[this_i], 2) + c * U_high[this_i] + d;
      f_prime     = 3 * a * pow(U_high[this_i], 2) + 2 * b *     U_high[this_i]     + c;
      iter  = 0;
      while(fabs(y_high[this_i] - f) > eps && iter < 1000) {
      iter ++;
      U_high[this_i] = U_high[this_i] + (y_high[this_i] - f) / f_prime;
      if (U_high[this_i] < kappa[bin_high[this_i]] || U_high[this_i] > kappa[bin_high[this_i] + 1]) {
      U1 = kappa[bin_high[this_i]];
      U2 = kappa[bin_high[this_i] + 1];
      U_high[this_i] = (U1 + U2) / 2;
      while(fabs(y_high[this_i] - f) > eps && iter < 1000) {
      iter ++;
      f = a * pow(U_high[this_i], 3) + b * pow(U_high[this_i], 2) + c * U_high[this_i] + d;
      if (f < y_high[this_i]) {U1 = U_high[this_i];}else{U2 = U_high[this_i];}
      U_high[this_i] = (U1 + U2) / 2;
      }
      }
      f           =     a * pow(U_high[this_i], 3) +     b * pow(U_high[this_i], 2) + c * U_high[this_i] + d;
      f_prime     = 3 * a * pow(U_high[this_i], 2) + 2 * b *     U_high[this_i]     + c;
      }
      if (iter == 1000) {U_high[this_i] = (kappa[bin_high[this_i]] + kappa[bin_high[this_i] + 1]) / 2; Rprintf("max_iter!");}
      }
      LAST_STATUS:
      if (status[this_i] == 3 || status[this_i] == 4) {
      if (low_thresh > y_low[this_i]) {
      sig_low = m_low_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_low += X[this_x_a + this_x_b * (i + *N * p)] * m_low_p[p];
      }
      z = (low_thresh - y_low[this_i]) / sig_low;
      if (*Pareto_low == 0) {
      U_low[this_i]         =  *tau_low * exp(-z);
      } else{
      U_low[this_i]         = *tau_low * (pow(1 +  exp(*xi_low) * z,-1 / exp(*xi_low)));
      }
      if (U_low[this_i] <= 0) {
      U_low[this_i] = 0.0000001;
      }
      goto TAU_DIFF;
      }
      if (high_thresh < y_low[this_i]) {
      sig_high = m_high_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      sig_high += X[this_x_a + this_x_b * (i + *N * p)] * m_high_p[p];
      }
      z = (y_low[this_i] - high_thresh) / sig_high;
      if (*Pareto_high == 0) {
      U_low[this_i]         =   *tau_high + (1 - *tau_high) * (1 - exp(-z));
      } else {
      U_low[this_i]         = *tau_high + (1 - *tau_high) * (1 - pow(1 +  exp(*xi_high) * z,-1 / exp(*xi_high)));
      }
      if (U_low[this_i] >= 1.0) {
      U_low[this_i] = 0.9999999;
      }
      goto TAU_DIFF;
      }
      for (m = 0; m < *M; m++) {
      w[m] = 0.0;
      dummy_low = 0.0;
      dummy_high = 0.0;
      for (p = 0; p < *P; p++) {
      z = X[this_x_a + this_x_b * (i + *N * p)];
      w[m]        += theta[this_model + this_obs * (m + *M * p)] * z;
      dummy_low   += BM[ bin_low[this_i]      + kappa_length * p] * z;
      dummy_high  += BM[(bin_low[this_i] + 1) + kappa_length * p] * z;
      }
      }
      if (dummy_low > y_low[this_i] || dummy_high < y_low[this_i]) {
      for (k = 0; k < kappa_length; k++) {
      breaks[k] = 0.0;
      for (p = 0; p < *P; p++) {
      breaks[k] += BM[k + kappa_length * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(&kappa_length, breaks, &y_low[this_i], &bin_low[this_i]);
      }
      a = 0.0;
      b = 0.0;
      c = 0.0;
      d = w[0];
      for (m = 1; m < *M; m++) {
      k = 4 * (m + *M * bin_low[this_i]);
      a += AAAA[0 + k] * w[m];
      b += AAAA[1 + k] * w[m];
      c += AAAA[2 + k] * w[m];
      d += AAAA[3 + k] * w[m];
      }
      if (kappa[bin_low[this_i]] > U_low[this_i] || kappa[bin_low[this_i] + 1] < U_low[this_i]) {
      U_low[this_i] =  (kappa[bin_low[this_i]] +   kappa[bin_low[this_i] + 1]) / 2;
      }
      f           =     a * pow(U_low[this_i], 3) +     b * pow(U_low[this_i], 2) + c * U_low[this_i] + d;
      f_prime     = 3 * a * pow(U_low[this_i], 2) + 2 * b *     U_low[this_i]     + c;
      iter  = 0;
      while(fabs(y_low[this_i] - f) > eps && iter < 1000) {
      iter ++;
      U_low[this_i] = U_low[this_i] + (y_low[this_i] - f) / f_prime;
      if (U_low[this_i] < kappa[bin_low[this_i]] || U_low[this_i] > kappa[bin_low[this_i] + 1]) {
      U1 = kappa[bin_low[this_i]];
      U2 = kappa[bin_low[this_i] + 1];
      U_low[this_i] = (U1 + U2) / 2;
      while(fabs(y_low[this_i] - f) > eps && iter < 1000) {
      iter ++;
      f = a * pow(U_low[this_i], 3) + b * pow(U_low[this_i], 2) + c * U_low[this_i] + d;
      if (f < y_low[this_i]) {
      U1 = U_low[this_i];
      }else{
      U2 = U_low[this_i];
      }
      U_low[this_i] = (U1 + U2) / 2;
      }
      }
      f           =     a * pow(U_low[this_i], 3) +     b * pow(U_low[this_i], 2) + c * U_low[this_i] + d;
      f_prime     = 3 * a * pow(U_low[this_i], 2) + 2 * b *     U_low[this_i]     + c;
      }
      if (iter == 1000) {U_low[this_i] = (kappa[bin_low[this_i]] + kappa[bin_low[this_i] + 1]) / 2; Rprintf("max_iter!");}
      }
      TAU_DIFF:
      log_like[this_i] =  log(U_high[this_i] - U_low[this_i]);
      U[this_i]      = (U_high[this_i] + U_low[this_i]) / 2;
      LL_SUM:
      *ll_sum += log_like[this_i];
      }
      Free(BM);
      Free(breaks);
      }
      
      void clike_spline_wrapper(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      double *w, double *M_low, double *M_high, double *I_low, double *I_high,
      int *g, int *j, int *h
      ) {
      
      int M_1 = *M - 1;
      int m;
      
      double *AAAA = (double *)Calloc((4 * *M  * M_1), double);
      double *B    = (double *)Calloc(((*M + 4)  * *M), double);
      
      void (*q_ptr)(double *, double *, double *, double *, double*);
      void (*p_ptr)(double *, double *, double *, double *, double*);
      void (*log_d_ptr)(double *, double *, double *, double *, double*);
      
      for (m = 0; m < 4 * *M * M_1; m++) {
      AAAA[m] = 0.0;
      }
      
      MakeCubicSplineCoefficients(M, kappa, AAAA);
      
      q_ptr               = &QNorm;
      p_ptr               = &PNorm;
      log_d_ptr           = &LogDNorm;
      //
      int u = *g + *G * *j + *G * *J * *h;
      int tail_index = *g + *G * *h;
      
      shape[*h] = -1.0;
      for (m = 0; m < (*M + 4)  * *M; m++) {
      B[m] = 0.0;
      }
      make_b(q_ptr, M, kappa, &shape[*h], B);
      
      clike_spline(G, J, H, M, P,
      N, &n_gjh[u], status, kappa,
      X, y, y_low, y_high, theta, &shape[*h],
      q_ptr, p_ptr, log_d_ptr, AAAA,
      B,
      bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
      Pareto_low, Pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
      w, M_low, M_high, I_low, I_high,
      g, j, h);
      
      Free(AAAA);
      Free(B);
      }
      
      /**** make sure that for interval censored y_high < y_low ***/
      void clike (int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      void q_ptr(double *, double *, double *, double *, double*),
      void p_ptr(double *, double *, double *, double *, double*),
      void log_d_ptr(double *, double *, double *, double *, double*),
      double *AAAA,
      double *B,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      double *w, double *M_low, double *M_high, double *I_low, double *I_high,
      int *g, int *j, int *h
      ) {
      
      int i, k, m, p, this_i;
      int this_model = *g + *G * *j + *G * *J * *h;
      int this_obs = *G * *J * *H;
      int this_x_a = *g + *G * *j;
      int this_x_b = *G * *J;
      
      int M_1 = *M - 1;
      
      double mn, s;
      double dummy_mean = 0;
      double dummy_scale = 1;
      
      double *BM = (double *)Calloc((*M * *P), double);
      double *breaks = (double *)Calloc((*M), double);/*breaks for each iteration*/
      double *q_kappa = (double *)Calloc((M_1), double);
      double *mean_low_p = (double *)Calloc(*P, double);
      double *sd_low_p = (double *)Calloc(*P, double);
      
      double *mean_high_p = (double *)Calloc(*P, double);
      double *sd_high_p = (double *)Calloc(*P, double);
      double sig_low, sig_high, low_thresh, high_thresh, mean_low, sd_low, mean_high, sd_high;
      double *dummy_vec = (double *)Calloc((M_1), double);
      double *dummy_mat = (double *)Calloc((M_1 * *P), double);
      
      double *mn_1 = (double *)Calloc(*P, double);
      double *mn_2 = (double *)Calloc(*P, double);
      double *mn_3 = (double *)Calloc(*P, double);
      
      for (m = 0; m < M_1; m++) {
      for (p = 0; p < *P; p++) {
      dummy_mat[m + M_1 * p] = theta[this_model + this_obs * ((m + 1) + *M * p)];
      }
      }
      for (m = 1; m < M_1; m++) {
      q_ptr(&kappa[m], &dummy_mean, &dummy_scale, shape, &q_kappa[m]);
      }
      if (*tau_low >  0) {
      I_low[0] = 1.0;
      m = 1;
      q_ptr(tau_low, &dummy_mean, &dummy_scale, &shape[*h], &I_low[m]);
      if (*M > 2) {
      for (m = 2; m < *M; m++) {
      I_low[m] = 0.0;
      }
      }
      for (p = 0; p < *P; p++) {
      mean_low_p[p] = theta[this_model + this_obs * (0 + *M * p)];
      sd_low_p[p]   = theta[this_model + this_obs * (1 + *M * p)];
      }
      }
      if (*tau_high < 1) {
      dummy_vec[0] = I_high[1];
      I_high[0] = 1.0;
      m = 1;
      if (*M == 2) {
      q_ptr(tau_high, &dummy_mean, &dummy_scale, &shape[*h], &I_high[m]);
      }
      if (*M > 2) {
      I_high[m] = q_kappa[m];
      for (m = 2; m < (*M - 1); m++) {
      I_high[m] = q_kappa[m] - q_kappa[m - 1];
      }
      q_ptr(tau_high, &dummy_mean, &dummy_scale, &shape[*h], &I_high[M_1]);
      I_high[M_1] -= q_kappa[(*M - 2)];
      }
      for (p = 0; p < *P; p++) {
      mn_1[p] = theta[this_model + this_obs * (0 + *M * p)];
      }
      if (*M == 2) {
      for (p = 0; p < *P; p++) {
      mean_high_p[p] = mn_1[p];
      }
      }
      if (*M > 2) {
      for (m = 0; m < (*M - 2); m++) {
      dummy_vec[m] = I_high[m + 1];
      }
      int T_A = 1;
      int T_B = 0;
      int dummy_one = 1;
      //      matrix_multiply(&T_A, &T_B, P, &dummy_one, &M_1, dummy_mat, dummy_vec, mn_2);
      matrix_multiply(&T_A, &T_B, &M_1, P, P, &dummy_one, dummy_mat, dummy_vec, mn_2);
      for (p = 0; p < *P; p++) {
      mean_high_p[p] = mn_1[p] + mn_2[p] - q_kappa[M_1 - 1] * theta[this_model + this_obs * ((*M - 1) + *M * p)];
      }
      }
      for (p = 0; p < *P; p++) {
      sd_high_p[p] = theta[this_model + this_obs * ((*M - 1) + *M * p)];
      }
      }
      /**********construct the break matrix************/
      for (i = 1; i < (*M - 1); i++) {
      for (p = 0; p < *P; p++) {
      BM[i + *M * p] = 0;
      for (k = 0; k < i; k++) {
      BM[i + *M * p] += B[i + *M * k] * theta[this_model + this_obs * ((k + 1) + *M * p)];
      }
      }
      }
      *ll_sum = 0;
      breaks[0] = -1000000.0;
      breaks[*M - 1] = 1000000.0;
      
      low_thresh  =  -100000000.0;
      high_thresh =  100000000.0;
      
      //*** status 1 is continuous *****//
      //*** status 2 is left censored (censored below) *****//
      //*** status 3 is right censored (censored above) *****//
      //*** status 4 is interval censored *****//
      for (i = 0; i < *n_gjh; i++) {
      this_i = this_model + this_obs * i;
      if (status[this_i] == -99) {
      log_like[this_i] = 0.0; goto LL_SUM;
      }
      if (*tau_low > 0) {
      mean_low = X[this_x_a + this_x_b * (i + *N * 0)] * mean_low_p[0];
      sd_low = X[this_x_a + this_x_b * (i + *N * 0)] * sd_low_p[0];
      for (p = 1; p < *P; p++) {
      mean_low += X[this_x_a + this_x_b * (i + *N * p)] * mean_low_p[p];
      sd_low   += X[this_x_a + this_x_b * (i + *N * p)] * sd_low_p[p];
      }
      q_ptr(tau_low, &mean_low, &sd_low, shape, &low_thresh);
      }
      if (*tau_high < 1) {
      mean_high = X[this_x_a + this_x_b * (i + *N * 0)] * mean_high_p[0];
      sd_high = X[this_x_a + this_x_b * (i + *N * 0)] * sd_high_p[0];
      for (p = 1; p < *P; p++) {
      mean_high += X[this_x_a + this_x_b * (i + *N * p)] * mean_high_p[p];
      sd_high   += X[this_x_a + this_x_b * (i + *N * p)] * sd_high_p[p];
      }
      q_ptr(tau_high, &mean_high, &sd_high, shape, &high_thresh);
      }
      if (status[this_i] == 0 || status[this_i] == 1) { /*uncensored observations*/
      if (low_thresh > y[this_i]) {
      log_d_ptr(&low_thresh, &mean_low, &sd_low, shape, &sig_low);
      sig_low = *tau_low * exp(-sig_low);
      if (*Pareto_low == 0) {
      log_like[this_i] = log(*tau_low) - log(sig_low) - ((low_thresh - y[this_i]) / sig_low);
      U[this_i] =  *tau_low - *tau_low * (1 - exp(-(low_thresh - y[this_i]) / sig_low));
      } else{
      log_like[this_i] = (log(*tau_low) - log(sig_low) - (1 / exp(*xi_low) + 1) * log(1 + exp(*xi_low) * (low_thresh - y[this_i])/ sig_low));
      U[this_i] = (*tau_low - *tau_low * ( 1-  pow(1 + exp(*xi_low) * (low_thresh - y[this_i]) / sig_low, -1 / exp(*xi_low))));
      }
      if (U[this_i] <= 0) {
      U[this_i] = 0.0000001;
      }
      goto LL_SUM;
      }
      if (y[this_i] > high_thresh) {
      log_d_ptr(&high_thresh, &mean_high, &sd_high, shape, &sig_high);
      sig_high = (1 - *tau_high ) * exp(-sig_high);
      if (*Pareto_high == 0) {
      log_like[this_i] = log(1 - *tau_high) - log(sig_high) - ((y[this_i] - high_thresh) / sig_high);
      U[this_i] =    *tau_high + (1 - *tau_high) * (1 - exp(-(y[this_i] - high_thresh) / sig_high));
      } else{
      log_like[this_i] = log(1 - *tau_high) - log(sig_high) -(1 / exp(*xi_high) + 1) * log(1 + exp(*xi_high) * (y[this_i] - high_thresh) / sig_high);
      U[this_i] = *tau_high + (1 - *tau_high) * (1 - pow(1 + exp(*xi_high) * (y[this_i] - high_thresh) / sig_high, -1 / exp(*xi_high)));
      }
      if (U[this_i] >= 1) {
      U[this_i] = 0.9999999;
      }
      goto LL_SUM;
      }
      for (m = 1; m < (*M - 1); m++) { /*element of breaks*/
      breaks[m] = 0;
      for (p = 0; p < *P; p++) { /*sum across all predictors*/
      breaks[m] +=   (BM[m + *M * p] + theta[this_model + 0 + this_obs * *M * p]) * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(M, breaks, &y[this_i], &bin[this_i]);
      s = 0;
      if (bin[this_i] == 0) {
      mn = 0;
      for (p = 0; p < *P; p++) {
      mn += theta[this_model + 0 + this_obs * *M * p]  * X[this_x_a + this_x_b * (i + *N * p)];
      s += theta[this_model + this_obs * (1 + *M * p)] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      } else{
      for (p = 0; p < *P; p++) {
      s += theta[this_model + this_obs * (bin[this_i] + 1) + this_obs * *M * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      mn = breaks[bin[this_i]] - s * q_kappa[bin[this_i]];
      }
      log_d_ptr(&y[this_i], &mn, &s, shape, &log_like[this_i]);
      p_ptr(&y[this_i], &mn, &s, shape, &U[this_i]);
      goto LL_SUM;
      }
      if (status[this_i] == 2 || status[this_i] == 4) {
      if (low_thresh > y_high[this_i]) {
      log_d_ptr(&low_thresh, &mean_low, &sd_low, shape, &sig_low);
      sig_low = *tau_low * exp(-sig_low);
      if (*Pareto_low == 0) {
      U_high[this_i] =  *tau_low - *tau_low * (1 - exp(-(low_thresh - y_high[this_i]) / sig_low));
      } else{
      U_high[this_i] = (*tau_low - *tau_low * ( 1-  pow(1 + exp(*xi_low) * (low_thresh - y_high[this_i]) / sig_low, -1 / exp(*xi_low))));
      }
      if (U_high[this_i] <= 0) {
      U_high[this_i] = 0.0000001;
      }
      goto LAST_STATUS;
      }
      if (y_high[this_i] > high_thresh) {
      log_d_ptr(&high_thresh, &mean_high, &sd_high, shape, &sig_high);
      sig_high = (1 - *tau_high ) * exp(-sig_high);
      if (*Pareto_high == 0) {
      U_high[this_i] =    *tau_high + (1 - *tau_high) * (1 - exp(-(y_high[this_i] - high_thresh) / sig_high));
      } else{
      U_high[this_i] = *tau_high + (1 - *tau_high) * ( 1- pow(1 + exp(*xi_high) * (y_high[this_i] - high_thresh) / sig_high, -1 / exp(*xi_high)));
      }
      if (U_high[this_i] >= 1) {
      U_high[this_i] = 0.9999999;
      }
      goto LAST_STATUS;
      }
      for (m = 1; m < (*M - 1); m++) { /*element of breaks*/
      breaks[m] = 0;
      for (p = 0; p < *P; p++) { /*sum across all predictors*/
      breaks[m] +=   (BM[m + *M * p] + theta[this_model + 0 + this_obs * *M * p]) * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(M, breaks, &y_high[this_i], &bin_high[this_i]);
      s = 0;
      if (bin_high[this_i] == 0) {
      mn = 0;
      for (p = 0; p < *P; p++) {
      mn += theta[this_model + 0 + this_obs * *M * p]  * X[this_x_a + this_x_b * (i + *N * p)];
      s += theta[this_model + this_obs * (1 + *M * p)] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      else{
      for (p = 0; p < *P; p++) {
      s += theta[this_model + this_obs * (bin_high[this_i] + 1) + this_obs * *M * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      mn = breaks[bin_high[this_i]] - s * q_kappa[bin_high[this_i]];
      }
      p_ptr(&y_high[this_i], &mn, &s, shape, &U_high[this_i]);
      if (U_high[this_i] >= 1) {
      U_high[this_i] = 0.9999999;
      }
      }
      LAST_STATUS:
      if (status[this_i] == 3 || status[this_i] == 4) {
      if (low_thresh > y_low[this_i]) {
      log_d_ptr(&low_thresh, &mean_low, &sd_low, shape, &sig_low);
      sig_low = *tau_low * exp(-sig_low);
      if (*Pareto_low == 0) {
      U_low[this_i] =  *tau_low - *tau_low * (1 - exp(-(low_thresh - y_low[this_i]) / sig_low));
      } else{
      U_low[this_i] = (*tau_low - *tau_low * ( 1-  pow(1 + exp(*xi_low) * (low_thresh - y_low[this_i]) / sig_low, -1 / exp(*xi_low))));
      }
      if (U_low[this_i] <= 0) {
      U_low[this_i] = 0.0000001;
      }
      goto TAU_DIFF;
      }
      if (y_low[this_i] > high_thresh) {
      log_d_ptr(&high_thresh, &mean_high, &sd_high, shape, &sig_high);
      sig_high = (1 - *tau_high ) * exp(-sig_high);
      if (*Pareto_high == 0) {
      U_low[this_i] =    *tau_high + (1 - *tau_high) * (1 - exp(-(y_low[this_i] - high_thresh) / sig_high));
      } else{U_low[this_i] = *tau_high + (1 - *tau_high) * ( 1- pow(1 + exp(*xi_high) * (y_low[this_i] - high_thresh) / sig_high, -1 / exp(*xi_high)));}
      if (U_low[this_i] >= 1.0) {
      U_low[this_i] = 0.9999999;
      }
      goto TAU_DIFF;
      }
      for (m = 1; m < (*M - 1); m++) { /*element of breaks*/
      breaks[m] = 0;
      for (p = 0; p < *P; p++) { /*sum across all predictors*/
      breaks[m] +=   (BM[m + *M * p] + theta[this_model + 0 + this_obs * *M * p]) * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(M, breaks, &y_low[this_i], &bin_low[this_i]);
      s = 0;
      if (bin_low[this_i] == 0) {
      mn = 0;
      for (p = 0; p < *P; p++) {
      mn += theta[this_model + 0 + this_obs * *M * p]  * X[this_x_a + this_x_b * (i + *N * p)];
      s += theta[this_model + this_obs * (1 + *M * p)] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      } else{
      for (p = 0; p < *P; p++) {
      s += theta[this_model + this_obs * (bin_low[this_i] + 1) + this_obs * *M * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      mn = breaks[bin_low[this_i]] - s * q_kappa[bin_low[this_i]];
      }
      p_ptr(&y_low[this_i], &mn, &s, shape, &U_low[this_i]);
      if (U_low[this_i] <= 0) {
      U_low[this_i] = 0.0000001;
      }
      if (U_low[this_i] >= 1) {
      U_low[this_i] = 0.9999999;
      }
      }
      TAU_DIFF:
      log_like[this_i] =  log(U_high[this_i] - U_low[this_i]);
      U[this_i]      = (U_high[this_i] + U_low[this_i]) / 2;
      LL_SUM:
      *ll_sum += log_like[this_i];
      }
      Free(BM);
      Free(breaks);
      Free(mean_low_p);
      Free(sd_low_p);
      Free(mean_high_p);
      Free(sd_high_p);
      Free(q_kappa);
      Free(mn_1);
      Free(mn_2);
      Free(mn_3);
      
      Free(dummy_vec);
      Free(dummy_mat);
      }
      
      void clike_wrapper(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      int *basis_ind,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      double *w, double *M_low, double *M_high, double *I_low, double *I_high,
      int *g, int *j, int *h
      ) {
      
      void (*q_ptr)(double *, double *, double *, double *, double*);
      void (*p_ptr)(double *, double *, double *, double *, double*);
      void (*log_d_ptr)(double *, double *, double *, double *, double*);
      
      
      if (basis_ind[*h] == 0) {
      q_ptr               = &QNorm;
      p_ptr               = &PNorm;
      log_d_ptr           = &LogDNorm;
      }
      else if (basis_ind[*h] == 1) {
      q_ptr               = &QT;
      p_ptr               = &PT;
      log_d_ptr           = &LogDT;
      }
      else if (basis_ind[*h] == 2) {
      q_ptr               = &QLogistic;
      p_ptr               = &PLogistic;
      log_d_ptr           = &LogDLogistic;
      }
      else if (basis_ind[*h] == 3) {
      q_ptr               = &QAsymmetricLaplace;
      p_ptr               = &PAsymmetricLaplace;
      log_d_ptr           = &LogDAsymmetricLaplace;
      }
      else if (basis_ind[*h] == 4) {
      q_ptr               = &QWeibull;
      p_ptr               = &PWeibull;
      log_d_ptr           = &LogDWeibull;
      }
      else{
      q_ptr               = &QGamma;
      p_ptr               = &PGamma;
      log_d_ptr           = &LogDGamma;
      }
      
      double AAAA = 0.0;
      double *B  = (double *)Calloc(*M * (*M - 1), double); /*basis matrix*/
      make_b(q_ptr, M, kappa, &shape[*h], B);
      
      int u = *g + *G * *j + *G * *J * *h;
      ll_sum[u] = 0;
      clike(G, J, H, M, P, N, &n_gjh[u], status,
      kappa, X, y, y_low, y_high, theta, &shape[*h],
      q_ptr, p_ptr, log_d_ptr,
      &AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
      Pareto_low, Pareto_high, tau_low, tau_high, &xi_low[*g + *G * *h], &xi_high[*g + *G * *h],
      w, M_low, M_high, I_low, I_high,
      g, j, h);
      
      Free(B);
      }
      
      void quantile_function(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      void q_ptr(double *, double *, double *, double *, double*),
      void p_ptr(double *, double *, double *, double *, double*),
      void log_d_ptr(double *, double *, double *, double *, double*),
      double *AAAA,
      double *B,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      int *g, int *j, int *h
      ) {
      
      int i, k, m, p, this_i;
      int this_model = *g + *G * *j + *G * *J * *h;
      int this_obs = *G * *J * *H;
      int this_x_a = *g + *G * *j;
      int this_x_b = *G * *J;
      
      int M_1 = *M - 1;
      
      double mn, s;
      double dummy_mean = 0;
      double dummy_scale = 1;
      
      double *BM = (double *)Calloc((*M * *P), double);
      double breaks[*M]; /*breaks for each iteration*/
      
      double *q_kappa = (double *)Calloc((M_1), double);
      double *mean_low_p = (double *)Calloc(*P, double);
      double *sd_low_p = (double *)Calloc(*P, double);
      
      double *mean_high_p = (double *)Calloc(*P, double);
      double *sd_high_p = (double *)Calloc(*P, double);
      
      double sig_low, sig_high, low_thresh, high_thresh, mean_low, sd_low, mean_high, sd_high;
      
      double *w = (double *)Calloc((*M), double);
      double *I_low = (double *)Calloc((*M), double);
      double *I_high = (double *)Calloc((*M), double);
      double *dummy_vec = (double *)Calloc((M_1), double);
      double *dummy_mat = (double *)Calloc((M_1 * *P), double);
      
      double *mn_1 = (double *)Calloc(*P, double);
      double *mn_2 = (double *)Calloc(*P, double);
      double *mn_3 = (double *)Calloc(*P, double);
      
      
      for (m = 0; m < M_1; m++) {
      for (p = 0; p < *P; p++) {
      dummy_mat[m + M_1 * p] = theta[this_model + this_obs * ((m + 1) + *M * p)];
      }
      }
      
      for (m = 1; m < M_1; m++) {q_ptr(&kappa[m], &dummy_mean, &dummy_scale, shape, &q_kappa[m]);}
      
      if (*tau_low >  0) {
      I_low[0] = 1.0;
      m = 1;
      q_ptr(tau_low, &dummy_mean, &dummy_scale, &shape[*h], &I_low[m]);
      if (*M > 2) {
      for (m = 2; m < *M; m++) {
      I_low[m] = 0.0;
      }
      }
      for (p = 0; p < *P; p++) {
      mean_low_p[p] = theta[this_model + this_obs * (0 + *M * p)];
      sd_low_p[p]   = theta[this_model + this_obs * (1 + *M * p)];
      }
      }
      
      if (*tau_high < 1) {
      dummy_vec[0] = I_high[1];
      I_high[0] = 1;
      m = 1;
      if (*M == 2) {
      q_ptr(tau_high, &dummy_mean, &dummy_scale, &shape[*h], &I_high[m]);
      }
      if (*M > 2) {
      I_high[m] = q_kappa[m];
      for (m = 2; m < (*M - 1); m++) {
      I_high[m] = q_kappa[m] - q_kappa[m - 1];
      }
      q_ptr(tau_high, &dummy_mean, &dummy_scale, &shape[*h], &I_high[M_1]);
      I_high[M_1] -= q_kappa[(*M - 2)];
      }
      for (p = 0; p < *P; p++) {
      mn_1[p] = theta[this_model + this_obs * (0 + *M * p)];
      }
      if (*M == 2) {
      for (p = 0; p < *P; p++) {
      mean_high_p[p] = mn_1[p];
      }
      }
      if (*M > 2) {
      for (m = 0; m < (*M - 2); m++) {
      dummy_vec[m] = I_high[m + 1];
      }
      int T_A = 1;
      int T_B = 0;
      int dummy_one = 1;
      //      matrix_multiply(&T_A, &T_B, P, &dummy_one, &M_1, dummy_mat, dummy_vec, mn_2);
      matrix_multiply(&T_A, &T_B, &M_1, P, P, &dummy_one, dummy_mat, dummy_vec, mn_2);
      for (p = 0; p < *P; p++) {
      mean_high_p[p] = mn_1[p] + mn_2[p] - q_kappa[M_1 - 1] * theta[this_model + this_obs * ((*M - 1) + *M * p)];
      }
      }
      for (p = 0; p < *P; p++) {
      sd_high_p[p] = theta[this_model + this_obs * ((*M - 1) + *M * p)];
      }
      }
      /**********construct the break matrix************/
      for (i = 1; i < (*M - 1); i++) {
      for (p = 0; p < *P; p++) {
      BM[i + *M * p] = 0;
      for (k = 0; k < i; k++) {
      BM[i + *M * p] += B[i + *M * k] * theta[this_model + this_obs * ((k + 1) + *M * p)];
      }
      }
      }
      breaks[0] = -1000000.0;
      breaks[*M - 1] = 1000000.0;
      
      low_thresh  =  -100000000.0;
      high_thresh =  100000000.0;
      
      for (i = 0; i < *n_gjh; i++) {
      this_i = this_model + this_obs * i;
      if (status[this_i] == 0) { /*missing observations*/
      if (U[this_i] < *tau_low) {
      mean_low = X[this_x_a + this_x_b * (i + *N * 0)] * mean_low_p[0];
      sd_low = X[this_x_a + this_x_b * (i + *N * 0)] * sd_low_p[0];
      for (p = 1; p < *P; p++) {
      mean_low += X[this_x_a + this_x_b * (i + *N * p)] * mean_low_p[p];
      sd_low   += X[this_x_a + this_x_b * (i + *N * p)] * sd_low_p[p];
      }
      q_ptr(tau_low, &mean_low, &sd_low, shape, &low_thresh);
      log_d_ptr(&low_thresh, &mean_low, &sd_low, shape, &sig_low);
      sig_low = *tau_low * exp(-sig_low);
      if (*Pareto_low == 0) {
      y[this_i] = low_thresh + log(sig_low) * log(U[this_i] / *tau_low);
      } else {
      y[this_i] = low_thresh - (sig_low / exp(*xi_low)) *  (pow((U[this_i] / *tau_low), -exp(*xi_low)) - 1);
      }
      } else if (*tau_high < U[this_i]) {
      mean_high = X[this_x_a + this_x_b * (i + *N * 0)] * mean_high_p[0];
      sd_high = X[this_x_a + this_x_b * (i + *N * 0)] * sd_high_p[0];
      for (p = 1; p < *P; p++) {
      mean_high += X[this_x_a + this_x_b * (i + *N * p)] * mean_high_p[p];
      sd_high   += X[this_x_a + this_x_b * (i + *N * p)] * sd_high_p[p];
      }
      q_ptr(tau_high, &mean_high, &sd_high, shape, &high_thresh);
      log_d_ptr(&high_thresh, &mean_high, &sd_high, shape, &sig_high);
      sig_high = (1 - *tau_high ) * exp(-sig_high);
      if (*Pareto_high == 0) {
      y[this_i] = high_thresh - sig_high * log((1 - U[this_i])/(1 - *tau_high));
      } else{
      y[this_i] = high_thresh  + (sig_high / exp(*xi_high)) * (pow((1 - U[this_i])/(1 - *tau_high), -exp(*xi_high)) - 1);
      }
      } else{
      for (m = 1; m < (*M - 1); m++) { /*element of breaks*/
      breaks[m] = 0;
      for (p = 0; p < *P; p++) { /*sum across all predictors*/
      breaks[m] +=   (BM[m + *M * p] + theta[this_model + 0 + this_obs * *M * p]) * X[this_x_a + this_x_b * (i + *N * p)];
      }
      }
      find_interval(M, kappa, &U[this_i], &bin[this_i]);
      s = 0;
      if (bin[this_i] == 0) {
      mn = 0;
      for (p = 0; p < *P; p++) {
      mn += theta[this_model + 0 + this_obs * *M * p]  * X[this_x_a + this_x_b * (i + *N * p)];
      s += theta[this_model + this_obs * (1 + *M * p)] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      } else{
      for (p = 0; p < *P; p++) {
      s += theta[this_model + this_obs * (bin[this_i] + 1) + this_obs * *M * p] * X[this_x_a + this_x_b * (i + *N * p)];
      }
      mn = breaks[bin[this_i]] - s * q_kappa[bin[this_i]];
      }
      q_ptr(&U[this_i], &mn, &s, shape, &y[this_i]);
      }
      }
      }
      
      Free(q_kappa);
      Free(mean_low_p);
      Free(sd_low_p);
      Free(mean_high_p);
      Free(sd_high_p);
      
      Free(w);
      Free(I_low);
      Free(I_high);
      Free(dummy_vec);
      Free(dummy_mat);
      
      Free(mn_1);
      Free(mn_2);
      Free(mn_3);
      
      }
      
      void quantile_function_wrapper(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      int *basis_ind,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      int *g, int *j, int *h
      ) {
      void (*q_ptr)(double *, double *, double *, double *, double*);
      void (*p_ptr)(double *, double *, double *, double *, double*);
      void (*log_d_ptr)(double *, double *, double *, double *, double*);
      
      
      if (basis_ind[*h] < 1) {
      q_ptr               = &QNorm;
      p_ptr               = &PNorm;
      log_d_ptr           = &LogDNorm;
      }
      else if (basis_ind[*h] == 1) {
      q_ptr               = &QT;
      p_ptr               = &PT;
      log_d_ptr           = &LogDT;
      }
      else if (basis_ind[*h] == 2) {
      q_ptr               = &QLogistic;
      p_ptr               = &PLogistic;
      log_d_ptr           = &LogDLogistic;
      }
      else if (basis_ind[*h] == 3) {
      q_ptr               = &QAsymmetricLaplace;
      p_ptr               = &PAsymmetricLaplace;
      log_d_ptr           = &LogDAsymmetricLaplace;
      }
      else if (basis_ind[*h] == 4) {
      q_ptr               = &QWeibull;
      p_ptr               = &PWeibull;
      log_d_ptr           = &LogDWeibull;
      }
      else{
      q_ptr               = &QGamma;
      p_ptr               = &PGamma;
      log_d_ptr           = &LogDGamma;
      }
      
      double AAAA = 0.0;
      
      double *B  = (double *)Calloc(*M * (*M - 1), double); /*basis matrix*/
      
      make_b(q_ptr, M, kappa, &shape[*h], B);
      
      int u = *g + *G * *j + *G * *J * *h;
      ll_sum[u] = 0;
      quantile_function(G, J, H, M, P, N, &n_gjh[u], status,
      kappa, X, y, y_low, y_high, theta, &shape[*h],
      q_ptr, p_ptr, log_d_ptr,
      &AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
      Pareto_low, Pareto_high, tau_low, tau_high, &xi_low[*g + *G * *h], &xi_high[*g + *G * *h],
      g, j, h);
      Free(B);
      }
      
      
      void quantile_function_spline(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      void q_ptr(double *, double *, double *, double *, double*),
      void p_ptr(double *, double *, double *, double *, double*),
      void log_d_ptr(double *, double *, double *, double *, double*),
      double *AAAA,
      double *B,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      int *g, int *j, int *h
      ) {
      
      double z, a, b, c, d;
      double dummy_low = 0.0;
      double dummy_high = 0.0;
      int i, k, m, p, this_i;
      int this_model = *g + *G * *j + *G * *J * *h;
      int this_obs = *G * *J * *H;
      int this_x_a = *g + *G * *j;
      int this_x_b = *G * *J;
      
      int M_1 = *M - 1;
      int kappa_length = *M + 4;
      int dummy_bin = 0;
      
      double *BM = (double *)Calloc((kappa_length * *P), double);
      double *breaks = (double *)Calloc(kappa_length, double);
      
      double *I_low_p = (double *)Calloc(*P, double);
      double *m_low_p = (double *)Calloc(*P, double);
      double *I_high_p = (double *)Calloc(*P, double);
      double *m_high_p = (double *)Calloc(*P, double);
      
      double *w = (double *)Calloc((*M), double);
      double *M_low = (double *)Calloc((*M), double);
      double *M_high = (double *)Calloc((*M), double);
      double *I_low = (double *)Calloc((*M), double);
      double *I_high = (double *)Calloc((*M), double);
      double sig_low, sig_high, low_thresh, high_thresh;
      
      M_low[0] = 0.0;
      M_high[0] = 0.0;
      I_low[0] = 1.0;
      I_high[0] = 1.0;
      
      for (m = 0; m < M_1; m++) {
      CubicMSpline(tau_low,   &kappa_length, kappa, &m,  &dummy_bin, &M_low[m + 1]);
      M_low[m + 1] *= *tau_low;
      CubicMSpline(tau_high,   &kappa_length, kappa, &m, &dummy_bin, &M_high[m + 1]);
      M_high[m + 1] *= (1 - *tau_high);
      CubicISpline(tau_low,   &kappa_length, kappa, &m, &dummy_bin, &I_low[m + 1]);
      CubicISpline(tau_high,  &kappa_length, kappa, &m, &dummy_bin, &I_high[m + 1]);
      }
      for (p = 0; p < *P; p++) {
      I_low_p[p] = 0.0;
      m_low_p[p] = 0.0;
      I_high_p[p] = 0.0;
      m_high_p[p] = 0.0;
      for (m = 0; m < *M; m++) {
      k = this_model + this_obs * (m + *M * p);
      m_low_p[p]    += theta[k] * M_low[m];
      I_low_p[p]      += theta[k] * I_low[m];
      m_high_p[p]   += theta[k] * M_high[m];
      I_high_p[p]     += theta[k] * I_high[m];
      }
      }
      for (i = 0; i < kappa_length; i++) {
      for (p = 0; p < *P; p++) {
      BM[i + kappa_length * p] = 0.0;
      for (m = 0; m < *M; m++) {
      BM[i + kappa_length * p] += B[i + kappa_length * m] * theta[this_model + this_obs * (m + *M * p)];
      }
      }
      }
      for (i = 0; i < *n_gjh; i++) {
      this_i = this_model + this_obs * i;
      if (status[this_i] == 0) {
      if (*tau_low > U[this_i]) {
      low_thresh = I_low_p[0] * X[this_x_a + this_x_b * i];
      sig_low = m_low_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      low_thresh += X[this_x_a + this_x_b * (i + *N * p)] * I_low_p[p];
      sig_low += X[this_x_a + this_x_b * (i + *N * p)] * m_low_p[p];
      }
      if (*Pareto_low == 0) {
      y[this_i] = low_thresh + log(sig_low) * log(U[this_i] / *tau_low);
      } else {
      y[this_i] = low_thresh - (sig_low / exp(*xi_low)) *  (pow((U[this_i] / *tau_low), -exp(*xi_low)) - 1);
      }
      } else if (*tau_high < U[this_i]) {
      high_thresh = I_high_p[0] * X[this_x_a + this_x_b * i];
      sig_high = m_high_p[0] * X[this_x_a + this_x_b * i];
      for (p = 1; p < *P; p++) {
      high_thresh += X[this_x_a + this_x_b * (i + *N * p)] * I_high_p[p]; sig_high += X[this_x_a + this_x_b * (i + *N * p)] * m_high_p[p];
      }
      if (*Pareto_high == 0) {
      y[this_i] = high_thresh - sig_high * log((1 - U[this_i])/(1 - *tau_high));
      } else {
      y[this_i] = high_thresh  + (sig_high / exp(*xi_high)) * (pow((1 - U[this_i])/(1 - *tau_high), -exp(*xi_high)) - 1);
      }
      } else{
      for (m = 0; m < *M; m++) {
      w[m] = 0.0;
      dummy_low = 0.0;
      dummy_high = 0.0;
      for (p = 0; p < *P; p++) {
      z = X[this_x_a + this_x_b * (i + *N * p)];
      w[m]        += theta[this_model + this_obs * (m + *M * p)] * z;
      dummy_low   += BM[ bin[this_i]      + kappa_length * p] * z;
      dummy_high  += BM[(bin[this_i] + 1) + kappa_length * p] * z;
      }
      }
      find_interval(&kappa_length, kappa, &U[this_i], &bin[this_i]);
      a = 0.0;
      b = 0.0;
      c = 0.0;
      d = w[0];
      for (m = 0; m < M_1; m++) {
      k = 4 * (m + M_1 * bin[this_i]);
      a += AAAA[0 + k] * w[m + 1];
      b += AAAA[1 + k] * w[m + 1];
      c += AAAA[2 + k] * w[m + 1];
      d += AAAA[3 + k] * w[m + 1];
      }
      y[this_i]           =     a * pow(U[this_i], 3) +     b * pow(U[this_i], 2) + c * U[this_i] + d;
      }
      }
      }
      Free(BM);
      Free(breaks);
      Free(w);
      Free(M_low);
      Free(M_high);
      Free(I_low);
      Free(I_high);
      }
      
      void quantile_function_spline_wrapper(int *G, int *J, int *H, int *M, int *P,
      int *N, int *n_gjh, int *status, double *kappa,
      double *X, double *y, double *y_low, double *y_high,
      double *theta,
      double *shape,
      int *bin, int *bin_low, int *bin_high, double *ll_sum, double *log_like, double *U, double *U_low, double *U_high,
      int *Pareto_low, int *Pareto_high, double *tau_low, double *tau_high, double *xi_low, double *xi_high,
      int *g, int *j, int *h
      ) {
      int M_1 = *M - 1;
      int m;
      
      double *AAAA = (double *)Calloc((4 * *M  * M_1), double);
      double *B = (double *)Calloc(((*M + 4)  * *M), double);
      
      void (*q_ptr)(double *, double *, double *, double *, double*);
      void (*p_ptr)(double *, double *, double *, double *, double*);
      void (*log_d_ptr)(double *, double *, double *, double *, double*);
      
      for (m = 0; m < 4 * *M * M_1; m++) {
      AAAA[m] = 0.0;
      }
      MakeCubicSplineCoefficients(M, kappa, AAAA);
      
      q_ptr               = &QNorm;
      p_ptr               = &PNorm;
      log_d_ptr           = &LogDNorm;
      
      for (m = 0; m < (*M + 4) * *M; m++) {
      B[m] = 0.0;
      }
      make_b(q_ptr, M, kappa, &shape[*h], B);
      
      int u = *g + *G * *j + *G * *J * *h;
      int tail_index = *g + *G * *h;
      
      shape[*h] = -1.0;
      
      quantile_function_spline(G, J, H, M, P,
      N, &n_gjh[u], status, kappa,
      X, y, y_low, y_high, theta, &shape[*h],
      q_ptr, p_ptr, log_d_ptr, AAAA,
      B,
      bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
      Pareto_low, Pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
      g, j, h);
      
      Free(AAAA);
      Free(B);
      }
      
      
      /********  regression functions *********/
      /*function that for mth basis function and pth predictor projects polynomial onto theta and ensures mth row is proper*/
      void threshold(int *G, int *J, int *D_vec, int *D, int *H, int *M, int *P, int *g, int *h, int *m, int *p, double *timepoints, double *thetastar, double *theta) {
      int d, j, k, l, p_int;
      double neg_part;
      double *X_D = (double *)Calloc((*J * *D), double);
      int this_theta;
      int theta_index = *G * *J;
      int thetastar_index = *G * *D;
      
      for (p_int = 0; p_int < *P; p_int++) {
      this_theta = *h + *H * *m + *H * *M * p_int;
      /*identity mapping*/
      if (D_vec[*p] == 1) {
      k = *G * *J * this_theta;
      l = *G * *D * this_theta;
      for (j = 0; j < *J; j++) {
      theta[*g + *G * j + k] = thetastar[*g + l];
      }
      }
      else{
      for (j = 0; j < *J; j++) {
      X_D[j] = 1.0;
      X_D[j + *J * 1] = timepoints[j];
      }
      for (j = 0; j < *J; j++) {
      k = *g + *G * j + theta_index * this_theta;
      theta[k] = 0;
      for (d = 0; d < D_vec[*p]; d++) {
      theta[k] += X_D[j + *J * d] * thetastar[*g + *G * d + thetastar_index * this_theta];
      }
      }
      }
      if (*m > 0) {
      for (j = 0; j < *J; j++) {
      k = *g + *G * j + theta_index * (*h + *H * *m); /*p = 0 case*/
      if (theta[k] <= 0) {
      theta[k] = 0.001;
      }
      neg_part = theta[k];
      for (l = 1; l < *P; l++) {
      k = *g + *G * j + theta_index * (*h + *H * *m + *H * *M * l);
      if (theta[k] < 0) {
      neg_part +=  theta[k];
      }
      else{
      neg_part -= theta[k];
      }
      }
      if (neg_part < 0) {
      for (l = 1; l < *P; l++) {
      theta[*g + *G * j + theta_index * (*h + *H * *m + *H * *M * l)] = 0;
      }
      }
      }
      }
      }
      Free(X_D);
      }
      
      /********* prior functions ***********/
      void gibbs_mu(int *G, int *D, int *H, int *M, int *P, int *d, int *m, int *p,
      double *thetastar, double *prec_g, double *prec_h, double *z_mu,
      double *mu_0, double *prec_mu_0, double *mu, double *mu_array
      ) {
      
      double *dummy_mean_mu       = (double *)Calloc(*H, double);
      double *dummy_prec_h        = (double *)Calloc(*H * *H, double);
      double *col_sums            = (double *)Calloc(*G, double);
      double *this_theta          = (double *)Calloc(*G * *H, double);
      double *cov_mu              = (double *)Calloc(*H * *H, double);
      double *dummy_matrix        = (double *)Calloc(*G * *G * *H, double);
      double *mean_mu             = (double *)Calloc(*H, double);
      double *chol_mu             = (double *)Calloc(*H * *H, double);
      
      int prior_index = *d + *D * *m + *D * *M * *p;
      int DMP = *D * *M * *P;
      int thetastar_index = *m + *M * *p;
      int this_thetastar = *G * *D * *H * thetastar_index;
      int g, h, k, k1, k2;
      int dummy_one = 1;
      int T_A = 0;
      int T_B = 0;
      int dim_ll = *G * *H;
      
      double dummy_log_det = 0.0;
      
      double sumprec_g = 0;
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *G; h++) {
      sumprec_g += prec_g[prior_index + DMP * (g + *G * h)];
      }
      }
      /*sumprec_g is right*/
      for (g = 0; g < *G; g++) {
      col_sums[g] = 0;
      for (h = 0; h < *G; h++) {
      col_sums[g] += prec_g[prior_index + DMP * (h + *G * g)]; /*keep the column fixed*/
      }
      }
      /*col_sums are right*/
      for (g = 0; g < *H; g++) {
      for (h = 0; h < *H; h++) {
      k = g + *H * h;
      dummy_prec_h[k] = prec_h[prior_index + DMP * k];
      chol_mu[k]  = sumprec_g * dummy_prec_h[k];
      }
      }
      /*precision matrix is correct*/
      for (h = 0; h < *H; h++) {
      chol_mu[h + *H * h] += prec_mu_0[prior_index];
      }
      /*invert the precision*/
      InvertSymmetricMatrix(H, chol_mu, cov_mu, &dummy_log_det);
      /*perform first matrix mulitiplication for the mean*/
      Kronecker(&dummy_one, G, H, H, col_sums, dummy_prec_h, dummy_matrix);
      
      /*find the correct theta*/
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      k1 = g + *G * *d + *G * *D * h;
      k2 = h + *H * g;
      this_theta[k2] = thetastar[k1 + this_thetastar];
      }
      }
      /*multiply the matrix by theta*/
      matrix_multiply(&T_A, &T_B, H, &dim_ll, &dim_ll, &dummy_one, dummy_matrix, this_theta, dummy_mean_mu);
      //  matrix_multiply(&T_A, &T_B, H, &dummy_one, &dim_ll, dummy_matrix, this_theta, dummy_mean_mu);
      
      /*add the precision*/
      for (h = 0; h < *H; h++) {
      dummy_mean_mu[h] += mu_0[prior_index];
      }
      /*multiply the location by the covariance matrix cov_mu*/
      matrix_multiply(&T_A, &T_B, H, H, H, &dummy_one, cov_mu, dummy_mean_mu, mean_mu);
      //  matrix_multiply(&T_A, &T_B, H, &dummy_one, H, cov_mu, dummy_mean_mu, mean_mu);
      /*cholesky the covariance matrix*/
      CholMatrix(H, cov_mu, chol_mu, &dummy_log_det);
      
      for (h = 0; h < *H; h++) {
      k = prior_index + DMP * h;
      mu[k] = mean_mu[h];
      for (g = 0; g < (h + 1); g++) {
      mu[k] += chol_mu[h + *H * g] * z_mu[g];
      }
      }
      for (h = 0; h < *H; h++) {
      k = *G * *d + *G * *D * h + this_thetastar;
      for (g = 0; g < *G; g++) {
      mu_array[g + k] = mu[prior_index + DMP * h];
      }
      }
      
      Free(dummy_mean_mu);
      Free(dummy_prec_h);
      Free(col_sums);
      Free(this_theta);
      Free(cov_mu);
      Free(dummy_matrix);
      Free(mean_mu);
      Free(chol_mu);
      }
      
      void make_sai(int *G, int *D, int *H, int *M, int *P, int *d, int *m, int *p,
      double *thetastar, double *prec_g, double *mu_array, double *sai_h, double *sai_h_0
      ) {
      
      double *g1_vec       = (double *)Calloc(*G, double);
      double *g2_vec       = (double *)Calloc(*G, double);
      double *dummy_prec_g = (double *)Calloc(*G * *G, double);
      
      int T_A = 0;
      int T_B = 1;
      int dummy_one = 1;
      int this_sai;
      int g, h, k1, k2, u, v;
      
      int prior_index = *d + *D * *m + *D * *M * *p;
      int DMP         = *D * *M * *P;
      int this_thetastar = *G * *D * *H * (*m + *M * *p);
      
      for (u = 0; u < *H * *H; u++) {sai_h[u] = 0;}
      
      for (u = 0; u < *H; u++) {
      for (v = 0; v < u + 1; v++) {
      this_sai = u + *H * v;
      for (g = 0; g < *G; g++) {
      k1 = g + *G * (*d + *D * u) + this_thetastar;
      k2 = g + *G * (*d + *D * v) + this_thetastar;
      g1_vec[g] = thetastar[k1] - mu_array[k1];
      g2_vec[g] = thetastar[k2] - mu_array[k2];
      }
      //      matrix_multiply(&T_A, &T_B, G, G, &dummy_one, g1_vec, g2_vec, dummy_prec_g);
      matrix_multiply(&T_A, &T_B, G, &dummy_one, G, &dummy_one, g1_vec, g2_vec, dummy_prec_g);
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *G; h++) {
      sai_h[this_sai] += dummy_prec_g[g + *G * h] * prec_g[prior_index + DMP * (g + *G * h)];
      }
      }
      sai_h[this_sai] += sai_h_0[prior_index + DMP * this_sai];
      }
      }
      for (u = 0; u < *H; u++) {
      for (v = u + 1; v < *H; v++) {
      sai_h[u + *H * v] = sai_h[v + *H * u];
      }
      }
      Free(g1_vec);
      Free(g2_vec);
      Free(dummy_prec_g);
      }
      
      /*function that updates the precision array omega*/
      void make_omega(int *G, int *D, int *H, int *M, int *P, int *d, int *m, int *p,
      double *prec_g, double *prec_h, double *omega) {
      int g, h, k;
      double *omega_g = (double *)Calloc(*G * *G, double);
      double *omega_h = (double *)Calloc(*H * *H, double);
      double *dummy_omega = (double *)Calloc( (*G * *H) * (*G * *H), double);
      
      int this_omega =  *d + *D * *m + *D * *M * *p;
      int omega_index = *D * *M * *P;
      
      if (*G > 1) {
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *G; h++) {
      k = g + *G * h;
      omega_g[k] = prec_g[this_omega + omega_index * k];
      }
      }
      for (g = 0; g < *H; g++) {
      for (h = 0; h < *H; h++) {
      k = g + *H * h;
      omega_h[k] = prec_h[this_omega + omega_index * k];
      }
      }
      Kronecker (G, G, H, H, omega_g, omega_h, dummy_omega);
      for (g = 0; g < *G * *H; g++) {
      for (h = 0; h < *G * *H; h++) {
      k = g + *G * *H * h;
      omega[this_omega + omega_index * k] = dummy_omega[k];
      }
      }
      }
      if (*G == 1) {
      omega[this_omega] = prec_h[this_omega];
      }
      Free(omega_g);
      Free(omega_h);
      Free(dummy_omega);
      }
      
      /*computes (theta-1*mu)%*%AR(PREC)(theta-1*mu) where mu is a scalar*/
      /*easy to extend later to mu being a function of g*/
      /*OMEGA is D x M x P x G*H x G*H */
      void rho_quad(int *G, int *D, int *H, int *M, int *P, int *d, int *m, int *p, double *omega, double *thetastar, double *mu_array, double *rhoquad) {
      int g, h, k;
      double *x = (double *)Calloc((*G * *H), double);
      int this_omega =  *d + *D * *m + *D * *M * *p;
      int dim_omega  =  *D * *M * *P;
      int this_thetastar = *G * *D * *H * (*m + *M * *p);
      int GH = *G * *H;
      double dummy_scalar;
      *rhoquad = 0.0;
      
      if (*G > 1) {
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      k = g + *G * *d  + *G * *D * h + this_thetastar;
      x[h + *H * g] = thetastar[k] - mu_array[k];
      }
      }
      for (g = 0; g < GH; g++) {
      for (h = 0; h < GH; h++) {
      *rhoquad += x[g] * x[h] * omega[this_omega + dim_omega * (g + GH * h)];
      }
      }
      }
      if (*G == 1) {
      for (h = 0; h < *H; h++) {
      k = *d  + *G * *D * h + this_thetastar;
      dummy_scalar = thetastar[k] - mu_array[k];
      *rhoquad += dummy_scalar * dummy_scalar;
      }
      *rhoquad *= omega[this_omega];
      }
      Free(x);
      }
      
      /****** copula functions ******/
      /*****************************
      function that calculates the log density of N MVN random variables of dimension J with mean 0 and covariance Sigma
      V is N x J matrix of Gaussian RVs
      Z is N x J x Q matrix of predictors
      Sigma is N x J x J array of covariance matrices
      sig_alpha is Q x 1 array of random effect precisions
      sig_eps is 1 x 1 array of random error precision
      ********************************/
      
      //cop_ll_sum should be of length G
      void copula_ll(int *G, int *J, int *H, int *Q_int, int *n_g, int *N,
      int *g,
      double *W2, double *Z, double *eta, double *DD, double *DD_inv, double *gamma, double *copula_loglike, double *cop_ll_sum) {
      
      int h, j, k, n;
      int JH = *J * *H;
      int HQ = *H * *Q_int;
      double qf, marginal_qf;
      double log_det;
      double *x              = (double *)Calloc(JH, double);
      
      *cop_ll_sum = 0;
      
      if (n_g[*g] > 0) {
      for (n = 0; n < n_g[*g]; n++) {
      qf = 0.0;
      marginal_qf = 0.0;
      log_det = 0.0;
      for (j = 0; j < JH; j++) {
      x[j] = 0;
      k = *g + *G * j + *G * JH * n;
      for (h = 0; h < HQ; h++) {
      x[j] += Z[*g + *G * (j +  JH * n + JH * *N * h)] * gamma[*g + *G * (n + *N * h)];
      }
      x[j] += eta[*g + *G * (j + JH * n)];
      x[j] *= -DD[k];
      x[j] +=  W2[k];
      qf += x[j] * x[j] * DD_inv[k] * DD_inv[k];
      log_det -= log(DD[k]);
      marginal_qf += W2[k] * W2[k];
      }
      /*include a factor of 2 in log_det because variance of W2 is D^2*/
      copula_loglike[*g + *G * n] = 0.5 * (2 * log_det - qf + marginal_qf);
      *cop_ll_sum += copula_loglike[*g + *G * n];
      }
      }
      Free(x);
      }
      
      void make_D(int *G, int *J, int *H, int *Q_int, int *n_g, int *N,
      double *Z, double *delta, double *W_cov, int *g, int *status, double *DD_inv, double *DD) {
      
      int j, h, k, n;
      int JH = *J * *H;
      int HQ = *H * *Q_int;
      double zzz = 0;
      
      if (n_g[*g] > 0) {
      for (n = 0; n < n_g[*g]; n++) {
      for (j = 0; j < JH; j++) {
      k = *g + *G * (j + JH * n);
      DD_inv[k] = 0.0;
      if (status[k] == -99 /*dead people*/) {
      DD[k] = 1.0;
      } else{
      for (h = 0; h < HQ; h++) {
      zzz = Z[*g + *G * (j +  JH * n + JH * *N * h)];
      DD_inv[k] += zzz * zzz * delta[*g + *G * h];
      }
      DD_inv[k] += W_cov[*g + *G * (j + JH * j)];
      DD_inv[k] = sqrt(DD_inv[k]);
      DD[k]     = pow(DD_inv[k], -1); /*need to see if I divide by 0 for dead/missing values here*/
      }
      }
      }
      }
      }
      
      void longitudinal_log_det(int *J, double *AR, double *AR_chol, double *log_det_matrix) {
      char uplo = 'L';
      int j, u, v, k1, k2, info;
      for (j = 0; j < *J; j++) {
      for (u = 0; u < *J; u++) {
      k1 = (*J - j + u) % *J;
      for (v = 0; v < *J; v++) {
      k2 = (*J - j + v) % *J;
      AR_chol[k1 + *J * k2] = AR[u + *J * v];
      }
      }
      F77_CALL(dpotrf)(&uplo, J, AR_chol, J, &info);
      if (info) {
      Rprintf("Error with chol(mat_chol): info = %d\n", info);
      }
      for (u = 0; u < *J; u++) {
      log_det_matrix[j + *J * u] = -2 * log(AR_chol[u + *J * u]);
      }
      }
      }
      
      /********* tail functions *********/
      /*this function draws xi from the prior*/
      void update_xi_exponential(int *G, int *H, double *xi_prec_g, double *xi_sig_h, double *xi_mu_array, double *xi) {
      int g, h;
      int GH = *G * *H;
      int G2 = *G * *G;
      int H2 = *H * *H;
      double dummy_log_det = 0.0;
      double *dummy_sig_g    = (double *)Calloc(G2, double);
      double *dummy_chol_g   = (double *)Calloc(G2, double);
      double *dummy_chol_h   = (double *)Calloc(H2, double);
      double *dummy_chol_gh  = (double *)Calloc(GH * GH, double);
      double *z              = (double *)Calloc(GH, double);
      
      for (g = 0; g < G2; g++) {
      dummy_sig_g[g] = 0.0;
      dummy_chol_g[g] = 0.0;
      }
      for (h = 0; h < H2; h++) {
      dummy_chol_h[h] = 0.0;
      }
      for (g = 0; g < GH * GH; g++) {
      dummy_chol_gh[g] = 0.0;
      }
      for (g = 0; g < GH; g++) {
      z[g] = rnorm(0,1);
      }
      /*invert the G x G precision matrix*/
      InvertSymmetricMatrix(G, xi_prec_g, dummy_sig_g, &dummy_log_det);
      /*cholesky the G x G covariance matrix*/
      CholMatrix(G, dummy_sig_g, dummy_chol_g, &dummy_log_det);
      /*cholesky the H x H covariance matrix*/
      CholMatrix(H, xi_sig_h, dummy_chol_h, &dummy_log_det);
      /*take the Kronecker product of the choleskys*/
      Kronecker(G, G, H, H,  dummy_chol_g, dummy_chol_h, dummy_chol_gh);
      /*update xi*/
      for (g = 0; g < GH; g++) {
      xi[g] = xi_mu_array[g];
      for (h = 0; h < (g + 1); h++) {
      xi[g] += dummy_chol_gh[g + GH * h] * z[h];
      }
      }
      Free(dummy_sig_g);
      Free(dummy_chol_g);
      Free(dummy_chol_h);
      Free(dummy_chol_gh);
      Free(z);
      }
      
      
      
      /********* data manipulation functions ***************/
      
      void make_array(int *G, int *J, int *H, int *N, int *P, int *Q_int, int *n_g, int *n_g_cumsum, int *N_row,
      double *qreg_Y, double *qreg_Y_low, double *qreg_Y_high, double *qreg_X, double *qreg_Z, int *qreg_status,
      double *Y, double *Y_low, double *Y_high, double *X, double *Z_red, double *Z, int *status) {
      
      int g, h, j, k1, k2, n, p, q, n_obs;
      int GJ = *G * *J;
      int GJH = *G * *J * *H;
      int GJHN = *G * *J * *H * *N;
      double *overall_mean = (double *)Calloc(*H, double);
      double h_mean;
      
      for (h = 0; h < *H; h++) {
      overall_mean[h] = 0.0;
      n_obs = 0;
      for (g = 0; g < *G; g++) {
      for (j = 0; j < *J; j++) {
      for (n = 0; n < n_g[g]; n++) {
      k1 = g + *G * j + *G * *J * (h + *H * n);
      k2 = (n_g_cumsum[g] + n) * *J + j + *N_row * h;
      if (qreg_status[k2] > 0) {
      overall_mean[h] += Y[k1];
      n_obs ++;
      }
      }
      }
      }
      overall_mean[h] /= n_obs;
      }
      for (g = 0; g < *G; g++) {
      for (n = 0; n < n_g[g + 1]; n++) {
      for (h = 0; h < *H; h++) {
      h_mean = 0.0;
      n_obs = 0;
      for (j = 0; j < *J; j++) {
      k1 = g + *G * j + GJ * (h + *H * n);
      k2 = (n_g_cumsum[g] + n) * *J + j + *N_row * h;
      Y[k1] = qreg_Y[k2];
      Y_low[k1] = qreg_Y_low[k2];
      Y_high[k1] = qreg_Y_high[k2];
      status[k1] = qreg_status[k2];
      if (qreg_status[k2] > 0) {
      h_mean += Y[k1];
      n_obs ++;
      }
      for (p = 0; p < *P; p++) {
      X[g + *G * j + GJ * (n + *N * p)] = qreg_X[(n_g_cumsum[g] + n) * *J + j + *N_row * p];
      }
      for (q = 0; q < *Q_int; q++) {
      Z_red[g + *G * j + GJ * (n + *N * q)] = qreg_Z[(n_g_cumsum[g] + n) * *J + j + *N_row * q];
      }
      }
      /*if all observations are missing for an individual use population mean for initial value*/
      if (n_obs == 0) {
      for (j = 0; j < *J; j++) {
      Y[g + *G * j + GJ * (h + *H * n)] = overall_mean[h];
      }
      } else if (n_obs < *J) {
      /*if some observations are present for an individual use individual mean for initial value*/
      h_mean /= n_obs;
      for (j = 0; j < *J; j++) {
      k2 = (n_g_cumsum[g] + n) * *J + j + *N_row * h;
      if (qreg_status[k2] < 1) {
      Y[g + *G * j + GJ * (h + *H * n)] = h_mean;
      }
      }
      }
      }
      for (h = 0; h < *H; h++) {
      for (j = 0; j < *J; j++) {
      for (q = 0; q < *Q_int; q++) {
      Z[g + *G * (h * *J + j) + GJH * n + GJHN * (h * *Q_int + q)] = Z_red[g + *G * j + GJ * (n + *N * q)];
      }
      }
      }
      }
      }
      Free(overall_mean);
      }
      
      /*
      n1 is the number of observations that were censored below
      n2 - n1 is the number of uncensored observations
      n - n2 is the number of observations that were censored above
      
      P is the number of predictors
      P1 is the subset of predictors that affect the shape
      M is the number of basis functions
      
      X is the design matrix (n) x P
      y is the observed data
      
      beta are the location parameters (P x 1 vector)
      alpha are the scale parameters (M x P matrix)
      kappa are the breakpoints
      shape is the shape parameter for the basis functions
      
      q_ptr is the quantile function
      log_d_ptr is the log density function
      
      B is the matrix of basis functions
      
      ll_sum is the summed log likelihood
      bin is the vector of length N containing which bin each observation is in
      
      BM is the break matrix
      */
      void MCMCOld(int *MCMC_ints, int *basis_ind, double *kappa,
      int *G, int *J, int *D_vec, int *D, int *H, int *M, int *P, int *P1, int *Q_int,
      int *N, int *n_gjh, int *status,
      double *X, double *Y_array,
      
      double *thetastar,
      double *mu_0, double *prec_mu_0,
      double *nu_h, double *sai_h_0,
      double *rho_a, double *rho_b, int *correlation, double *dist,
      
      double *shape_prior,
      
      int *cop, double *Z,
      double *sai_lambda_0, double *nu_lambda, double *zeta, double *dist_J, double *timepoints,
      
      double *pareto_probs, double *tau_low, double *tau_high,
      double *xi_mu_0,  double *xi_prec_mu_0,  double *xi_sai_0,  double *xi_nu,
      
      double *THETASTAR, double *MU, double *SIGMA, double *RHO,
      double *SHAPE,
      double *MIC,
      double *post_q,
      double *tau_q,
      int *N_tau_q,
      double *GAMMA, double *ETA, double *ALPHA, double *LAMBDA, double *DELTA,
      int *PARETO, double *XI,  double *XI_MU,  double *XI_SIG,  double *XI_RHO) {
      
      
      ///***************** Part I: allocate memory, initialize stuff *********************/
      int G2             = *G * *G;
      int H2             = *H * *H;
      int HQ             = *H * *Q_int;
      int JH             = *J * *H;
      int J2             = *J * *J;
      int dim_Z2_sum     = *G * *N * HQ * HQ;  if (*cop == 0) {dim_Z2_sum = 1;}
      int GJ             = *G * *J;
      int GH             = *G * *H;
      int GJH            = *G * *J * *H;
      int DMP            = *D * *M * *P;
      int GDH            = *G * *D * *H;
      int GDHP           = *G * *D * *H * *P;
      
      int dim_ll         = *G * *J * *H;
      int dim_thetastar  = *G * *D * *H * *M * *P;
      int dim_theta      = *G * *J * *H * *M * *P;
      int dim_Y          = *G * *J * *H * *N;
      int dim_omega      =           *D * *M * *P * (*G  * *H)  * (*G  * *H);
      int dim_prec_g     =           *D * *M * *P * *G * *G;
      int dim_prec_h     =           *D * *M * *P * *H * *H;
      int dim_xi         = *G * *H;
      int dim_eta_cov = *G * JH * JH;
      int dim_this_AR = *G * *J * *J;
      int dim_gamma   = *G * *N * HQ;
      int dim_eta     = *G * *N * JH;
      int dim_lambda  = *G * *H * *H;
      int dim_delta   = *G * HQ;
      
      int M_1 = *M - 1;
      int kappa_length = *M + 4;
      int i_1000;
      int dim_CPO = dim_Y;
      if (*cop > 0) {
      dim_CPO = *G * *N;
      }
      
      double mu_mean, mu_var, mu_sd;
      double ll_total, can_ll_total, rhoquadsum, canrhoquadsum, dummy_lp, dummy_can_lp, canrho, canrhoquad, canlog_det, can_lp, dummy_log_det, rho_lp,
      alpha_ll, can_alpha_ll, lambda_ll ,can_lambda_ll, delta_ll, can_delta_ll, lambda_lp, can_lambda_lp, delta_lp, can_delta_lp, qf1, qf2, t1, t2, can_lp_0, lp_0;
      ll_total = 0.0;
      can_ll_total = 0.0;
      dummy_log_det = 0.0;
      
      int this_thetastar, thetastar_index;
      int m_init = 0;
      int update, this_theta, prior_index_0;
      int cop_index = 0;
      int prior_index = 0 + cop_index; /*compiler wanted this*/
      int k0 = 0;
      int d, g, h, i, i2, j, k, k1, l, m, n, p, q_int, u, v;
      int dummy_bin = 0;
      int zero = 0;
      v = 0; /*why does the compiler want this?*/
      
      int burn = MCMC_ints[0];
      int keep = MCMC_ints[1];
      int verbose = MCMC_ints[2];
      int thin = MCMC_ints[3];
      int sweeps = burn + keep;
      int num_keep = keep / thin;
      
      int missing_sum = 0;
      
      double canshape, logit_shape, can_logit_shape;
      double scale = 1.0;
      double mn = 0.0;
      double LPML = 0.0;
      double DIC = 0.0;
      double WAIC = 0.0;
      double epostlogpy = 0.0;
      double this_dens = 0.0;
      double lppd_sum = 0.0;
      double pwaic2 = 0.0;
      double elppd = 0.0;
      double llpm = 0.0;
      double pdic = 0.0;
      
      
      int Dummy_one = 1;
      int tail_index;
      
      double pareto_prob_low = pareto_probs[0];
      double pareto_prob_high = pareto_probs[1];
      int pareto_low = 0;
      int can_pareto = 0;
      double random_uniform = runif (0,1);
      double prob_pareto_stay;
      if (random_uniform < pareto_prob_low) {
      pareto_low = 1;
      }
      int update_pareto_low = 1;
      if (pareto_prob_low == 0.0) {
      update_pareto_low = 0;
      }
      if (pareto_prob_low == 1.0) {
      update_pareto_low = 0;
      }
      
      int pareto_high = 0;
      random_uniform = runif (0,1);
      if (random_uniform < pareto_prob_high) {
      pareto_high = 1;
      }
      int update_pareto_high = 1;
      if (pareto_prob_high == 0.0) {
      update_pareto_high = 0;
      }
      if (pareto_prob_high == 1.0) {
      update_pareto_high = 0;
      }
      
      double xi_low_mu_0, xi_low_log_det_g, xi_low_log_det_h, xi_low_rhoquad, xi_low_lp, xi_low_sig_a, xi_low_sig_b, xi_low_prec_mu_0, can_xi_low_lp, can_xi_low_rhoquad;
      double xi_high_mu_0, xi_high_log_det_g, xi_high_log_det_h, xi_high_rhoquad, xi_high_lp, xi_high_sig_a, xi_high_sig_b, xi_high_prec_mu_0, can_xi_high_lp, can_xi_high_rhoquad;
      
      double xi_low_rho = 0.5;
      double xi_high_rho = 0.5;
      
      double xi_low_nu    = xi_nu[0];
      double xi_high_nu   = xi_nu[1];
      
      int basis_correlation = *correlation;
      int xi_low_correlation = *correlation;
      int xi_high_correlation = *correlation;
      /*** data and regression variables ***/
      int *n_g                    = (int *)Calloc(*G, int);
      
      int *status_2            = (int *)Calloc(dim_Y, int);
      
      double *ll_sum              = (double *)Calloc(dim_ll, double);
      double *can_ll_sum          = (double *)Calloc(dim_ll, double);
      double *log_like            = (double *)Calloc(dim_Y, double);
      double *can_log_like        = (double *)Calloc(dim_Y, double);
      double *CPO                 = (double *)Calloc(dim_CPO, double);
      double *lppd_vec            = (double *)Calloc(dim_CPO, double);
      double *log_dens_mean       = (double *)Calloc(dim_CPO, double);  /**** mean of log densities ****/
      double *log_dens_var          = (double *)Calloc(dim_CPO, double);  /**** variance of log densities ****/
      
      double *y                   = (double *)Calloc(dim_Y, double);
      double *y_low               = (double *)Calloc(dim_Y, double);
      double *y_high              = (double *)Calloc(dim_Y, double);
      double *Y_mean              = (double *)Calloc(dim_Y, double);
      double *Y_var               = (double *)Calloc(dim_Y, double);
      double *U_mean            = (double *)Calloc(dim_Y, double);
      double *U_var             = (double *)Calloc(dim_Y, double);
      
      double *U                 = (double *)Calloc(dim_Y, double); /* quantiles */
      double *can_U             = (double *)Calloc(dim_Y, double);
      double *U_low             = (double *)Calloc(dim_Y, double);
      double *U_high            = (double *)Calloc(dim_Y, double);
      double *qnorm_U           = (double *)Calloc(dim_Y, double); /*Z values*/
      int *bin                    = (int *)Calloc(dim_Y, int);
      int *bin_low                = (int *)Calloc(dim_Y, int);
      int *bin_high               = (int *)Calloc(dim_Y, int);
      
      double *cop_ll_sum             = (double *)Calloc(*G ,    double);
      double *can_cop_ll_sum         = (double *)Calloc(*G ,    double);
      
      double *DD                     = (double *)Calloc(dim_Y,    double);
      double *DD_inv                 = (double *)Calloc(dim_Y,    double);
      double *can_DD                 = (double *)Calloc(dim_Y,    double);
      double *can_DD_inv             = (double *)Calloc(dim_Y,    double);
      
      int *missing_index             = (int *)Calloc(*G * *N, int);
      int *dead_index                = (int *)Calloc(*G * *N, int); /*used to update alpha*/
      int *dead_sum                  = (int *)Calloc(*G     , int); /*used to update lambda*/
      int *start_index               = (int *)Calloc(*G * *N, int); /*used to update alpha*/
      
      double *theta           = (double *)Calloc(dim_theta, double);
      double *cantheta        = (double *)Calloc(dim_theta, double);
      double *canthetastar    = (double *)Calloc(dim_thetastar, double);
      double *tuning_theta    = (double *)Calloc(dim_thetastar, double);
      int *acc_thetastar      = (int *)Calloc(dim_thetastar, int);
      int *att_thetastar      = (int *)Calloc(dim_thetastar, int);
      double *meantheta       = (double *)Calloc(dim_theta, double);
      double *X_D             = (double *)Calloc((*J * *D), double);
      double *dummy_seq       = (double *)Calloc((*J), double);
      
      double *mn_H             = (double *)Calloc(     *H, double);
      double *sd_H             = (double *)Calloc(     *H, double);
      
      double *thetastar_keep   = (double *)Calloc(     dim_thetastar, double);
      double *theta_keep       = (double *)Calloc(     dim_theta, double);
      
      double *meanthetastar = (double *)Calloc(GDHP * 3, double);
      double *varthetastar  = (double *)Calloc(GDHP * 2, double);
      double *corrthetastar = (double *)Calloc(GDHP, double);
      double *z             = (double *)Calloc(2, double);
      
      /*** basis matrix variables ***/
      int length_B = *M * (*M - 1);
      if (basis_ind[0] == -1) {
      length_B = kappa_length * *M;
      }
      
      double *B  = (double *)Calloc(length_B, double); /*basis matrix*/
      double *can_B  = (double *)Calloc(length_B, double); /*candidate basis matrix*/
      double *B_q  = (double *)Calloc(*N_tau_q * *M, double); /*basis matrix*/
      
      /*spline bases*/
      double *AAAA = (double *)Calloc((4 * *M  * M_1), double);
      
      /*** prior variables ***/
      double *mu              = (double *)Calloc(     DMP * *H, double);
      double *sig_h           = (double *)Calloc(     DMP * H2, double);
      double *rho             = (double *)Calloc(     DMP, double);
      double *sig_a           = (double *)Calloc(     DMP, double);
      double *sig_b           = (double *)Calloc(     DMP, double);
      
      double *omega = (double *)Calloc(dim_omega, double);
      double *canomega = (double *)Calloc(dim_omega, double);
      
      double *prec_g         = (double *)Calloc(dim_prec_g, double);
      double *can_prec_g     = (double *)Calloc(dim_prec_g, double);
      double *dummy_prec_g   = (double *)Calloc(G2, double);
      
      double *prec_h         = (double *)Calloc(dim_prec_h, double);
      double *dummy_sig_h    = (double *)Calloc(H2, double);
      double *dummy_prec_h   = (double *)Calloc(H2, double);
      
      double *lp = (double *)Calloc(DMP, double);
      double *rhoquad  = (double *)Calloc(DMP, double);
      double *canrhoquad_array  = (double *)Calloc(DMP, double);
      double *sumalpha  = (double *)Calloc(DMP, double);
      
      double *log_det_g  = (double *)Calloc(DMP, double);
      double *log_det_h  = (double *)Calloc(DMP, double);
      
      double *mu_array                  = (double *)Calloc(dim_thetastar, double);
      double *dummy_thetastar_gh        = (double *)Calloc(dim_thetastar, double);
      
      double *g1_vec                    = (double *)Calloc(*G, double);
      double *g2_vec                    = (double *)Calloc(*G, double);
      double *col_sums                  = (double *)Calloc(*G, double);
      double *dummy_mean_matrix         = (double *)Calloc(*G * *G * *H, double);
      double *dummy_mean_mu             = (double *)Calloc(GH, double);
      double *mean_mu                   = (double *)Calloc(*H, double);
      double *chol_mu                   = (double *)Calloc(H2, double);
      double *dummy_chol_mu             = (double *)Calloc(H2, double);
      double *z_mu                      = (double *)Calloc(*H, double);
      double *sai_h                     = (double *)Calloc(H2, double);
      double *W_h                       = (double *)Calloc(H2, double);
      
      /*** shape variables ***/
      double *shape = (double *)Calloc(*H, double);
      double *shape_mean = (double *)Calloc(*H, double);
      double *shape_var = (double *)Calloc(*H, double);
      int *acc_shape = (int *)Calloc(*H, int);
      double *tuning_shape = (double *)Calloc(*H, double);
      double *shapeprec = (double *)Calloc(*H, double);
      double *meanshape = (double *)Calloc(*H, double);
      
      /*** copula variables ***/
      double *W                       = (double *)Calloc(dim_Y, double);
      double *can_W                   = (double *)Calloc(dim_Y, double);
      double *W2                      = (double *)Calloc(dim_Y, double);
      double *can_W2                  = (double *)Calloc(dim_Y, double);
      
      double *copula_loglike          = (double *)Calloc(*G * *N, double);
      double *can_copula_loglike      = (double *)Calloc(*G * *N, double);
      
      /* gamma variables*/
      double *gamma                        = (double *)Calloc(dim_gamma, double);
      double *omega_gamma                  = (double *)Calloc(HQ,           double);
      double *Omega_gamma                  = (double *)Calloc(HQ * HQ,      double);
      double *Omega_gamma_inv              = (double *)Calloc(HQ * HQ,      double);
      double *Omega_gamma_inv_chol         = (double *)Calloc(HQ * HQ,      double);
      
      double *dummy_vec_HQ           = (double *)Calloc(HQ,           double);
      double *z_HQ                   = (double *)Calloc(HQ,           double);
      double *Z2_sum                 = (double *)Calloc(dim_Z2_sum,   double);
      
      double *meangamma                  = (double *)Calloc(dim_gamma, double);
      
      /*eta variables*/
      double *eta                    = (double *)Calloc(dim_eta, double);
      double *omega_eta              = (double *)Calloc(JH,      double);
      double *Omega_eta_inv          = (double *)Calloc(*G * JH * JH, double);
      double *Omega_eta_inv_chol     = (double *)Calloc(*G * JH * JH, double);
      
      double *dummy_vec_JH           = (double *)Calloc(JH,           double);
      double *z_JH                   = (double *)Calloc(JH,           double);
      
      double *eta_cov                = (double *)Calloc(dim_eta_cov,           double);
      double *eta_prec               = (double *)Calloc(dim_eta_cov,           double);
      double *W_cov                  = (double *)Calloc(dim_eta_cov,           double);
      
      double *can_eta_cov            = (double *)Calloc(dim_eta_cov,           double);
      double *can_eta_prec           = (double *)Calloc(dim_eta_cov,           double);
      double *can_W_cov              = (double *)Calloc(dim_eta_cov,           double);
      
      double *cov_JH                 = (double *)Calloc(J2 * H2,           double);
      double *prec_JH                = (double *)Calloc(J2 * H2,           double);
      
      double *meaneta                  = (double *)Calloc(dim_eta, double);
      
      /*alpha variables*/
      double *alpha                  = (double *)Calloc(*G,           double);
      double *can_alpha              = (double *)Calloc(*G,           double);
      
      double *AR                     = (double *)Calloc(dim_this_AR,  double);
      double *AR_prec                = (double *)Calloc(dim_this_AR,  double);
      double *AR_chol                = (double *)Calloc(dim_this_AR,  double);
      double *alpha_log_det           = (double *)Calloc(*G,           double);
      
      double *can_AR                 = (double *)Calloc(dim_this_AR,  double);
      double *can_AR_prec            = (double *)Calloc(dim_this_AR,  double);
      double *can_AR_chol            = (double *)Calloc(dim_this_AR,  double);
      double *can_alpha_log_det       = (double *)Calloc(*G,           double);
      
      double *cov_J                  = (double *)Calloc(J2,           double);
      double *prec_J                 = (double *)Calloc(J2,           double);
      
      double *meanalpha              = (double *)Calloc(*G,           double);
      double *dummy_mat_J2           = (double *)Calloc(*G * J2,      double);
      double *alpha_det_matrix       = (double *)Calloc(*G * J2,      double);
      
      /*lambda variables*/
      double *lambda                 = (double *)Calloc(dim_lambda,   double);
      double *lambda_prec            = (double *)Calloc(dim_lambda,   double);
      double *lambda_chol            = (double *)Calloc(dim_lambda,   double);
      double *lambda_log_det          = (double *)Calloc(*G,           double);
      
      double *can_lambda             = (double *)Calloc(dim_lambda,   double);
      double *can_lambda_prec        = (double *)Calloc(dim_lambda,   double);
      double *can_lambda_chol        = (double *)Calloc(dim_lambda,   double);
      double *can_lambda_log_det      = (double *)Calloc(*G,           double);
      double *tuning_nu              = (double *)Calloc(*G,           double);
      int *acc_lambda                = (int *)Calloc(*G,           int);
      
      double *cov_H                  = (double *)Calloc(H2,           double);
      double *prec_H                 = (double *)Calloc(H2,           double);
      
      double *meanlambda                 = (double *)Calloc(dim_lambda,   double);
      
      /*delta variables*/
      double *delta                  = (double *)Calloc(dim_delta,           double);
      double *delta_prec             = (double *)Calloc(dim_delta,           double);
      double *can_delta              = (double *)Calloc(dim_delta,           double);
      double *can_delta_prec         = (double *)Calloc(dim_delta,           double);
      
      double *tuning_delta           = (double *)Calloc(dim_delta,           double);
      int *acc_delta                 = (int *)Calloc(dim_delta,           int);
      
      double *zeta_A                 = (double *)Calloc(dim_delta,           double);
      double *zeta_B                 = (double *)Calloc(dim_delta,           double);
      
      double *meandelta              = (double *)Calloc(dim_delta,           double);
      
      /*** tail variables ***/
      double *xi_low              = (double *)Calloc(dim_xi, double);
      double *xi_high             = (double *)Calloc(dim_xi, double);
      
      double *can_xi_low          = (double *)Calloc(dim_xi, double);
      double *can_xi_high         = (double *)Calloc(dim_xi, double);
      
      double *xi_low_mu           = (double *)Calloc(*H, double);
      double *xi_high_mu          = (double *)Calloc(*H, double);
      double *xi_low_mu_array     = (double *)Calloc(dim_xi, double);
      double *xi_high_mu_array    = (double *)Calloc(dim_xi, double);
      
      double *xi_low_sai_0        = (double *)Calloc(H2, double);
      double *xi_high_sai_0       = (double *)Calloc(H2, double);
      double *xi_low_prec_g       = (double *)Calloc(G2, double);
      double *xi_low_sig_h        = (double *)Calloc(H2, double);
      double *xi_low_prec_h       = (double *)Calloc(H2, double);
      double *xi_low_omega        = (double *)Calloc(G2 * H2, double);
      double *can_xi_low_omega    = (double *)Calloc(G2 * H2, double);
      double *xi_high_prec_g      = (double *)Calloc(G2, double);
      double *xi_high_sig_h       = (double *)Calloc(H2, double);
      double *xi_high_prec_h      = (double *)Calloc(H2, double);
      double *xi_high_omega       = (double *)Calloc(G2 * H2, double);
      double *can_xi_high_omega   = (double *)Calloc(G2 * H2, double);
      
      double *meanxilow           = (double *)Calloc(dim_xi, double);
      double *tuning_xi_low       = (double *)Calloc(dim_xi, double);
      int *acc_xi_low             = (int *)Calloc(dim_xi, int);
      double *meanxihigh          = (double *)Calloc(dim_xi, double);
      double *tuning_xi_high      = (double *)Calloc(dim_xi, double);
      int *acc_xi_high            = (int *)Calloc(dim_xi, int);
      
      double *I_low              = (double *)Calloc(*M, double);
      double *M_low              = (double *)Calloc(*M, double);
      double *I_high             = (double *)Calloc(*M, double);
      double *M_high             = (double *)Calloc(*M, double);
      double *w                  = (double *)Calloc(*M, double);
      
      /*** data variables ***/
      for (g = 0; g < *G; g++) {
      n_g[g] = 0;
      for (j = 0; j < *J; j++) {
      for (h = 0; h < *H; h++) {
      if (n_g[g]  < n_gjh[g + *G * (j + *J * h)]) {n_g[g] = n_gjh[g + *G * (j + *J * h)];}
      }
      }
      }
      
      k = 0;
      if (basis_ind[0] == -1) {
      k = 3;
      }
      for (i = 0; i < dim_Y; i++) {
      y[i]            = Y_array[0 + 7 * i];
      y_low[i]        = Y_array[1 + 7 * i];
      y_high[i]       = Y_array[2 + 7 * i];
      Y_mean[i]       = 0.0;
      Y_var[i]        = 0.0;
      U_mean[i]     = 0.0;
      U_var[i]      = 0.0;
      bin[i]          = k;
      bin_low[i]      = k;
      bin_high[i]     = k;
      log_like[i]     = 0.0;
      can_log_like[i] = 0.0;
      U[i]          = 0.5;
      can_U[i]      = 0.5;
      U_low[i]      = 0.0;
      U_high[i]     = 1.0;
      qnorm_U[i]    = 0.0;
      status_2[i]     = 0;
      }
      for (i = 0; i < dim_CPO; i++) {
      CPO[i] = 0.0;
      lppd_vec[i] = 0.0;
      log_dens_mean[i] = 0.0;
      log_dens_var[i] = 0.0;
      }
      /*** regression variables ***/
      for (i = 0; i < dim_thetastar; i++) {
      canthetastar[i] = thetastar[i];
      tuning_theta[i] = 0.01;
      acc_thetastar[i] = 0;
      att_thetastar[i] = 0;
      }
      
      /*initialize theta = X_D %*% thetastar*/
      for (j = 0; j < *J; j++) {
      dummy_seq[j] = (j + 1) /  (double) (*J);
      }
      for (k = 0; k < dim_theta; k++) {
      theta[k] = 0.0;
      }
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      for (p = 0; p < *P; p++) {
      for (m = 0; m < *M; m++) {
      this_theta = GJH * (m + *M * p);
      if (D_vec[p] == 1) {
      for (j = 0; j < *J; j++) {
      theta[g + *G * j + *G * *J * h + this_theta] = thetastar[g + *G * 0 + *G * *D * (h + *H * m + *H * *M * p)];
      }
      } else{
      for (j = 0; j < *J; j++) {
      for (d = 0; d < D_vec[p]; d++) {
      X_D[j + *J * d] = pow(dummy_seq[j],d);
      }
      }
      for (j = 0; j < *J; j++) {
      k = g + *G * j + *G * *J * h + this_theta;
      theta[k] = 0;
      for (d = 0; d < D_vec[p]; d++) {
      theta[k] += X_D[j + *J * d] * thetastar[g + *G * d + *G * *D * h + *G * *D * *H * (m + *M * p)];
      }
      }
      }
      }
      }
      }
      }
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      for (p = 0; p < *P; p++) {
      for (m = 0; m < *M; m++) {
      update = 1;
      if (m == 0 && basis_ind[h] > 3) {
      update = 0;
      }/*for gamma and weibull only threshold for m > 0*/
      if (update == 1) {
      threshold(G, J, D_vec, D, H, M, P, &g, &h, &m, &p, timepoints, thetastar, theta);
      }
      }
      }
      }
      }
      for (i = 0; i < dim_theta; i++) {
      cantheta[i] = theta[i];
      meantheta[i] = 0;
      }
      /*** basis matrix variables ***/
      /*initialize B*/
      for (i = 0; i < length_B; i++) {
      B[i] = 0.0;
      }
      /*construct pointers to the quantile function, log cdf and log density*/
      void (*q_ptr[*H])(double *, double *, double *, double *, double*); /*array of pointers to the quantile function*/
      void (*p_ptr[*H])(double *, double *, double *, double *, double*); /*array of pointers to the CDF*/
      void (*log_d_ptr[*H])(double *, double *, double *, double *, double*); /*array of pointers to the log density*/
      
      for (h = 0; h < *H; h++) {
      if (basis_ind[h] <= 0) {
      q_ptr[h]           = &QNorm;
      p_ptr[h]           = &PNorm;
      log_d_ptr[h]       = &LogDNorm;
      } else if (basis_ind[h] == 1) {
      q_ptr[h]           = &QT;
      p_ptr[h]           = &PT;
      log_d_ptr[h]       = &LogDT;
      } else if (basis_ind[h] == 2) {
      q_ptr[h]           = &QLogistic;
      p_ptr[h]           = &PLogistic;
      log_d_ptr[h]       = &LogDLogistic;
      } else if (basis_ind[h] == 3) {
      q_ptr[h]           = &QAsymmetricLaplace;
      p_ptr[h]           = &PAsymmetricLaplace;
      log_d_ptr[h]       = &LogDAsymmetricLaplace;
      } else if (basis_ind[h] == 4) {
      q_ptr[h]           = &QWeibull;
      p_ptr[h]           = &PWeibull;
      log_d_ptr[h]       = &LogDWeibull;
      } else{
      q_ptr[h]           = &QGamma;
      p_ptr[h]           = &PGamma;
      log_d_ptr[h]       = &LogDGamma;
      }
      }
      /*construct pointers to the likelihoods and quantile functions*/
      void (*clike_ptr)(int *, int *, int *, int *, int *,
      int *, int *, int *, double *,
      double *, double *, double *, double *,
      double *,
      double *,
      void (double *, double *, double *, double *, double*),
      void (double *, double *, double *, double *, double*),
      void (double *, double *, double *, double *, double*),
      double *, double *,
      int *, int *, int *, double *, double *, double *, double *, double *,
      int *, int *, double *, double *, double *, double *,
      double *, double *, double *, double *, double *,
      int *, int *, int *);
      
      void (*qf_ptr)(int *, int *, int *, int *, int *,
      int *, int *, int *, double *,
      double *, double *, double *, double *,
      double *,
      double *,
      void (double *, double *, double *, double *, double*),
      void (double *, double *, double *, double *, double*),
      void (double *, double *, double *, double *, double*),
      double *,
      double *,
      int *, int *, int *, double *, double *, double *, double *, double *,
      int *, int *, double *, double *, double *, double *,
      int *, int *, int *);
      
      if (basis_ind[0] == -1) {
      clike_ptr = &clike_spline;
      qf_ptr = &quantile_function_spline;
      }
      if (basis_ind[0] >=  0) {
      clike_ptr = &clike;
      qf_ptr = &quantile_function;
      }
      /*spline basis*/
      if (basis_ind[0] == -1) {
      for (m = 0; m < 4 * *M * M_1; m++) {
      AAAA[m] = 0.0;
      }
      MakeCubicSplineCoefficients(M, kappa, AAAA);
      }
      /*** prior variables ***/
      for (i = 0; i <      DMP * H2; i++) {
      sig_h[i] = 0.0;
      }
      for (d = 0; d < *D; d++) {
      for (m = 0; m < *M; m++) {
      for (p = 0; p < *P; p++) {
      for (h = 0; h < *H; h++) {
      mu[d + *D * m + *D * *M * p + DMP * h] = mu_0[d + *D * m + *D * *M * p];
      for (v = 0; v < *H; v++) {
      k = d + *D * m + *D * *M * p + DMP * (h + *H * v);
      sig_h[k] = sai_h_0[k];
      }
      }
      }
      }
      }
      
      for (i = 0; i < DMP; i++) {
      rho[i] = 0.5;
      }
      for (i = 0; i < dim_omega; i++) {
      omega[i] = 0.0; canomega[i] = 0.0;
      }
      for (i = 0; i < G2; i++) {
      dummy_prec_g[i] = 0.0;
      }
      for (i = 0; i < *G; i++) {
      dummy_prec_g[i + *G * i] = 1.0;
      }
      for (i = 0; i < H2; i++) {
      dummy_sig_h[i] = 0.0; dummy_prec_h[i] = 0.0;
      }
      for (i = 0; i < *H; i++) {
      dummy_prec_h[i + *H * i] = 1.0;
      }
      for (i = 0; i < dim_prec_g; i++) {
      prec_g[i] = 0.0;
      }
      for (i = 0; i < dim_prec_h; i++) {
      prec_h[i] = 0.0;
      }
      /*recenter mu_0 for mcmc - no need to update every time*/
      for (d = 0; d < *D; d++) {
      for (m = 0; m < *M; m++) {
      for (p = 0; p < *P; p++) {
      prior_index = d + *D * m + *D * *M * p;
      mu_0[prior_index] = mu_0[prior_index] * prec_mu_0[prior_index];
      }
      }
      }
      
      void (*make_prec_ptr)(int *, double *, double *, double *, double *); /*pointer to the correlation type across models*/
      if (*G == 1) {
      basis_correlation = 0; xi_low_correlation = 0; xi_high_correlation = 0;
      }
      if (basis_correlation == 0) {
      make_prec_ptr = &MakeDiagonalPrecision;
      }
      if (basis_correlation == 1) {
      make_prec_ptr = &MakeAutoregressivePrecision;
      }
      if (basis_correlation == 2) {
      make_prec_ptr = &MakeSpatialPrecision;
      }
      /*initialize prior arrays*/
      for (d = 0; d < *D; d++) {
      for (m = 0; m < *M; m++) {
      for (p = 0; p < *P; p++) {
      prior_index = d + *D * m + *D * *M * p;
      sig_a[prior_index] = nu_h[prior_index];
      sig_b[prior_index] = sai_h_0[prior_index];
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      k = d + *D * h + *D * *H * (m + *M * p);
      mu_array[g + *G * k] = mu[d + *D * m + *D * *M * p + DMP * h];
      }
      }
      make_prec_ptr(G, &rho[prior_index], dist, dummy_prec_g, &log_det_g[prior_index]);
      for (i = 0; i < *G; i++) {
      for (j = 0; j < *G; j++) {
      prec_g[prior_index + DMP * (i + *G * j)] = dummy_prec_g[i + *G * j];
      }
      }
      for (i = 0; i < *H; i++) {
      for (j = 0; j < *H; j++) {
      dummy_sig_h[i + *H * j] = sig_h[prior_index + DMP * (i + *H * j)];
      }
      }
      InvertSymmetricMatrix(H, dummy_sig_h, dummy_prec_h, &log_det_h[prior_index]);
      for (i = 0; i < *H; i++) {
      for (j = 0; j < *H; j++) {
      prec_h[prior_index  + DMP * (i + *H * j)] = dummy_prec_h[i + *H * j];
      }
      }
      make_omega(G, D, H, M, P, &d, &m, &p, prec_g,prec_h, omega);
      make_omega(G, D, H, M, P, &d, &m, &p, prec_g,prec_h, canomega);
      rho_quad  (G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
      lp[prior_index] = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - rhoquad[prior_index]);
      if (verbose) {
      Rprintf("initial lp = %f\n", lp[prior_index]);
      }
      }
      }
      }
      for (k = 0; k < DMP; k++) {
      canrhoquad_array[k] = 0.0;
      }
      for (i = 0; i < dim_prec_g; i++) {
      can_prec_g[i] = prec_g[i];
      }
      /*** shape variables ***/
      for (h = 0; h < *H; h++) {
      shape_mean[h] = shape_prior[0 + 4 * h];
      shape_var[h]  = shape_prior[1 + 4 * h];
      shape[h] = shape_mean[h];
      acc_shape[h] = 0; shapeprec[h] = 1 / (2 * shape_var[h]); tuning_shape[h] = 1.0; meanshape[h] = 0.0;
      mn_H[h] = shape_prior[2 + 4 * h];
      sd_H[h] = shape_prior[3 + 4 * h];
      }
      /******* tail variables ****/
      for (k = 0; k < dim_xi; k++) {
      xi_low[k] = runif (0,1) * (-1.0 - -3.0) + -3.0;
      xi_high[k] = runif (0,1) * (-1.0 - -3.0) + -3.0;
      can_xi_low[k] = xi_low[k];
      can_xi_high[k] = xi_high[k];
      acc_xi_low[k] = 0;
      tuning_xi_low[k] = 0.1;
      meanxilow[k] = 0.0;
      acc_xi_high[k] = 0;
      tuning_xi_high[k] = 0.1;
      meanxihigh[k] = 0.0;
      }
      
      for (k = 0; k < G2; k++) {
      xi_low_prec_g[k] = 0.0;
      }
      for (k = 0; k < *G; k++) {
      xi_low_prec_g[k + *G * k] = 1.0;
      }
      for (k = 0; k < H2; k++) {
      xi_low_prec_h[k] = 0.0; xi_low_sig_h[k] = 0.0;
      }
      
      for (k = 0; k < G2; k++) {
      xi_high_prec_g[k] = 0.0;
      }
      for (k = 0; k < *G; k++) {
      xi_high_prec_g[k + *G * k] = 1.0;
      }
      for (k = 0; k < H2; k++) {
      xi_high_prec_h[k] = 0.0;
      xi_high_sig_h[k] = 0.0;
      }
      
      for (h = 0; h < *H; h++) {
      xi_low_mu[h] = xi_mu_0[0];
      xi_high_mu[h] = xi_mu_0[1];
      }
      xi_low_mu_0 = xi_mu_0[0];
      xi_high_mu_0 = xi_mu_0[1];
      xi_low_prec_mu_0                = xi_prec_mu_0[0];
      xi_high_prec_mu_0               = xi_prec_mu_0[1];
      
      for (u = 0; u < *H; u++) {
      for (v = 0; v < *H; v++) {
      xi_low_sai_0[u + *H * v]        = xi_sai_0[u + *H * v + *H * *H * 0];
      xi_high_sai_0[u + *H * v]       = xi_sai_0[u + *H * v + *H * *H * 1];
      }
      }
      /*recenter mu_0 for mcmc - no need to update every time*/
      xi_low_mu_0 *= xi_low_prec_mu_0;
      xi_high_mu_0 *= xi_high_prec_mu_0;
      for (h = 0; h < H2; h++) {
      xi_low_sig_h[h] = xi_low_sai_0[h];
      xi_high_sig_h[h] = xi_high_sai_0[h];
      }
      xi_low_sig_a = xi_low_nu;
      xi_low_sig_b = xi_low_sai_0[0];
      xi_high_sig_a = xi_high_nu;
      xi_high_sig_b = xi_high_sai_0[0];
      
      I_low[0] = 1.0;
      M_low[0] = 0.0;
      I_high[0] = 1.0;
      M_high[0] = 0.0;
      w[0] = 0.0;
      for (m = 0; m < M_1; m++) {
      w[m + 1] = 0.0;
      CubicMSpline(tau_low, &kappa_length, kappa, &m, &dummy_bin, &M_low[m + 1]);
      M_low[m + 1] *= *tau_low;
      CubicISpline(tau_low, &kappa_length, kappa, &m, &dummy_bin, &I_low[m + 1]);
      CubicMSpline(tau_high, &kappa_length, kappa, &m, &dummy_bin, &M_high[m + 1]);
      M_high[m + 1] *= (1 - *tau_high);
      CubicISpline(tau_high, &kappa_length, kappa, &m, &dummy_bin, &I_high[m + 1]);
      }
      
      for (i = 0; i < DMP * H2; i++) {
      sig_h[i] = 0.0;
      }
      for (d = 0; d < *D; d++) {
      for (m = 0; m < *M; m++) {
      for (p = 0; p < *P; p++) {
      for (h = 0; h < *H; h++) {
      for (v = 0; v < *H; v++) {
      k = d + *D * m + *D * *M * p + DMP * (h + *H * v);
      sig_h[k] = sai_h_0[k];
      }
      }
      }
      }
      }
      
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      xi_low_mu_array[g + *G * h] = xi_low_mu[h];
      xi_high_mu_array[g + *G * h] = xi_high_mu[h];
      }
      }
      make_prec_ptr     (G, &xi_low_rho, dist, xi_low_prec_g, &xi_low_log_det_g);
      InvertSymmetricMatrix(H, xi_low_sig_h,     xi_low_prec_h, &xi_low_log_det_h);
      make_prec_ptr     (G, &xi_high_rho, dist, xi_high_prec_g, &xi_high_log_det_g);
      InvertSymmetricMatrix(H, xi_high_sig_h,     xi_high_prec_h, &xi_high_log_det_h);
      
      d = 0;
      m = 0;
      p = 0;
      
      make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_prec_g, xi_low_prec_h, xi_low_omega);
      rho_quad  (G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
      xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - xi_low_rhoquad);
      make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_prec_g, xi_high_prec_h, xi_high_omega);
      rho_quad  (G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
      xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - xi_high_rhoquad);
      
      /*** copula variables ***/
      for (i = 0; i < *G * *N; i++) {
      copula_loglike[i] = 0.0;
      can_copula_loglike[i] = 0.0;
      }
      /*gamma variables*/
      for (h = 0; h < dim_gamma; h++) {
      gamma[h] = 0.0;
      meangamma[h] = 0.0;
      }
      for (h = 0; h < HQ * HQ; h++) {
      Omega_gamma[h] = 0.0;
      Omega_gamma_inv[h] = 0.0;
      Omega_gamma_inv_chol[h] = 0.0;
      }
      /*eta variables*/
      for (j = 0; j < dim_eta; j++) {
      eta[j] = 0.0;
      meaneta[j] = 0.0;
      }
      for (j = 0; j < *G * JH * JH; j++) {
      Omega_eta_inv[j] = 0.0;
      Omega_eta_inv_chol[j] = 0.0;
      }
      for (j = 0; j < dim_eta_cov; j++) {
      eta_cov[j]  = 0.0;
      eta_prec[j] = 0.0;
      W_cov[j]    = 0.0;
      }
      for (j = 0; j < JH; j++) {
      dummy_vec_JH[j] = 0.0;
      }
      
      if (*cop > 0) {
      for (g = 0; g < *G; g++) {
      for (n = 0; n < n_g[g]; n++) {
      for (u = 0; u < HQ; u++) {
      for (v = 0; v < HQ; v++) {
      dummy_log_det = 0;
      for (j = 0; j < JH; j++) {
      dummy_log_det += Z[g + *G * (j + JH * n + JH * *N * u)] * Z[g + *G * (j + JH * n + JH * *N * v)];
      }
      Z2_sum[g + *G * n + *G * *N * (u + HQ * v)] = dummy_log_det;
      }
      }
      }
      }
      }
      /*alpha variables*/
      for (j = 0; j < JH * JH; j++) {
      cov_JH[j] = 0.0;
      prec_JH[j] = 0.0;
      }
      for (j = 0; j < J2; j++) {
      cov_J[j] = 0.0;
      prec_J[j] = 0.0;
      dummy_mat_J2[j] = 0.0;
      }
      
      for (i = 0; i < *G * J2; i++) {
      AR[i] = 0.0;
      AR_prec[i] = 0.0;
      AR_chol[i] = 0.0;
      can_AR[i] = 0.0;
      can_AR_prec[i] = 0.0;
      can_AR_chol[i] = 0.0;
      alpha_det_matrix[i] = 0.0;
      }
      for (g = 0; g < *G; g++) {
      meanalpha[g] = 0.0;
      alpha[g] = 0.0;
      alpha_log_det[g] = 0.0;
      for (j = 0; j < *J; j++) {
      k = g + *G * (j + *J * j);
      AR[k] = 1.0;
      AR_chol[k] = 1.0;
      AR_prec[k] = 1.0;
      }
      if (*cop > 1) {
      alpha[g] = 0.5;
      for (u = 0; u < *J; u++) {
      for (v = 0; v < *J; v++) {
      k = u + *J * v;
      cov_J[k] = pow(alpha[g], dist_J[u + *J * v]);
      }
      }
      CholMatrix(J, cov_J, prec_J, &dummy_log_det);
      for (u = 0; u < *J; u++) {
      for (v = 0; v < *J; v++) {
      k = u + *J * v;
      AR_chol[g + *G * k] = prec_J[k];
      }
      }
      longitudinal_log_det(J, cov_J, prec_J, dummy_mat_J2);
      for (j = 0; j < J2; j++) {
      alpha_det_matrix[g + *G * j] = dummy_mat_J2[j];
      }
      InvertSymmetricMatrix(J, cov_J, prec_J, &alpha_log_det[g]);
      for (u = 0; u < *J; u++) {
      for (v = 0; v < *J; v++) {
      k = u + *J * v;
      AR[g + *G * k] = cov_J[k];
      AR_prec[g + *G * k] = prec_J[k];
      }
      }
      }
      can_alpha[g]        = alpha[g];
      can_alpha_log_det[g] = alpha_log_det[g];
      for (j = 0; j < J2; j++) {
      k = g + *G * j;
      can_AR[k]      = AR[k];
      can_AR_chol[k] = AR_chol[k];
      can_AR_prec[k] = AR_prec[k];
      }
      }
      /*lambda variables*/
      for (i = 0; i < dim_lambda; i++) {
      lambda[i] = 0.0;
      lambda_prec[i] = 0.0;
      lambda_chol[i] = 0.0;
      meanlambda[i] = 0.0;
      }
      for (g = 0; g < *G; g++) {
      for (k = 0; k < *H; k++) {
      lambda[g + *G * (k + *H * k)] = 1.0;
      }
      }
      for (h = 0; h < H2; h++) {
      sai_h[h] = 0.0;
      }
      for (h = 0; h < H2; h++) {
      cov_H[h] = 0.0;
      prec_H[h] = 0.0;
      W_h[h] = 0.0;
      }
      
      for (g = 0; g < *G; g++) {
      for (u = 0; u < *H; u++) {
      for (v = 0; v < *H; v++) {
      cov_H[u + *H * v] = lambda[g + *G * (u + *H * v)];
      }
      }
      CholMatrix(H, cov_H, prec_H, &dummy_log_det);
      for (u = 0; u < *H; u++) {
      for (v = 0; v < *H; v++) {
      lambda_chol[g + *G * (u + *H * v)] = prec_H[u + *H * v];
      }
      }
      InvertSymmetricMatrix(H, cov_H, prec_H, &lambda_log_det[g]);
      for (u = 0; u < *H; u++) {
      for (v = 0; v < *H; v++) {
      lambda_prec[g + *G * (u + *H * v)] = prec_H[u + *H * v];
      }
      }
      for (u = 0; u < *J; u++) {
      for (v = 0; v < *J; v++) {
      k = u + *J * v;
      cov_J[k] = AR[g + *G * k];
      prec_J[k] = AR_prec[g + *G * k];
      }
      }
      Kronecker(J, J, H, H, cov_J, cov_H, cov_JH);
      Kronecker(J, J, H, H, prec_J, prec_H, prec_JH);
      
      for (u = 0; u < JH; u++) {
      for (v = 0; v < JH; v++) {
      k = u + JH * v;
      eta_cov[g + *G * k] = cov_JH[k];
      eta_prec[g + *G * k] = prec_JH[k];
      W_cov[g + *G * k] = eta_cov[g + *G * k];
      can_W_cov[g + *G * k] = W_cov[g + *G * k];
      }
      }
      for (u = 0; u < JH; u++) {
      W_cov[g + *G * (u + JH * u)] ++;
      can_W_cov[g + *G * (u + JH * u)] ++;
      cov_JH[u + JH * u] ++;
      prec_JH[u + JH * u] ++;
      }
      InvertSymmetricMatrix(&JH, prec_JH, cov_JH, &dummy_log_det);
      for (u = 0; u < JH * JH; u++) {
      Omega_eta_inv[g + *G * u] = cov_JH[u];
      }
      CholMatrix(&JH, cov_JH, prec_JH, &dummy_log_det);
      for (u = 0; u < JH * JH; u++) {
      Omega_eta_inv_chol[g + *G * u] = prec_JH[u];
      }
      }
      for (g = 0; g < *G; g++) {
      can_lambda_log_det[g] = lambda_log_det[g];
      tuning_nu[g] = *H + 2.0;
      acc_lambda[g] = 0;
      }
      for (i = 0; i < dim_lambda; i++) {
      can_lambda[i] = lambda[i];
      can_lambda_chol[i] = lambda_chol[i];
      can_lambda_prec[i] = lambda_prec[i];
      }
      
      /*delta variables*/
      for (g = 0; g < *G * HQ; g++) {
      delta[g]            = 1.0;
      delta_prec[g]       = 1 / delta[g];
      can_delta[g]        = delta[g];
      can_delta_prec[g]   = delta_prec[g];
      tuning_delta[g]     = 1.0;
      acc_delta[g]        = 0;
      meandelta[g]        = 0.0;
      }
      
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      for (q_int = 0; q_int < *Q_int; q_int++) {
      k = q_int + *Q_int * h;
      zeta_A[g + *G * k] = zeta[g + *G * (k + HQ * 0)];
      zeta_B[g + *G * k] = zeta[g + *G * (k + HQ * 1)];
      }
      }
      }
      /*regression variables*/
      for (i = 0; i < dim_Y; i++) {
      W[i] = 0.0;
      can_W[i] = 0.0;
      W2[i] = 0.0;
      can_W2[i] = 0.0;
      Y_mean[i] = 0.0;
      Y_var[i] = 0.0;
      U_mean[i] = 0.0;
      U_var[i] = 0.0;
      DD[i] = 1.0;
      DD_inv[i] = 0.0;
      }
      
      for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
      for (n = 0; n < *N; n++) {
      for (j = 0; j < *J; j++) {
      status_2[g + *G * (j * *H + h + JH * n)] = status[g + *G * j + *G * *J * h + GJH * n];
      }
      }
      }
      }
      
      if (*cop > 0) {
      for (g = 0; g < *G; g++) {
      make_D(G, J, H, Q_int, n_g, N, Z, delta, W_cov, &g, status, DD_inv, DD);
      }
      }
      for (i = 0; i < dim_Y; i++) {
      can_DD[i] = DD[i];
      can_DD_inv[i] = DD_inv[i];
      }
      
      for (g = 0; g < *G; g++) {
      dead_sum[g] = 0;
      for (n = 0; n < *N; n++) {
      k = g + *G * n;
      missing_index[k] = 0; /*missing people*/
      dead_index[k] = JH;
      for (j = 0; j < JH; j++) {
      u = g + *G * (j + JH * n);
      if (status_2[u] == 0) {
      missing_index[k] ++;
      missing_sum ++;
      }
      /*in longitudinal code set U = 0.5 instead of W = 0*/
      if (status_2[u] == -99) {
      if (*cop > 0) {
      W[u] = 0.0;
      W2[u] = 0.0;
      /*note that I change Z here - make sure you don't use Z from this function in R*/
        for (q_int = 0; q_int < HQ; q_int++) {
          Z[u + *G * JH * *N * q_int] = 0.0;
        }
      }
  dead_index[k] --;
  dead_sum[g] ++;
      }
      }
dead_index[k] /= *H;
start_index[k] = -1;
j = 0;
while(start_index[k] == -1 && j < (*J * *H)) {
  u = g + *G * (j + JH * n);
  if (status_2[u] != -99) {
    start_index[k] = j / *H;
  }
  j++;
}
      }
dead_sum[g] /= *H;
      }
/*initialize at random spot, whether or not there is a copula*/
  if (missing_sum > 0) {
    for (g = 0; g < *G; g++) {
      for (n = 0; n < n_g[g]; n++) {
        if (missing_index[g + *G * n] > 0) {
          for (j = 0; j < JH; j++) {
            k = g + *G * j + *G * JH * n;
            if (status[k] == 0) {
              U[k] = runif (0,1);
              can_U[k] = U[k];
            }
          }
        }
      }
    }
    for (g = 0; g < *G; g++) {
      for (j = 0; j < *J; j++) {
        for (h = 0; h < *H; h++) {
          tail_index = g + *G * h;
          u = g + *G * (j + *J * h);
          
          qf_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                 kappa, X, y, y_low, y_high, theta, &shape[h],
                 q_ptr[h], p_ptr[h], log_d_ptr[h],
                 AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
                 &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                 &g, &j, &h);
        }
      }
    }
  }

for (k = 0; k < GDHP * 3; k++) {
  meanthetastar[k] = 0.0;
}
for (k = 0; k < GDHP * 2; k++) {
  varthetastar[k] = 0.0;
}
for (k = 0; k < GDHP; k++) {
  corrthetastar[k] = 0.0;
}
z[0] = 0.0;
z[1] = 0.0;

/*   initialize the likelihood */
  for (g = 0; g < *G; g++) {
    cop_ll_sum[g] = 0.0;
    for (h = 0; h < *H; h++) {
      k = g + *G * h;
      tail_index = g + *G * h;
      make_b(q_ptr[h], M, kappa, &shape[h], B);
      for (j = 0; j < *J; j++) {
        u = g + *G * j + *G * *J * h;
        ll_sum[u] = 0;
        clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                  kappa, X, y, y_low, y_high, theta, &shape[h],
                  q_ptr[h], p_ptr[h], log_d_ptr[h],
                  AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
                  &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                  w, M_low, M_high, I_low, I_high,
                  &g, &j, &h);
        can_ll_sum[u] = ll_sum[u];
        if (verbose) {
          Rprintf("initial ll = %f\n", ll_sum[u]);
        }
        for (n = 0; n < *N; n++) {
          k = g + *G * j + *G * *J * h + GJH * n;
          W[k] = qnorm(U[k], 0, 1, 1, 0);
        }
      }
    }
  }

if (*cop > 0) {
  for (g = 0; g < *G; g++) {
    for (h = 0; h < *H; h++) {
      for (n = 0; n < *N; n++) {
        for (j = 0; j < *J; j++) {
          W2[g + *G * (j * *H + h + JH * n)] = W[g + *G * j + *G * *J * h + *G * *J * *H * n];
        }
      }
    }
    copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]);
    for (n = 0; n < *N; n++) {
      can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
    }
    can_cop_ll_sum[g] = cop_ll_sum[g];
    if (verbose) {
      Rprintf("initial cop_ll_sum = %f\n", cop_ll_sum[g]);
    }
  }
}
for (k = 0; k < dim_Y; k++) {
  can_log_like[k] = log_like[k];
  can_U[k] = U[k];
  can_W[k] = W[k];
  can_W2[k] = W2[k];
}
///***************** Part II: MCMC *********************/
  if (verbose) {
    Rprintf("burn, burn, burn... \n");
  }
GetRNGstate();
for (i = 0; i < sweeps; i++) {
  if (verbose == 1 && i == burn) {
    Rprintf("Burn-in Finished. \n");
  }
  if (verbose) {
    if (i % 1000 == 0 && i > burn) {
      Rprintf("Keepers %d %% done.\n", 100 * (i - burn) /  (sweeps - burn) );
    }
  }
  /** update hyper parameters **/
    if (GH > 1) {
      for (d = 0; d < *D; d++) {
        for (m = 0; m < *M; m++) {
          for (p = 0; p < *P; p++) {
            prior_index = d + *D * m + *D * *M * p;
            thetastar_index = m + *M * p;
            this_thetastar = *G * *D * *H * thetastar_index;
            if (*H == 1) {
              h = 0;
              /*** update mean via Gibbs sampling ***/
                mu_var = 0.0;
                for (u = 0; u < *G; u++) {
                  for (v = 0; v < *G; v++) {
                    mu_var += prec_g[prior_index + DMP * (u + *G * v)];
                  }
                }
                mu_var *= prec_h[prior_index];
                mu_var += prec_mu_0[prior_index];
                mu_mean = 0.0;
                for (u = 0; u < *G; u++) {
                  for (v = 0; v < *G; v++) {
                    mu_mean += 1 * thetastar[v + *G * d + 0 + *G * *D * *H * (m + *M * p)] * prec_g[prior_index + DMP * (u + *G * v)];
                  }
                }
                mu_mean *= prec_h[prior_index];
                mu_mean += mu_0[prior_index];
                mu_var = 1 / mu_var;
                mu_mean *= mu_var;
                mu_sd = sqrt(mu_var);
                mu[prior_index] = mu_sd * rnorm(0,1) + mu_mean;
                for (g = 0; g < *G; g++) {
                  mu_array[g + *G * d + *G * *D * 0 + this_thetastar] = mu[prior_index];
                }
                rho_quad(G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
                /*** update prec_h via Gibbs sampling ***/
                  dummy_log_det = rhoquad[prior_index] / prec_h[prior_index]; /*remove precision term from quadratic form*/
                  prec_h[prior_index] = rgamma(*G / 2 + sig_a[prior_index], (2 * sig_b[prior_index])/(sig_b[prior_index] * dummy_log_det + 2));
                sig_h[prior_index] = 1 / prec_h[prior_index];
                make_omega(G, D, H, M, P, &d, &m, &p, prec_g, prec_h, omega);
                rho_quad(G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
                log_det_h[prior_index] = log(prec_h[prior_index]);
                lp[prior_index] = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - rhoquad[prior_index]);
            }
            if (*G == 1) {
              /*** update mean via Gibbs sampling ***/
                mu_mean = thetastar[d + *G * *D * 0 + GDH * thetastar_index];
                for (h = 1; h < *H; h++) {
                  mu_mean += thetastar[d + *G * *D * h + GDH * thetastar_index];
                }
                mu_mean *= omega[prior_index];
                mu_mean += mu_0[prior_index];
                mu_var = *H * omega[prior_index] + prec_mu_0[prior_index];
                mu_var = 1 / mu_var;
                mu_mean *= mu_var;
                mu_sd = sqrt(mu_var);
                mu[prior_index] = mu_sd * rnorm(0,1) + mu_mean;
                for (g = 0; g < *G; g++) {
                  for (h = 0; h < *H; h++) {
                    mu_array[g + *G * d + *G * *D * h + this_thetastar] = mu[prior_index];
                  }
                }
                rho_quad(G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
                /*** update prec_h via Gibbs sampling ***/
                  /*nu_h is shape parameter*/
                  /*sai_h_0 is scale parameter*/
                  dummy_log_det = rhoquad[prior_index] / omega[prior_index]; /*remove precision term from quadratic form*/
                  prec_h[prior_index] = rgamma(*H / 2 + sig_a[prior_index], (2 * sig_b[prior_index])/(sig_b[prior_index] * dummy_log_det + 2));
                sig_h[prior_index] = 1 / prec_h[prior_index];
                omega[prior_index] = prec_h[prior_index];
                rho_quad(G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
                log_det_h[prior_index] = log(omega[prior_index]);
                lp[prior_index] = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - rhoquad[prior_index]);
            }
            if (*G > 1 && *H > 1) {
              /*** update mean via Gibbs sampling ***/
                for (h = 0; h < *H; h++) {
                  z_mu[h] = rnorm(0,1);
                }
              gibbs_mu(G, D, H, M, P, &d, &m, &p, thetastar, prec_g, prec_h, z_mu, mu_0, prec_mu_0, mu, mu_array);
              /*** update sig_h via Gibbs sampling ***/
                /*if *G > 1 then we have replication so we can use inverse wishart*/
                make_sai(G, D, H, M, P, &d, &m, &p, thetastar, prec_g, mu_array, sai_h, sai_h_0);
              RandomInverseWishart(H, &nu_h[prior_index], sai_h, W_h);
              InvertSymmetricMatrix(H, W_h, dummy_prec_h, &log_det_h[prior_index]);
              for (u = 0; u < *H; u++) {
                for (v = 0; v < *H; v++) {
                  k = u + *H * v;
                  sig_h [prior_index + DMP * k] = W_h[k];
                  prec_h[prior_index + DMP * k] = dummy_prec_h[k];
                }
              }
              make_omega(G, D, H, M, P, &d, &m, &p, prec_g, prec_h, omega);
              rho_quad(G, D, H, M, P, &d, &m, &p, omega, thetastar, mu_array, &rhoquad[prior_index]);
              lp[prior_index] = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - rhoquad[prior_index]);
            }
          }
        }
      }
    }
  /*** independent metropolis update the between model correlation ***/
    if (basis_correlation > 0) {
      rhoquadsum = 0.0;
      canrhoquadsum = 0.0;
      canrho = runif (0,1) * (*rho_b - *rho_a) + *rho_a;
      make_prec_ptr(G, &canrho, dist, dummy_prec_g, &canlog_det);
      for (d = 0; d < *D; d++) {
        for (m = 0; m < *M; m++) {
          for (p = 0; p < *P; p++) {
            prior_index = d + *D * m + *D * *M * p;
            for (g = 0; g < *G; g++) {
              for (h = 0; h < *G; h++) {
                can_prec_g[prior_index + DMP * (g + *G * h)] = dummy_prec_g[g + *G * h];
              }
            }
            make_omega(G, D, H, M, P, &d, &m, &p, can_prec_g, prec_h, canomega);
            rho_quad  (G, D, H, M, P, &d, &m, &p, canomega, thetastar, mu_array, &canrhoquad);
            canrhoquadsum += canrhoquad;
            rhoquadsum += rhoquad[prior_index];
            canrhoquad_array[prior_index] = canrhoquad;
          }
        }
      }
      can_xi_low_lp = xi_low_lp;
      can_xi_high_lp = xi_high_lp;
      if (xi_low_correlation > 0) {
        d = 0; m = 0; p = 0;
        make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, dummy_prec_g, xi_low_prec_h, can_xi_low_omega);
        rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, can_xi_low_omega, xi_low, xi_low_mu_array, &can_xi_low_rhoquad);
        can_xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * canlog_det - can_xi_low_rhoquad);
      }
      if (xi_high_correlation > 0) {
        d = 0; m = 0; p = 0;
        make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, dummy_prec_g, xi_high_prec_h, can_xi_high_omega);
        rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, can_xi_high_omega, xi_high, xi_high_mu_array, &can_xi_high_rhoquad);
        can_xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * canlog_det - can_xi_high_rhoquad);
      }
      rho_lp = 0.5 * (DMP * (*H * log_det_g[prior_index]) - rhoquadsum) + xi_low_lp + xi_high_lp;
      can_lp = 0.5 * (DMP * (*H * canlog_det) - canrhoquadsum)  + can_xi_low_lp + can_xi_high_lp;
      
      if (rexp(1) > (rho_lp - can_lp)) {
        for (d = 0; d < *D; d++) {
          for (m = 0; m < *M; m++) {
            for (p = 0; p < *P; p++) {
              prior_index = d + *D * m + *D * *M * p;
              rho[prior_index] = canrho;
              log_det_g[prior_index] = canlog_det;
              rhoquad[prior_index] = canrhoquad_array[prior_index];
              lp[prior_index] = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - rhoquad[prior_index]);
              for (g = 0; g < *G; g++) {
                for (h = 0; h < *G; h++) {
                  k = g + *G * h;
                  prec_g[prior_index + DMP * k] = dummy_prec_g[k];
                }
              }
              for (g = 0; g < *G * *H; g++) {
                for (h = 0; h < *G * *H; h++) {
                  k = prior_index + DMP * (g + *G * *H * h);
                  omega[k] = canomega[k];
                }
              }
            }
          }
        }
        if (xi_low_correlation > 0) {
          xi_low_rho       = canrho;
          xi_low_log_det_g  = canlog_det;
          xi_low_rhoquad   = can_xi_low_rhoquad;
          xi_low_lp        = can_xi_low_lp;
          for (g = 0; g < G2; g++) {
            xi_low_prec_g[g] = dummy_prec_g[g];
          }
          for (g = 0; g < G2 * H2; g++) {
            xi_low_omega[g] = can_xi_low_omega[g];
          }
        }
        if (xi_high_correlation > 0) {
          xi_high_rho       = canrho;
          xi_high_log_det_g  = canlog_det;
          xi_high_rhoquad   = can_xi_high_rhoquad;
          xi_high_lp        = can_xi_high_lp;
          for (g = 0; g < G2; g++) {
            xi_high_prec_g[g] = dummy_prec_g[g];
          }
          for (g = 0; g < G2 * H2; g++) {
            xi_high_omega[g] = can_xi_high_omega[g];
          }
        }
      }
    }
  /******update regression parameters*****/
    for (g = 0; g < *G; g++) {
      for (h = 0; h < *H; h++) {
        cop_index = g + *G * h;
        tail_index = g + *G * h;
        make_b(q_ptr[h], M, kappa, &shape[h], B);
        for (d = 0; d < *D; d++) {
          for (m = m_init; m < *M; m++) {
            for (p = 0; p < *P; p++) {
              
              prior_index = d + *D * m + *D * *M * p;
              prior_index_0 = d + *D * 0 + *D * *M * p;
              
              k0 = g + *G * d + *G * *D * h + *G * *D * *H * (0 + *M * p);
              k = g + *G * d + *G * *D * h + *G * *D * *H * (m + *M * p);
              
              update = 1;
              if (m == 0 && basis_ind[h] > 3) {
                update = 0;
              } /*don't update gamma or weibull cbf*/
              if (m > 0  && p >= *P1) {
              update = 0;
            } /*don't update scale parameters for mean only effects*/
              if (update == 1) {
                z[0] = rnorm(0,1);
                canthetastar[k] = tuning_theta[k] * z[0] + thetastar[k];
                if((m == 1)  & (i >= burn)){
                  z[1] = rnorm(0,1);
                  canthetastar[k0] = thetastar[k0];
                  canthetastar[k0] += tuning_theta[k0] * sqrt(1.0 - corrthetastar[g + *G * d + *G * *D * (h + *H * p)] * (1.0 - corrthetastar[g + *G * d + *G * *D * (h + *H * p)])) * z[1]
                  + tuning_theta[k] * corrthetastar[g + *G * d + *G * *D * (h + *H * p)] * z[0];
                }
                att_thetastar[k] ++;
                threshold(G, J, D_vec, D, H, M, P, &g, &h, &m, &p, timepoints, canthetastar, cantheta);
                if((m == 1)  & (i >= burn)){
                  threshold(G, J, D_vec, D, H, M, P, &g, &h, &zero, &p, timepoints, canthetastar, cantheta);
                  att_thetastar[k0] ++;
                }
                /*evaluate the candidate likelihood*/
                  ll_total = 0;
                can_ll_total = 0;
                for (j = 0; j < *J; j++) {
                  u = g + *G * j + *G * *J * h;
                  clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                            kappa, X, y, y_low, y_high, cantheta, &shape[h],
                            q_ptr[h], p_ptr[h], log_d_ptr[h],
                            AAAA, B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                            &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                            w, M_low, M_high, I_low, I_high,
                            &g, &j, &h);
                  ll_total += ll_sum[u];
                  can_ll_total += can_ll_sum[u];
                }
                if (*cop > 0) {
                  for (n = 0; n < *N; n++) {
                    for (j = 0; j < *J; j++) {
                      q_int = g + *G * j + *G * *J * h + GJH * n;
                      can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                      can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
                    }
                  }
                  copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
                }
                ll_total     +=     cop_ll_sum[g];
                can_ll_total += can_cop_ll_sum[g];
                rho_quad  (G, D, H, M, P, &d, &m, &p, omega, canthetastar, mu_array, &canrhoquad);
                can_lp = 0.5 * (*G * log_det_h[prior_index] + *H * log_det_g[prior_index] - canrhoquad);
                can_lp_0 = 0.0;
                lp_0 = 0.0;
                if((m == 1)  & (i >= burn)){
                  rho_quad (G, D, H, M, P, &d, &zero, &p, omega, canthetastar, mu_array, &canrhoquad);
                  can_lp_0 = 0.5 * (*G * log_det_h[prior_index_0] + *H * log_det_g[prior_index_0] - canrhoquad);
                  lp_0 = lp[prior_index_0];
                }
                
                if (rexp(1) > (ll_total - can_ll_total  + lp[prior_index] + lp_0 - can_lp - can_lp_0)) {
                  thetastar[k] = canthetastar[k];
                  acc_thetastar[k] ++ ;
                  for (n = 0; n < dim_theta; n++) {
                    theta[n] = cantheta[n];
                  }
                  for (n = 0; n < *N; n++) {
                    copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
                    for (j = 0; j < *J; j++) {
                      u               = g + *G * j + *G * *J * h;
                      ll_sum[u]       = can_ll_sum[u];
                      q_int           = u + GJH * n;
                      U[q_int]        = can_U[q_int];
                      W[q_int]        = can_W[q_int];
                      log_like[q_int] = can_log_like[q_int];
                      v               = g + *G * (j * *H + h + JH * n);
                      W2[v]           = can_W2[v];
                    }
                  }
                  cop_ll_sum[g] = can_cop_ll_sum[g];
                  lp[prior_index] = can_lp;
                  rho_quad (G, D, H, M, P, &d, &m, &p, omega, canthetastar, mu_array, &canrhoquad);
                  rhoquad[prior_index] = canrhoquad;
                  if((m == 1)  & (i >= burn)) {
                    thetastar[k0] = canthetastar[k0];
                    acc_thetastar[k0] ++;
                    lp[prior_index_0] = can_lp_0;
                    rho_quad (G, D, H, M, P, &d, &zero, &p, omega, canthetastar, mu_array, &canrhoquad);
                    rhoquad[prior_index] = canrhoquad;
                  }
                } else{ /*restore cantheta to theta*/
                    canthetastar[k] = thetastar[k];
                    if((m == 1)  & (i >= burn)) {
                      canthetastar[k0] = thetastar[k0];
                    }
                    for (n = 0; n < dim_theta; n++) {
                      cantheta[n] = theta[n];
                    }
                    for (n = 0; n < *N; n++) {
                      can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
                      for (j = 0; j < *J; j++) {
                        u                   = g + *G * j + *G * *J * h;
                        can_ll_sum[u]       = ll_sum[u];
                        q_int               = u + GJH * n;
                        can_U[q_int]        = U[q_int];
                        can_W[q_int]        = W[q_int];
                        can_log_like[q_int] = log_like[q_int];
                        v                   = g + *G * (j * *H + h + JH * n);
                        can_W2[v]           = W2[v];
                      }
                    }
                    can_cop_ll_sum[g] = cop_ll_sum[g];
                }
          }
              if (m == 1) {
                i_1000 = (i+1) % 1000;
                if (i_1000 == 0) {
                  i_1000 = 1000;
                }
                j = g + *G * d + *G * *D * (h + *H * p);
                n =  g + *G * d + *G * *D * h + *G * *D * *H * (0 + *M * p);
                k1 = g + *G * d + *G * *D * h + *G * *D * *H * (1 + *M * p);
                
                if (i_1000 % 1000 > 1) {
                  varthetastar[j + GDHP * 0] = ((i_1000 - 2) * varthetastar[j + GDHP * 0] / (i_1000 - 1) + pow(thetastar[n] - meanthetastar[j + GDHP * 0],2) / i_1000);
                  varthetastar[j + GDHP * 1] = ((i_1000 - 2) * varthetastar[j + GDHP * 1] / (i_1000 - 1) + pow(thetastar[k1] - meanthetastar[j + GDHP * 1],2) / i_1000);
                }
                meanthetastar[j + GDHP * 0] = (meanthetastar[j + GDHP * 0] * (i_1000 - 1) +  thetastar[n]) / i_1000; /*mean for m = 0*/
                  meanthetastar[j + GDHP * 1] = (meanthetastar[j + GDHP * 1] * (i_1000 - 1) +  thetastar[k1]) / i_1000; /*mean for m = 1*/
                  meanthetastar[j + GDHP * 2] = (meanthetastar[j + GDHP * 2] * (i_1000 - 1) +  thetastar[n] * thetastar[k1]) / i_1000; /*mean for cross product*/
              }
          }
          }
      }
    }
    }
  if ( ((i+1) % 1000 == 0) ) {
    for (k = 0; k < GDHP; k++) {
      corrthetastar[k] = meanthetastar[k + GDHP * 2] -  meanthetastar[k + GDHP * 0] * meanthetastar[k + GDHP * 1];
      if((varthetastar[k + GDHP * 0] > 0.0001) & (varthetastar[k + GDHP * 1] > 0.0001)){
        corrthetastar[k] /= sqrt(varthetastar[k + GDHP * 0] * varthetastar[k + GDHP * 1]);
      } else {
        corrthetastar[k] = 0.0;
      }
      corrthetastar[k] *= (1000.0 / 999.0);
    }
    for (k = 0; k < GDHP * 3; k++) {
      meanthetastar[k] = 0.0;
    }
    for (k = 0; k < GDHP * 2; k++) {
      varthetastar[k] = 0.0;
    }
  }
  /***** update copula parameters ***/
    if (*cop > 0) {
      for (g = 0; g < *G; g++) {
        if (*cop == 1 || *cop == 3) { /*use if there are random intercepts or slopes*/
            for (n = 0; n < n_g[g]; n++) {
              for (h = 0; h < HQ; h++) {
                dummy_vec_HQ[h] = 0.0;
              }
              for (j = 0; j < JH; j++) {
                k = g + *G * (j + JH * n);
                dummy_vec_JH[j]      = DD_inv[k] * W2[k] - eta[k];
              }
              for (h = 0; h < HQ; h++) {
                for (j = 0; j < JH; j++) {
                  dummy_vec_HQ[h] += Z[g + *G * (j + JH * n + JH * *N * h)] * dummy_vec_JH[j];
                }
              }
              for (u = 0; u < HQ; u++) {
                for (v = 0; v < HQ; v++) {
                  k = u + HQ * v;
                  Omega_gamma[k] = Z2_sum[g + *G * (n + *N * k)];
                }
              }
              for (u = 0; u < HQ; u++) {
                Omega_gamma[u + HQ * u] += delta_prec[g + *G * u];
              }
              InvertSymmetricMatrix(&HQ, Omega_gamma, Omega_gamma_inv, &dummy_log_det);
              CholMatrix(&HQ, Omega_gamma_inv, Omega_gamma_inv_chol, &dummy_log_det);
              for (h = 0; h < HQ; h++) {
                omega_gamma[h] = 0;
                z_HQ[h] = rnorm(0,1);
                for (j = 0; j < HQ; j++) {
                  omega_gamma[h] += Omega_gamma_inv[h + HQ * j] * dummy_vec_HQ[j];
                }
              }
              for (h = 0; h < HQ; h++) {
                k = g + *G * n + *G * *N * h;
                gamma[k] = omega_gamma[h];
                for (u = 0; u < (h + 1); u++) {
                  gamma[k] += Omega_gamma_inv_chol[h + HQ * u] * z_HQ[u];
                }
              }
            }
          copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]);
          for (n = 0; n < *N; n++) {
            can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
          }
          can_cop_ll_sum[g] = cop_ll_sum[g];
          
          /*update delta*/
            for (q_int = 0; q_int < *Q_int; q_int++) {
              for (h = 0; h < *H; h++) {
                k = q_int + *Q_int * h;
                can_delta_prec[g + *G * k] = exp(log(delta_prec[g + *G * k]) + tuning_delta[g + *G * k] * rnorm(0,1));
                can_delta[g + *G * k] = 1 / can_delta_prec[g + *G * k];
                make_D(G, J, H, Q_int, n_g, N, Z, can_delta, W_cov, &g, status, can_DD_inv, can_DD);
                copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, can_DD, can_DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
                
                delta_ll =          cop_ll_sum[g];
                can_delta_ll =  can_cop_ll_sum[g];
                
                for (n = 0; n < n_g[g]; n++) {
                  delta_ll            +=     0.5 * (log(    delta_prec[g + *G * k]) -      delta_prec[g + *G * k] * (pow(gamma[g + *G * (n + *N * k)],2)));
                  can_delta_ll        +=     0.5 * (log(can_delta_prec[g + *G * k]) -  can_delta_prec[g + *G * k] * (pow(gamma[g + *G * (n + *N * k)],2)));
                }
                delta_lp        = (zeta_A[g + *G * k] - 1) * log(    delta_prec[g + *G * k]) -     delta_prec[g + *G * k] / zeta_B[g + *G * k];
                can_delta_lp    = (zeta_A[g + *G * k] - 1) * log(can_delta_prec[g + *G * k]) - can_delta_prec[g + *G * k] / zeta_B[g + *G * k];
                
                if (rexp(1) > delta_ll + delta_lp - can_delta_ll - can_delta_lp) {
                  delta_prec[g + *G * k] = can_delta_prec[g + *G * k];
                  delta[g + *G * k] = can_delta[g + *G * k];
                  acc_delta[g + *G * k] ++;
                  cop_ll_sum[g] = can_cop_ll_sum[g];
                  for (n = 0; n < n_g[g]; n++) {
                    copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
                    for (j = 0; j < JH; j++) {
                      k = g + *G * j + *G * JH * n;
                      DD[k] = can_DD[k];
                      DD_inv[k] = can_DD_inv[k];
                    }
                  }
                } else{
                  can_delta_prec[g + *G * k] = delta_prec[g + *G * k];
                  can_delta[g + *G * k] = delta[g + *G * k];
                  can_cop_ll_sum[g] = cop_ll_sum[g];
                  for (n = 0; n < n_g[g]; n++) {
                    can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
                    for (j = 0; j < JH; j++) {
                      k = g + *G * j + *G * JH * n;
                      can_DD[k] = DD[k];
                      can_DD_inv[k] = DD_inv[k];
                    }
                  }
                }
              }
            }
        }
        if (*cop > 1) {
          /*update eta*/
            for (n = 0; n < n_g[g]; n++) {
              for (j = 0; j < JH; j++) {
                k = g + *G * (j + JH * n);
                dummy_vec_JH[j] = 0.0;
                for (h = 0; h < HQ; h++) {
                  dummy_vec_JH[j] -= Z[k + *G * JH * *N * h] * gamma[g + *G * n + *G * *N * h];
                }
                dummy_vec_JH[j] += DD_inv[k] * W2[k];
              }
              for (j = 0; j < JH; j++) {
                omega_eta[j] = 0;
                z_JH[j] = rnorm(0,1);
                for (k = 0; k < JH; k++) {
                  omega_eta[j] += Omega_eta_inv[g + *G * (j + JH * k)] * dummy_vec_JH[k];
                }
              }
              /*I just forced mean of eta to be 0 for each obs
              if I want to change later get rid of dummy_log_det in this part of code*/
                dummy_log_det = 0.0;
                for (j = 0; j < JH; j++) {
                  k = g + *G * (j + JH * n);
                  if (status_2[k] == -99) {
                    eta[k] = 0.0;
                  } else{
                    eta[k] = omega_eta[j];
                    for (u = 0; u < (j + 1); u++) {
                      eta[k] += Omega_eta_inv_chol[g + *G * (j + JH * u)] * z_JH[u];
                    }
                    if (*cop == 3) {
                      dummy_log_det += eta[k];
                    }
                  }
                }
                /*get rid of below me later if this does not work*/
                  dummy_log_det /= JH;
                  for (j = 0; j < JH; j++) {
                    k = g + *G * (j + JH * n);
                    eta[k] -= dummy_log_det;
                  }
                  /*get rid of above me later if this does not work*/
            }
          copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]);
          for (n = 0; n < *N; n++) {
            can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
          }
          can_cop_ll_sum[g] = cop_ll_sum[g];
          /*** update alpha ***/
            if (*J > 1) {
              can_alpha[g] = runif (0,1);
              for (u = 0; u < *J; u++) {
                for (v = 0; v < *J; v++) {
                  k = u + *J * v;
                  cov_J[k] = pow(can_alpha[g], dist_J[u + *J * v]);
                }
              }
              CholMatrix(J, cov_J, prec_J, &dummy_log_det);
              for (u = 0; u < *J; u++) {
                for (v = 0; v < *J; v++) {
                  k = u + *J * v;
                  can_AR_chol[g + *G * k] = prec_J[k];
                }
              }
              if (dead_sum > 0) {
                longitudinal_log_det(J, cov_J, prec_J, dummy_mat_J2);
              }
              InvertSymmetricMatrix(J, cov_J, prec_J, &can_alpha_log_det[g]);
              for (u = 0; u < *J; u++) {
                for (v = 0; v < *J; v++) {
                  k = u + *J * v;
                  can_AR[g + *G * k]      = cov_J[k];
                  can_AR_prec[g + *G * k] = prec_J[k];
                }
              }
              for (h = 0; h < H2; h++) {
                cov_H[h] = lambda[g + *G * h];
                prec_H[h] = lambda_prec[g + *G * h];
              }
              Kronecker(J, J, H, H, cov_J, cov_H, cov_JH);
              Kronecker(J, J, H, H, prec_J, prec_H, prec_JH);
              for (u = 0; u < JH; u++) {
                for (v = 0; v < JH; v++) {
                  k = u + JH * v;
                  can_eta_cov[g + *G * k] = cov_JH[k];
                  can_eta_prec[g + *G * k] = prec_JH[k];
                  can_W_cov[g + *G * k] = can_eta_cov[g + *G * k];
                }
              }
              for (j = 0; j < JH; j++) {
                can_W_cov[g + *G * (j + JH * j)] += 1.0;
              }
              make_D(G, J, H, Q_int, n_g, N, Z, delta, can_W_cov, &g, status, can_DD_inv, can_DD);
              copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, can_DD, can_DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
              alpha_ll = cop_ll_sum[g];
              can_alpha_ll = can_cop_ll_sum[g];
              
              for (n = 0; n < n_g[g]; n++) {
                qf1 = 0;
                qf2 = 0;
                for (j = 0; j < JH; j++) {
                  for (k = 0; k < JH; k++) {
                    qf1 += eta[g + *G * (j + JH * n)] * eta[g + *G * (k + JH * n)] *     eta_prec[g + *G * (j + JH * k)];
                    qf2 += eta[g + *G * (j + JH * n)] * eta[g + *G * (k + JH * n)] * can_eta_prec[g + *G * (j + JH * k)];
                  }
                }
                if (dead_index[g + *G * n] == *J) {
                  alpha_ll        += 0.5 * (*H *     alpha_log_det[g] - qf1);
                  can_alpha_ll    += 0.5 * (*H * can_alpha_log_det[g] - qf2);
                } else{
                  alpha_ll        -= 0.5 * qf1;
                  can_alpha_ll    -= 0.5 * qf2;
                  qf1 = 0.0;
                  qf2 = 0.0;
                  /*need to change this to proper index*/
                    for (j = 0; j < dead_index[g + *G * n]; j++) {
                      qf1 +=     alpha_det_matrix[g + *G * (start_index[g + *G * n] + *J * j)];
                      qf2 +=     dummy_mat_J2[start_index[g + *G * n] + *J * j];
                    }
                  alpha_ll        += 0.5 * *H * qf1;
                  can_alpha_ll    += 0.5 * *H * qf2;
                }
              }
              if (rexp(1) > alpha_ll - can_alpha_ll) {
                alpha[g] = can_alpha[g];
                alpha_log_det[g] = can_alpha_log_det[g];
                for (u = 0; u < *J; u++) {
                  for (v = 0; v < *J; v++) {
                    k = g + *G * (u + *J * v);
                    AR[k] = can_AR[k];
                    AR_chol[k] = can_AR_chol[k];
                    AR_prec[k] = can_AR_prec[k];
                  }
                }
                for (j = 0; j < J2; j++) {
                  alpha_det_matrix[g + *G * j] = dummy_mat_J2[j];
                }
                for (u = 0; u < JH; u++) {
                  prec_JH[u + JH * u] ++;
                  for (v = 0; v < JH; v++) {
                    k = g + *G * (u + JH * v);
                    eta_cov[k] = can_eta_cov[k];
                    eta_prec[k] = can_eta_prec[k];
                    W_cov[k] = can_W_cov[k];
                  }
                }
                InvertSymmetricMatrix(&JH, prec_JH, cov_JH, &dummy_log_det);
                
                for (u = 0; u < JH * JH; u++) {
                  Omega_eta_inv[g + *G * u] = cov_JH[u];
                }
                CholMatrix(&JH, cov_JH, prec_JH, &dummy_log_det);
                for (u = 0; u < JH * JH; u++) {
                  Omega_eta_inv_chol[g + *G * u] = prec_JH[u];
                }
                cop_ll_sum[g] = can_cop_ll_sum[g];
                for (n = 0; n < n_g[g]; n++) {
                  copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
                  for (j = 0; j < JH; j++) {
                    k = g + *G * j + *G * JH * n;
                    DD[k] = can_DD[k];
                    DD_inv[k] = can_DD_inv[k];
                  }
                }
              } else{
                can_alpha[g] = alpha[g];
                can_alpha_log_det[g] = alpha_log_det[g];
                for (u = 0; u < *J; u++) {
                  for (v = 0; v < *J; v++) {
                    k = g + *G * (u + *J * v);
                    can_AR[k] = AR[k];
                    can_AR_chol[k] = AR_chol[k];
                    can_AR_prec[k] = AR_prec[k];
                  }
                }
                for (u = 0; u < JH; u++) {
                  for (v = 0; v < JH; v++) {
                    k = g + *G * (u + JH * v);
                    can_eta_cov[k] = eta_cov[k];
                    can_eta_prec[k] = eta_prec[k];
                    can_W_cov[k] = W_cov[k];
                  }
                }
                can_cop_ll_sum[g] = cop_ll_sum[g];
                for (n = 0; n < n_g[g]; n++) {
                  can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
                  for (j = 0; j < JH; j++) {
                    k = g + *G * j + *G * JH * n;
                    can_DD[k] = DD[k];
                    can_DD_inv[k] = DD_inv[k];
                  }
                }
              }
            }
          /*update lambda*/
            for (h = 0; h < H2; h++) {
              sai_h[h] = (tuning_nu[g]  - *H - 1) * lambda[g + *G * h];
            }
          RandomInverseWishart(H, &tuning_nu[g], sai_h, cov_H);
          CholMatrix(H, cov_H, prec_H, &dummy_log_det);
          for (h = 0; h < H2; h++) {
            can_lambda_chol[g + *G * h] = prec_H[u + *H * v];
          }
          InvertSymmetricMatrix(H, cov_H, prec_H, &can_lambda_log_det[g]);
          for (h = 0; h < H2; h++) {
            can_lambda[g + *G * h] = cov_H[h];
            can_lambda_prec[g + *G * h] = prec_H[h];
          }
          for (j = 0; j < J2; j++) {cov_J[j] = AR[g + *G * j]; prec_J[j] = AR_prec[g + *G * j];}
          Kronecker(J, J, H, H, cov_J, cov_H, cov_JH);
          Kronecker(J, J, H, H, prec_J, prec_H, prec_JH);
          for (j = 0; j < JH * JH; j++) {
            k = g + *G * j;
            can_eta_cov[k] = cov_JH[j];
            can_eta_prec[k] = prec_JH[j];
            can_W_cov[k] = can_eta_cov[k];
          }
          for (j = 0; j < JH; j++) {
            can_W_cov[g + *G * (j + JH * j)] ++;
          }
          make_D(G, J, H, Q_int, n_g, N, Z, delta, can_W_cov, &g, status, can_DD_inv, can_DD);
          copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, can_DD, can_DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
          
          lambda_ll       =     cop_ll_sum[g];
          can_lambda_ll   = can_cop_ll_sum[g];
          
          qf1 = 0;
          qf2 = 0;
          
          for (n = 0; n < n_g[g]; n++) {
            for (j = 0; j < JH; j++) {
              for (k = 0; k < JH; k++) {
                qf1 += eta[g + *G * (j + JH * n)] * eta[g + *G * (k + JH * n)] *     eta_prec[g + *G * (j + JH * k)];
                qf2 += eta[g + *G * (j + JH * n)] * eta[g + *G * (k + JH * n)] * can_eta_prec[g + *G * (j + JH * k)];
              }
            }
          }
          lambda_ll        -= 0.5 * qf1;
          can_lambda_ll    -= 0.5 * qf2;
          
          lambda_ll           += 0.5 * (*J * n_g[g] - dead_sum[g]) *     lambda_log_det[g];
          can_lambda_ll       += 0.5 * (*J * n_g[g] - dead_sum[g]) * can_lambda_log_det[g];
          
          qf1 = 0;
          qf2 = 0;
          for (h = 0; h < *H; h++) {
            for (u = 0; u < *H; u++) {
              qf1 += sai_lambda_0[g + *G * (h + *H * u)] *     lambda_prec[g + *G * (u + *H * h)];
              qf2 += sai_lambda_0[g + *G * (h + *H * u)] * can_lambda_prec[g + *G * (u + *H * h)];
            }
          }
          lambda_lp       = 0.5 * (nu_lambda[g] *     lambda_log_det[g] - qf1);
          can_lambda_lp   = 0.5 * (nu_lambda[g] * can_lambda_log_det[g] - qf2);
          
          qf1 = 0;
          qf2 = 0;
          for (h = 0; h < *H; h++) {
            for (u = 0; u < *H; u++) {
              qf1 +=     lambda[g + *G * (h + *H * u)] * can_lambda_prec[g + *G * (u + *H * h)];
              qf2 += can_lambda[g + *G * (h + *H * u)] *     lambda_prec[g + *G * (u + *H * h)];
            }
          }
          t2 = 0.5 * (tuning_nu[g] * (*H * log(tuning_nu[g] - *H - 1)  - lambda_log_det[g])
                      + (tuning_nu[g] + *H + 1) * can_lambda_log_det[g]
                      - (tuning_nu[g] - *H - 1) * qf1);
          t1 = 0.5 * (tuning_nu[g] * (*H * log(tuning_nu[g] - *H - 1)  - can_lambda_log_det[g])
                      + (tuning_nu[g] + *H + 1) * lambda_log_det[g]
                      - (tuning_nu[g] - *H - 1) * qf2);
          if (rexp(1) > lambda_ll - can_lambda_ll + lambda_lp - can_lambda_lp + t2 - t1) {
            lambda_log_det[g] = can_lambda_log_det[g];
            acc_lambda[g] ++;
            for (h = 0; h < H2; h++) {
              k = g + *G * h;
              lambda[k] = can_lambda[k];
              lambda_chol[k] = can_lambda_chol[k];
              lambda_prec[k] = can_lambda_prec[k];
            }
            for (j = 0; j < JH * JH; j++) {
              eta_cov[g + *G * j] = can_eta_cov[g + *G * j];
              eta_prec[g + *G * j] = can_eta_prec[g + *G * j];
              W_cov[g + *G * j] = can_W_cov[g + *G * j];
            }
            cop_ll_sum[g] = can_cop_ll_sum[g];
            for (n = 0; n < n_g[g]; n++) {
              copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
              for (j = 0; j < JH; j++) {
                k = g + *G * j + *G * JH * n;
                DD[k] = can_DD[k];
                DD_inv[k] = can_DD_inv[k];
              }
            }
            for (u = 0; u < JH; u++) {
              prec_JH[u + JH * u] ++;
            }
            InvertSymmetricMatrix(&JH, prec_JH, cov_JH, &dummy_log_det);
            for (u = 0; u < JH * JH; u++) {
              Omega_eta_inv[g + *G * u] = cov_JH[u];
            }
            CholMatrix(&JH, cov_JH, prec_JH, &dummy_log_det);
            for (u = 0; u < JH * JH; u++) {
              Omega_eta_inv_chol[g + *G * u] = prec_JH[u];
            }
          } else{
            can_lambda_log_det[g] = lambda_log_det[g];
            for (h = 0; h < H2; h++) {
              k = g + *G * h;
              can_lambda[k] = lambda[k];
              can_lambda_chol[k] = lambda_chol[k];
              can_lambda_prec[k] = lambda_prec[k];
            }
            for (j = 0; j < JH * JH; j++) {
              can_eta_cov[g + *G * j] = eta_cov[g + *G * j];
              can_eta_prec[g + *G * j] = eta_prec[g + *G * j];
              can_W_cov[g + *G * j] = W_cov[g + *G * j];
            }
            can_cop_ll_sum[g] = cop_ll_sum[g];
            for (n = 0; n < n_g[g]; n++) {
              can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
              for (j = 0; j < JH; j++) {
                k = g + *G * j + *G * JH * n;
                can_DD[k] = DD[k];
                can_DD_inv[k] = DD_inv[k];
              }
            }
          }
        }
      }
    }
  /************** update shape parameter for basis functions **************/
    for (h = 0; h < *H; h++) {
      if (basis_ind[h] == 1 || basis_ind[h] == 4 || basis_ind[h] == 5 ) {
        canshape =  exp(tuning_shape[h] * rnorm(0,1) +  log(shape[h]));
        /*construct the candidate basis matrix*/
          make_b(q_ptr[h], M, kappa, &canshape, can_B);
        /*candidate likelihood*/
          ll_total = 0;
          can_ll_total = 0;
          for (g = 0; g < *G; g++) {
            cop_index = g + *G * h;
            tail_index = g + *G * h;
            for (j = 0; j < *J; j++) {
              u = g + *G * j + *G * *J * h;
              clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                        kappa, X, y, y_low, y_high, theta, &canshape,
                        q_ptr[h], p_ptr[h], log_d_ptr[h],
                        AAAA, can_B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                        &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                        w, M_low, M_high, I_low, I_high,
                        &g, &j, &h);
              ll_total += ll_sum[u];
              can_ll_total += can_ll_sum[u];
            }
            if (*cop > 0) {
              for (n = 0; n < *N; n++) {
                for (j = 0; j < *J; j++) {
                  q_int = g + *G * j + *G * *J * h + GJH * n;
                  can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                  can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
                }
              }
              copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
            }
            ll_total        +=     cop_ll_sum[g];
            can_ll_total    += can_cop_ll_sum[g];
          }
          dummy_lp =     -0.5 * shapeprec[h] * pow(log(shape[h]) - shape_mean[h], 2);
          dummy_can_lp = -0.5 * shapeprec[h] * pow(log(canshape) - shape_mean[h], 2);
          if (rexp(1) > (ll_total - can_ll_total  + dummy_lp - dummy_can_lp)) {
            shape[h] = canshape;
            acc_shape[h] ++;
            for (g = 0; g < *G; g++) {
              cop_ll_sum[g] = can_cop_ll_sum[g];
              for (n = 0; n < *N; n++) {
                copula_loglike[g + *G * n]  = can_copula_loglike[g + *G * n];
                for (j = 0; j < *J; j++) {
                  u               = g + *G * j + *G * *J * h;
                  ll_sum[u]       = can_ll_sum[u];
                  q_int           = u + GJH * n;
                  U[q_int]      = can_U[q_int];
                  W[q_int]        = can_W[q_int];
                  log_like[q_int] = can_log_like[q_int];
                  v               = g + *G * (j * *H + h + JH * n);
                  W2[v]           = can_W2[v];
                }
              }
            }
          } else{
            for (g = 0; g < *G; g++) {
              can_cop_ll_sum[g] = cop_ll_sum[g];
              for (n = 0; n < *N; n++) {
                can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
                for (j = 0; j < *J; j++) {
                  u = g + *G * j + *G * *J * h;
                  can_ll_sum[u]       = ll_sum[u];
                  q_int               = u + GJH * n;
                  can_U[q_int]      = U[q_int];
                  can_W[q_int]        = W[q_int];
                  can_log_like[q_int] = log_like[q_int];
                  v                   = g + *G * (j * *H + h + JH * n);
                  can_W2[v]           = W2[v];
                }
              }
            }
          }
      }
      if (basis_ind[h] == 3) {
        logit_shape = log(shape[h]) - log(1 - shape[h]);
        can_logit_shape =  logit_shape + tuning_shape[h] * rnorm(0,1);
        canshape = 1 / (1 + exp(-can_logit_shape));
        /*construct the candidate basis matrix*/
          make_b(q_ptr[h], M, kappa, &canshape, can_B);
        ll_total = 0;
        can_ll_total = 0;
        for (g = 0; g < *G; g++) {
          cop_index = g + *G * h;
          tail_index = g + *G * h;
          for (j = 0; j < *J; j++) {
            u = g + *G * j + *G * *J * h;
            clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                      kappa, X, y, y_low, y_high, theta, &canshape,
                      q_ptr[h], p_ptr[h], log_d_ptr[h],
                      AAAA, can_B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                      &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                      w, M_low, M_high, I_low, I_high,
                      &g, &j, &h);
            ll_total += ll_sum[u];
            can_ll_total += can_ll_sum[u];
          }
          if (*cop > 0) {
            for (n = 0; n < *N; n++) {
              for (j = 0; j < *J; j++) {
                q_int = g + *G * j + *G * *J * h + GJH * n;
                can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
              }
            }
            copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
          }
          ll_total        +=     cop_ll_sum[g];
          can_ll_total    += can_cop_ll_sum[g];
        }
        dummy_lp =     -0.5 * shapeprec[h] * pow(log(   logit_shape) - shape_mean[h], 2);
        dummy_can_lp = -0.5 * shapeprec[h]* pow(log(can_logit_shape) - shape_mean[h], 2);
        if (rexp(1) > (ll_total - can_ll_total  + dummy_lp - dummy_can_lp)) {
          shape[h] = canshape;
          acc_shape[h] ++;
          for (g = 0; g < *G; g++) {
            cop_index = g + *G * h;
            cop_ll_sum[g] = can_cop_ll_sum[g];
            for (n = 0; n < *N; n++) {
              copula_loglike[g + *G * n]  = can_copula_loglike[g + *G * n];
              for (j = 0; j < *J; j++) {
                u               = g + *G * j + *G * *J * h;
                ll_sum[u]       = can_ll_sum[u];
                q_int           = u + GJH * n;
                U[q_int]      = can_U[q_int];
                W[q_int]        = can_W[q_int];
                log_like[q_int] = can_log_like[q_int];
                v               = g + *G * (j * *H + h + JH * n);
                W2[v]           = can_W2[v];
              }
            }
          }
        } else{
          can_cop_ll_sum[g] = cop_ll_sum[g];
          for (n = 0; n < *N; n++) {
            can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
            for (j = 0; j < *J; j++) {
              u                   = g + *G * j + *G * *J * h;
              can_ll_sum[u]       = ll_sum[u];
              q_int               = u + GJH * n;
              can_U[q_int]        = U[q_int];
              can_W[q_int]        = W[q_int];
              can_log_like[q_int] = log_like[q_int];
              v                   = g + *G * (j * *H + h + JH * n);
              can_W2[v]           = W2[v];
            }
          }
        }
      }
    }
  /*********************   update tails             ************/
    if (*tau_low > 0.0) {
      if (update_pareto_low) {
        can_pareto = 1 - pareto_low;
        ll_total = 0;
        can_ll_total = 0;
        for (g = 0; g < *G; g++) {
          for (h = 0; h < *H; h++) {
            k = g + *G * h;
            cop_index = g + *G * h;
            tail_index = g + *G * h;
            make_b(q_ptr[h], M, kappa, &shape[h], B);
            for (j = 0; j < *J; j++) {
              u = g + *G * j + *G * *J * h;
              clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                        kappa, X, y, y_low, y_high, theta, &shape[h],
                        q_ptr[h], p_ptr[h], log_d_ptr[h],
                        AAAA, B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                        &can_pareto, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                        w, M_low, M_high, I_low, I_high,
                        &g, &j, &h);
              ll_total += ll_sum[u];
              can_ll_total += can_ll_sum[u];
            }
            if (*cop > 0) {
              for (n = 0; n < *N; n++) {
                for (j = 0; j < *J; j++) {
                  q_int = g + *G * j + *G * *J * h + GJH * n;
                  can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                  can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
                }
              }
              copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
            }
            ll_total        +=     cop_ll_sum[g];
            can_ll_total    += can_cop_ll_sum[g];
          }
        }
        dummy_log_det = ll_total;
        if (can_ll_total < ll_total) {
          dummy_log_det = can_ll_total;
        }
        ll_total -= dummy_log_det;
        can_ll_total -= dummy_log_det;
        if (pareto_low == 0) {
          prob_pareto_stay = exp(ll_total) * (1 - pareto_prob_low) / (exp(ll_total) * (1 - pareto_prob_low) + exp(can_ll_total) * (pareto_prob_low));
        }else{
          prob_pareto_stay = exp(ll_total) * pareto_prob_low / (exp(ll_total) * (pareto_prob_low) + exp(can_ll_total) * (1 - pareto_prob_low));
        }
        random_uniform = runif (0,1);
        if (random_uniform > prob_pareto_stay) {
          pareto_low = 1 - pareto_low;
          for (g = 0; g < *G; g++) {
            cop_ll_sum[g] = can_cop_ll_sum[g];
            for (h = 0; h < *H; h++) {
              for (n = 0; n < *N; n++) {
                copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
                for (j = 0; j < *J; j++) {
                  u = g + *G * j + *G * *J * h;
                  ll_sum[u]       = can_ll_sum[u];
                  q_int           = u + GJH * n;
                  U[q_int]        = can_U[q_int];
                  W[q_int]        = can_W[q_int];
                  log_like[q_int] = can_log_like[q_int];
                  v               = g + *G * (j * *H + h + JH * n);
                  W2[v]           = can_W2[v];
                }
              }
            }
          }
        }else{
          for (g = 0; g < *G; g++) {
            can_cop_ll_sum[g] = cop_ll_sum[g];
            for (h = 0; h < *H; h++) {
              for (n = 0; n < *N; n++) {
                can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
                for (j = 0; j < *J; j++) {
                  u = g + *G * j + *G * *J * h;
                  can_ll_sum[u]       = ll_sum[u];
                  q_int               = u + GJH * n;
                  can_U[q_int]        = U[q_int];
                  can_W[q_int]        = W[q_int];
                  can_log_like[q_int] = log_like[q_int];
                  v                   = g + *G * (j * *H + h + JH * n);
                  can_W2[v]           = W2[v];
                }
              }
            }
          }
        }
      }
      /*update the prior for xi_low*/
        d = 0;
        m = 0;
        p = 0;
        thetastar_index = m + *M * p;
        prior_index = d + *D * m + *D * *M * p;
        this_thetastar = *G * *D * *H * thetastar_index;
        
        if (GH > 1) {
          if (*H == 1) {
            h = 0;
            /*** update mean via Gibbs sampling ***/
              mu_var = 0.0;
              for (u = 0; u < *G; u++) {
                for (v = 0; v < *G; v++) {
                  mu_var += xi_low_prec_g[u + *G * v];
                }
              }
              mu_var *= xi_low_prec_h[0];
              mu_var += xi_low_prec_mu_0;
              mu_mean = 0.0;
              for (u = 0; u < *G; u++) {
                for (v = 0; v < *G; v++) {
                  mu_mean += 1 * xi_low[v] * xi_low_prec_g[u + *G * v];
                }
              }
              mu_mean *= xi_low_prec_h[0];
              mu_mean += xi_low_mu_0;
              mu_var = 1 / mu_var;
              mu_mean *= mu_var;
              mu_sd = sqrt(mu_var);
              xi_low_mu[0] = mu_sd * rnorm(0,1) + mu_mean;
              for (g = 0; g < *G; g++) {
                xi_low_mu_array[g] = xi_low_mu[0];
              }
              rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
              xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - xi_low_rhoquad);
              /*** update prec_h via Gibbs sampling ***/
                dummy_log_det  = xi_low_rhoquad / xi_low_prec_h[0]; /*remove precision term from quadratic form*/
                xi_low_prec_h[0] = rgamma(*G / 2 + xi_low_sig_a, (2 * xi_low_sig_b)/(xi_low_sig_b * dummy_log_det + 2));
              xi_low_sig_h[0]  = 1 / xi_low_prec_h[0];
              make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_prec_g, xi_low_prec_h, xi_low_omega);
              rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
              xi_low_log_det_h = log(xi_low_prec_h[0]);
              xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - xi_low_rhoquad);
          }
          if (*G == 1) {
            /*** update mean via Gibbs sampling ***/
              mu_mean = xi_low[0];
              for (h = 1; h < *H; h++) {
                mu_mean += xi_low[h];
              }
              mu_mean *= xi_low_prec_h[0];
              mu_mean += xi_low_mu_0;
              mu_var = *H * xi_low_prec_h[0] + xi_low_prec_mu_0;
              mu_var = 1 / mu_var;
              mu_mean *= mu_var;
              mu_sd = sqrt(mu_var);
              xi_low_mu[0] = mu_sd * rnorm(0,1) + mu_mean;
              for (h = 0; h < *H; h++) {
                xi_low_mu_array[h] = xi_low_mu[0];
              }
              rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
              /*** update prec_h via Gibbs sampling ***/
                dummy_log_det = xi_low_rhoquad / xi_low_prec_h[0]; /*remove precision term from quadratic form*/
                xi_low_prec_h[0] = rgamma(*H / 2 + xi_low_sig_a, (2 * xi_low_sig_b)/(xi_low_sig_b * dummy_log_det + 2));
              xi_low_sig_h[0] = 1 / xi_low_prec_h[0];
              make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_prec_g, xi_low_prec_h, xi_low_omega);
              rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
              xi_low_log_det_h = log(xi_low_prec_h[0]);
              xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - xi_low_rhoquad);
          }
          if (*G > 1 && *H > 1) {
            /*** update mean via Gibbs sampling ***/
              for (h = 0; h < *H; h++) {
                z_mu[h] = rnorm(0,1);
              }
            gibbs_mu(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low, xi_low_prec_g, xi_low_prec_h, z_mu, &xi_low_mu_0, &xi_low_prec_mu_0, xi_low_mu, xi_low_mu_array);
            /*** update sig_h via Gibbs sampling ***/
              /*if *G > 1 then we have replication so we can use inverse wishart*/
              make_sai(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low, xi_low_prec_g, xi_low_mu_array, sai_h, xi_low_sai_0);
            RandomInverseWishart(H, &xi_low_nu, sai_h, xi_low_sig_h);
            InvertSymmetricMatrix(H, xi_low_sig_h, xi_low_prec_h, &xi_high_log_det_h);
            make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_prec_g, xi_low_prec_h, xi_low_omega);
            rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, xi_low, xi_low_mu_array, &xi_low_rhoquad);
            xi_low_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - xi_low_rhoquad);
          }
        }
        /**** update xi_low ******/
          /**** use only prior if pareto_low == 0 ******/
          if (pareto_low == 0) {
            update_xi_exponential(G, H, xi_low_prec_g, xi_low_sig_h, xi_low_mu_array, xi_low);
          }
        /**** use data and prior if pareto_low == 1 ******/
          if (pareto_low == 1) {
            for (g = 0; g < *G; g++) {
              for (h = 0; h < *H; h++) {
                tail_index = g + *G * h;
                cop_index = g + *G * h;
                make_b(q_ptr[h], M, kappa, &shape[h], B);
                can_xi_low[tail_index] =  tuning_xi_low[tail_index] * rnorm(0,1) +  xi_low[tail_index];
                if (can_xi_low[tail_index] < -10.0) {
                  can_xi_low[tail_index] = -10.0;
                } /*avoid numerical instability*/
                  if (can_xi_low[tail_index] >  2.0) {
                    can_xi_low[tail_index] =  2.0;
                  } /*avoid numerical instability*/
                  ll_total = 0;
                  can_ll_total = 0;
                  for (j = 0; j < *J; j++) {
                    u = g + *G * j + *G * *J * h;
                    clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                              kappa, X, y, y_low, y_high, theta, &shape[h],
                              q_ptr[h], p_ptr[h], log_d_ptr[h],
                              AAAA, B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                              &pareto_low, &pareto_high, tau_low, tau_high, &can_xi_low[tail_index], &xi_high[tail_index],
                              w, M_low, M_high, I_low, I_high,
                              &g, &j, &h);
                    
                    ll_total += ll_sum[u];
                    can_ll_total += can_ll_sum[u];
                  }
                  if (*cop > 0) {
                    for (n = 0; n < *N; n++) {
                      for (j = 0; j < *J; j++) {
                        q_int = g + *G * j + *G * *J * h + GJH * n;
                        can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                        can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
                      }
                    }
                    copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
                  }
                  ll_total        +=     cop_ll_sum[g];
                  can_ll_total    += can_cop_ll_sum[g];
                  rho_quad  (G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_low_omega, can_xi_low, xi_low_mu_array, &canrhoquad);
                  can_lp = 0.5 * (*G * xi_low_log_det_h + *H * xi_low_log_det_g - canrhoquad);
                  if (rexp(1) > (ll_total - can_ll_total  + xi_low_lp - can_lp)) {
                    xi_low_lp = can_lp;
                    xi_low[tail_index] = can_xi_low[tail_index];
                    acc_xi_low[tail_index] ++;
                    cop_ll_sum[g] = can_cop_ll_sum[g];
                    for (n = 0; n < *N; n++) {
                      copula_loglike[g + *G * n]  = can_copula_loglike[g + *G * n];
                      for (j = 0; j < *J; j++) {
                        u = g + *G * j + *G * *J * h;
                        ll_sum[u]       = can_ll_sum[u];
                        q_int           = u + GJH * n;
                        U[q_int]        = can_U[q_int];
                        W[q_int]        = can_W[q_int];
                        log_like[q_int] = can_log_like[q_int];
                        v               = g + *G * (j * *H + h + JH * n);
                        W2[v]           = can_W2[v];
                      }
                    }
                  } else{
                    can_xi_low[tail_index] = xi_low[tail_index];
                    can_cop_ll_sum[g] = cop_ll_sum[g];
                    for (n = 0; n < *N; n++) {
                      can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
                      for (j = 0; j < *J; j++) {
                        u = g + *G * j + *G * *J * h;
                        can_ll_sum[u]       = ll_sum[u];
                        q_int               = u + GJH * n;
                        can_U[q_int]      = U[q_int];
                        can_W[q_int]        = W[q_int];
                        can_log_like[q_int] = log_like[q_int];
                        v                   = g + *G * (j * *H + h + JH * n);
                        can_W2[v]           = W2[v];
                      }
                    }
                  }
              }
            }
          }
    }
  if (*tau_high < 1.0) {
    if (update_pareto_high) {
      can_pareto = 1 - pareto_high;
      ll_total = 0;
      can_ll_total = 0;
      for (g = 0; g < *G; g++) {
        for (h = 0; h < *H; h++) {
          k = g + *G * h;
          cop_index = g + *G * h;
          tail_index = g + *G * h;
          make_b(q_ptr[h], M, kappa, &shape[h], B);
          for (j = 0; j < *J; j++) {
            u = g + *G * j + *G * *J * h;
            clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                      kappa, X, y, y_low, y_high, theta, &shape[h],
                      q_ptr[h], p_ptr[h], log_d_ptr[h],
                      AAAA, B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                      &pareto_low, &can_pareto, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                      w, M_low, M_high, I_low, I_high,
                      &g, &j, &h);
            ll_total += ll_sum[u];
            can_ll_total += can_ll_sum[u];
          }
          if (*cop > 0) {
            for (n = 0; n < *N; n++) {
              for (j = 0; j < *J; j++) {
                q_int = g + *G * j + *G * *J * h + GJH * n;
                can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
              }
            }
            copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
          }
          ll_total        +=     cop_ll_sum[g];
          can_ll_total    += can_cop_ll_sum[g];
        }
      }
      dummy_log_det = ll_total;
      if (can_ll_total < ll_total) {
        dummy_log_det = can_ll_total;
      }
      ll_total -= dummy_log_det;
      can_ll_total -= dummy_log_det;
      if (pareto_high == 0) {
        prob_pareto_stay = exp(ll_total) * (1 - pareto_prob_high) / (exp(ll_total) * (1 - pareto_prob_high) + exp(can_ll_total) * (pareto_prob_high));
      }else{
        prob_pareto_stay = exp(ll_total) * pareto_prob_high/ (exp(ll_total) * (pareto_prob_high) + exp(can_ll_total) * (1 - pareto_prob_high));
      }
      random_uniform = runif (0,1);
      if (random_uniform > prob_pareto_stay) {
        pareto_high = 1 - pareto_high;
        for (g = 0; g < *G; g++) {
          cop_ll_sum[g] = can_cop_ll_sum[g];
          for (h = 0; h < *H; h++) {
            for (n = 0; n < *N; n++) {
              copula_loglike[g + *G * n] = can_copula_loglike[g + *G * n];
              for (j = 0; j < *J; j++) {
                u = g + *G * j + *G * *J * h;
                ll_sum[u]       = can_ll_sum[u];
                q_int           = u + GJH * n;
                U[q_int]        = can_U[q_int];
                W[q_int]        = can_W[q_int];
                log_like[q_int] = can_log_like[q_int];
                v               = g + *G * (j * *H + h + JH * n);
                W2[v]           = can_W2[v];
              }
            }
          }
        }
      }else{
        for (g = 0; g < *G; g++) {
          can_cop_ll_sum[g] = cop_ll_sum[g];
          for (h = 0; h < *H; h++) {
            for (n = 0; n < *N; n++) {
              can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
              for (j = 0; j < *J; j++) {
                u = g + *G * j + *G * *J * h;
                can_ll_sum[u]       = ll_sum[u];
                q_int               = u + GJH * n;
                can_U[q_int]        = U[q_int];
                can_W[q_int]        = W[q_int];
                can_log_like[q_int] = log_like[q_int];
                v                   = g + *G * (j * *H + h + JH * n);
                can_W2[v]           = W2[v];
              }
            }
          }
        }
      }
    }
    /*update the prior for xi_high*/
      d = 0;
      m = 0;
      p = 0;
      thetastar_index = m + *M * p;
      prior_index = d + *D * m + *D * *M * p;
      this_thetastar = *G * *D * *H * thetastar_index;
      
      if (GH > 1) {
        if (*H == 1) {
          h = 0;
          /*** update mean via Gibbs sampling ***/
            mu_var = 0.0;
            for (u = 0; u < *G; u++) {
              for (v = 0; v < *G; v++) {
                mu_var += xi_high_prec_g[u + *G * v];
              }
            }
            mu_var *= xi_high_prec_h[0];
            mu_var += xi_high_prec_mu_0;
            mu_mean = 0.0;
            for (u = 0; u < *G; u++) {
              for (v = 0; v < *G; v++) {
                mu_mean += 1 * xi_high[v] * xi_high_prec_g[u + *G * v];
              }
            }
            mu_mean *= xi_high_prec_h[0];
            mu_mean += xi_high_mu_0;
            mu_var = 1 / mu_var;
            mu_mean *= mu_var;
            mu_sd = sqrt(mu_var);
            xi_high_mu[0] = mu_sd * rnorm(0,1) + mu_mean;
            for (g = 0; g < *G; g++) {
              xi_high_mu_array[g] = xi_high_mu[0];
            }
            rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
            xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - xi_high_rhoquad);
            /*** update prec_h via Gibbs sampling ***/
              dummy_log_det  = xi_high_rhoquad / xi_high_prec_h[0]; /*remove precision term from quadratic form*/
              xi_high_prec_h[0] = rgamma(*G / 2 + xi_high_sig_a, (2 * xi_high_sig_b)/(xi_high_sig_b * dummy_log_det + 2));
            xi_high_sig_h[0]  = 1 / xi_high_prec_h[0];
            make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_prec_g, xi_high_prec_h, xi_high_omega);
            rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
            xi_high_log_det_h = log(xi_high_prec_h[0]);
            xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - xi_high_rhoquad);
        }
        if (*G == 1) {
          /*** update mean via Gibbs sampling ***/
            mu_mean = xi_high[0];
            for (h = 1; h < *H; h++) {
              mu_mean += xi_high[h];
            }
            mu_mean *= xi_high_prec_h[0];
            mu_mean += xi_high_mu_0;
            mu_var = *H * xi_high_prec_h[0] + xi_high_prec_mu_0;
            mu_var = 1 / mu_var;
            mu_mean *= mu_var;
            mu_sd = sqrt(mu_var);
            xi_high_mu[0] = mu_sd * rnorm(0,1) + mu_mean;
            for (h = 0; h < *H; h++) {
              xi_high_mu_array[h] = xi_high_mu[0];
            }
            rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
            /*** update prec_h via Gibbs sampling ***/
              dummy_log_det = xi_high_rhoquad / xi_high_prec_h[0]; /*remove precision term from quadratic form*/
              xi_high_prec_h[0] = rgamma(*H / 2 + xi_high_sig_a, (2 * xi_high_sig_b)/(xi_high_sig_b * dummy_log_det + 2));
            xi_high_sig_h[0] = 1 / xi_high_prec_h[0];
            make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_prec_g, xi_high_prec_h, xi_high_omega);
            rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
            xi_high_log_det_h = log(xi_high_prec_h[0]);
            xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - xi_high_rhoquad);
        }
        if (*G > 1 && *H > 1) {
          /*** update mean via Gibbs sampling ***/
            for (h = 0; h < *H; h++) {
              z_mu[h] = rnorm(0,1);
            }
          gibbs_mu(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high, xi_high_prec_g, xi_high_prec_h, z_mu, &xi_high_mu_0, &xi_high_prec_mu_0, xi_high_mu, xi_high_mu_array);
          /*** update sig_h via Gibbs sampling ***/
            /*if *G > 1 then we have replication so we can use inverse wishart*/
            make_sai(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high, xi_high_prec_g, xi_high_mu_array, sai_h, xi_high_sai_0);
          RandomInverseWishart(H, &xi_high_nu, sai_h, xi_high_sig_h);
          InvertSymmetricMatrix(H, xi_high_sig_h, xi_high_prec_h, &xi_high_log_det_h);
          make_omega(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_prec_g, xi_high_prec_h, xi_high_omega);
          rho_quad(G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, xi_high, xi_high_mu_array, &xi_high_rhoquad);
          xi_high_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - xi_high_rhoquad);
        }
        /*** independent metropolis update the between model correlation ***/
      }
      /**** update xi_high ******/
        /**** use only prior if pareto_high == 0 ******/
        if (pareto_high == 0) {
          update_xi_exponential(G, H, xi_high_prec_g, xi_high_sig_h, xi_high_mu_array, xi_high);
        }
      /**** use data and prior if pareto_high == 1 ******/
        if (pareto_high == 1) {
          /*update xi_high */
            for (g = 0; g < *G; g++) {
              for (h = 0; h < *H; h++) {
                tail_index = g + *G * h;
                cop_index = g + *G * h;
                make_b(q_ptr[h], M, kappa, &shape[h], B);
                can_xi_high[tail_index] =  tuning_xi_high[tail_index] * rnorm(0,1) +  xi_high[tail_index];
                if (can_xi_high[tail_index] < -10.0) {
                  can_xi_high[tail_index] = -10.0;
                }
                if (can_xi_high[tail_index] >  2.0) {
                  can_xi_high[tail_index] =  2.0;
                }
                ll_total = 0;
                can_ll_total = 0;
                
                for (j = 0; j < *J; j++) {
                  u = g + *G * j + *G * *J * h;
                  clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                            kappa, X, y, y_low, y_high, theta, &shape[h],
                            q_ptr[h], p_ptr[h], log_d_ptr[h],
                            AAAA, B, bin, bin_low, bin_high, &can_ll_sum[u], can_log_like, can_U, U_low, U_high,
                            &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &can_xi_high[tail_index],
                            w, M_low, M_high, I_low, I_high,
                            &g, &j, &h);
                  
                  ll_total += ll_sum[u];
                  can_ll_total += can_ll_sum[u];
                }
                if (*cop > 0) {
                  for (n = 0; n < *N; n++) {
                    for (j = 0; j < *J; j++) {
                      q_int = g + *G * j + *G * *J * h + GJH * n;
                      can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
                      can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
                    }
                  }
                  copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, eta, DD, DD_inv, gamma, can_copula_loglike, &can_cop_ll_sum[g]);
                }
                ll_total        +=     cop_ll_sum[g];
                can_ll_total    += can_cop_ll_sum[g];
                rho_quad  (G, &Dummy_one, H, &Dummy_one, &Dummy_one, &d, &m, &p, xi_high_omega, can_xi_high, xi_high_mu_array, &canrhoquad);
                can_lp = 0.5 * (*G * xi_high_log_det_h + *H * xi_high_log_det_g - canrhoquad);
                if (rexp(1) > (ll_total - can_ll_total  + xi_high_lp - can_lp)) {
                  xi_high[tail_index] = can_xi_high[tail_index];
                  acc_xi_high[tail_index] ++;
                  cop_ll_sum[g] = can_cop_ll_sum[g];
                  xi_high_lp = can_lp;
                  for (n = 0; n < *N; n++) {
                    copula_loglike[g + *G * n]  = can_copula_loglike[g + *G * n];
                    for (j = 0; j < *J; j++) {
                      u               = g + *G * j + *G * *J * h;
                      ll_sum[u]       = can_ll_sum[u];
                      q_int           = u + GJH * n;
                      U[q_int]        = can_U[q_int];
                      W[q_int]        = can_W[q_int];
                      log_like[q_int] = can_log_like[q_int];
                      v               = g + *G * (j * *H + h + JH * n);
                      W2[v]           = can_W2[v];
                    }
                  }
                } else{
                  can_xi_high[tail_index] = xi_high[tail_index];
                  can_cop_ll_sum[g] = cop_ll_sum[g];
                  for (n = 0; n < *N; n++) {
                    can_copula_loglike[g + *G * n]  = copula_loglike[g + *G * n];
                    for (j = 0; j < *J; j++) {
                      u                   = g + *G * j + *G * *J * h;
                      can_ll_sum[u]       = ll_sum[u];
                      q_int               = u + GJH * n;
                      can_U[q_int]        = U[q_int];
                      can_W[q_int]        = W[q_int];
                      can_log_like[q_int] = log_like[q_int];
                      v                   = g + *G * (j * *H + h + JH * n);
                      can_W2[v]           = W2[v];
                    }
                  }
                }
              }
            }
        }
  }
  if (missing_sum > 0) {
    if (*cop == 0) {
      for (g = 0; g < *G; g++) {
        for (n = 0; n < n_g[g]; n++) {
          if (missing_index[g + *G * n] > 0) {
            for (j = 0; j < JH; j++) {
              k = g + *G * j + *G * JH * n;
              if (status[k] == 0) {
                U[k] = runif (0,1);
                can_U[k] = U[k];
              }
            }
          }
        }
      }
    }
    /*definitely need to check this*/
      if (*cop > 0) {
        for (g = 0; g < *G; g++) {
          for (j = 0; j < *J; j++) {
            for (h = 0; h < *H; h++) {
              for (n = 0; n < n_g[g]; n++) {
                k = g + *G * j + *G * *J * (h + *H * n);
                v = g + *G * (j * *H + h + JH * n);
                if (status[k] == 0) {
                  dummy_log_det = 0.0;
                  for (q_int = 0; q_int < HQ; q_int++) {
                    dummy_log_det += Z[g + *G * (j + *J * h) + *G * JH * (n + *N * q_int)] * gamma[g + *G * (n + *N * q_int)];
                  }
                  dummy_log_det += eta[v] + rnorm(0,1);
                  dummy_log_det *= DD[v];
                  W2[v] = dummy_log_det;
                  can_W2[v] = W2[v];
                  W[k] = W2[v];
                  can_W[k] = W[k];
                  PNorm(&W[k], &mn, &scale, &scale, &U[k]);
                  can_U[k] = U[k];
                }
              }
            }
          }
        }
      }
    for (g = 0; g < *G; g++) {
      for (j = 0; j < *J; j++) {
        for (h = 0; h < *H; h++) {
          tail_index = g + *G * h;
          u = g + *G * (j + *J * h);
          qf_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                 kappa, X, y, y_low, y_high, theta, &shape[h],
                 q_ptr[h], p_ptr[h], log_d_ptr[h],
                 AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
                 &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                 &g, &j, &h);
          
          clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                    kappa, X, y, y_low, y_high, theta, &shape[h],
                    q_ptr[h], p_ptr[h], log_d_ptr[h],
                    AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
                    &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                    w, M_low, M_high, I_low, I_high,
                    &g, &j, &h);
          
          can_ll_sum[u] = ll_sum[u];
          for (n = 0; n < *N; n++) {
            can_log_like[u + GJH * n] = log_like[u + GJH * n];
          }
        }
      }
      if (*cop > 0) {
        copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]);
        can_cop_ll_sum[g] = cop_ll_sum[g];
        for (n = 0; n < *N; n++) {
          can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
        }
      }
    }
  } /*  adjust candidate standard deviations */
    if ( ((i+1) % 100 == 0) && (i < burn) ) {
      for (k = 0; k < dim_thetastar; k++) {
        if (acc_thetastar[k]  > (0.5 * att_thetastar[k])) {
          tuning_theta[k] *= 1.2;
        } else if (acc_thetastar[k]  < (0.3 * att_thetastar[k])) {
          tuning_theta[k] *= 0.8;
        }
        acc_thetastar[k] = 0;
        att_thetastar[k] = 0;
      }
      for (h = 0; h < *H; h++) {
        if (acc_shape[h] > 50) {
          tuning_shape[h] *= 1.2;
        } else if (acc_shape[h] < 30) {
          tuning_shape[h] *= 0.8;
        }
        acc_shape[h] = 0;
      }
      for (k = 0; k < GH; k++) {
        if (acc_xi_low[k] > 50) {
          tuning_xi_low[k] *= 1.2;
        } else if (acc_xi_low[k] < 30) {
          tuning_xi_low[k] *= 0.8;
        }
        acc_xi_low[k] = 0;
        if (acc_xi_high[k] > 50) {
          tuning_xi_high[k] *= 1.2;
        } else if (acc_xi_high[k] < 30) {
          tuning_xi_high[k] *= 0.8;
        }
        acc_xi_high[k] = 0;
      }
      for (g = 0; g < *G; g++) {
        if (acc_lambda[g]  > 50) {
          tuning_nu[g] *= 0.8;
        } else if (acc_lambda[g]  < 30) {
          tuning_nu[g] *= 1.2;
        }
        acc_lambda[g] = 0;
        if (tuning_nu[g] < *H + 2.0) {
          tuning_nu[g] = *H + 2.0;
        }
        for (h = 0; h < *H * *Q_int; h++) {
          k = g + *G * h;
          if (acc_delta[k]  > 50) {
            tuning_delta[k] *= 1.2;
          } else if (acc_delta[k]  < 30) {
            tuning_delta[k] *= 0.8;
          }
          acc_delta[k] = 0;
        }
      }
    }
  //reset acceptance probabilities
  if (i == burn) {
    m_init = 1;
    
    for (k = 0; k < dim_thetastar; k++) {
      acc_thetastar[k] = 0;
      att_thetastar[k] = 0;
    }
    for (h = 0; h < *H; h++) {
      acc_shape[h] = 0;
    }
    for (k = 0; k < *G * *H; k++) {
      acc_xi_low[k] = 0;
      acc_xi_high[k] = 0;
    }
    for (g = 0; g < *G; g++) {
      acc_lambda[g] = 0;
    }
    for (k = 0; k < dim_delta; k++) {
      acc_delta[k] = 0;
    }
  }
  /*** store variables ***/
    if ((i >= keep)) {
      for (h = 0; h < *H; h++) {
        for (g = 0; g < *G; g++) {
          for (m = 0; m < *M; m++) {
            for (p = 0; p < *P; p++) {
              for (j = 0; j < *J; j++) {
                k = g + *G * j + *G * *J * h + *G * *J * *H * (m + *M * p);
                theta_keep[k] = sd_H[h] * theta[k];
              }
              for (d = 0; d < *D; d++) {
                k = g + *G * d + *G * *D * h + *G * *D * *H * (m + *M * p);
                thetastar_keep[k] = sd_H[h] * thetastar[k];
              }
            }
          }
          for (j = 0; j < *J; j++) {
            theta_keep[g + *G * j + *G * *J * h + *G * *J * *H * (0 + *M * 0)] += mn_H[h];
          }
          thetastar_keep[g + *G * 0 + *G * *D * h + *G * *D * *H * (0 + *M * 0)] += mn_H[h];
        }
      }
      if((i + 1 - burn) % thin == 0) {
        i2 = (i + 1 - burn) / thin - 1;
        PARETO[i2 + num_keep * 0] = pareto_low;
        PARETO[i2 + num_keep * 1] = pareto_high;
        XI_RHO[i2 + num_keep * 0] = xi_low_rho;
        XI_RHO[i2 + num_keep * 1] = xi_high_rho;
        for (h = 0; h < *H; h++) {
          XI_MU[i2 + num_keep * (h + *H * 0)]  = xi_low_mu[h];
          XI_MU[i2 + num_keep * (h + *H * 1)]  = xi_high_mu[h];
          SHAPE[i2 + num_keep * h] = shape[h];
          for (v = 0; v < *H; v++) {
            XI_SIG[i2 + num_keep * (h + *H * v) + num_keep * H2 * 0] = xi_low_sig_h[(h + *H * v)];
            XI_SIG[i2 + num_keep * (h + *H * v) + num_keep * H2 * 1] = xi_high_sig_h[(h + *H * v)];
          }
          MakeBasis(&basis_ind[h], M, kappa, N_tau_q, tau_q, &shape[h], B_q);
          for (g = 0; g < *G; g++) {
            l = g + *G * h;
            XI[i2 + num_keep * l + num_keep * *G * *H * 0]    =   xi_low[l];
            XI[i2 + num_keep * l + num_keep * *G * *H * 1]    =   xi_high[l];
            for (j = 0; j < *J; j++) {
              u = g + *G * j + *G * *J * h;
              for (p = 0; p < *P; p++) {
                for (n = 0; n < *N_tau_q; n++) {
                  k = g + *G * j + *G * *J * h + GJH * p;
                  if (tau_q[n] < *tau_low) {
                    if (pareto_low == 1) {
                      dummy_log_det = pow(tau_q[n] / *tau_low, - exp(xi_low[l])) - 1;
                      dummy_log_det /= exp(xi_low[l]);
                    } else{
                      dummy_log_det = log(tau_q[n] / *tau_low);
                    }
                    post_q[i2 + num_keep * (n + *N_tau_q * k)] = 1.0 * theta_keep[u + GJH * (0 + *M * p)];
                    for (m = 0; m < M_1; m++) {
                      post_q[i2 + num_keep * (n + *N_tau_q * k)] += (I_low[m + 1] + M_low[m + 1] * dummy_log_det) * theta_keep[u + GJH * ((m + 1) + *M * p)];
                    }
                  } else if (tau_q[n] > *tau_high) {
                    if (pareto_high == 1) {
                      dummy_log_det = pow((1 - tau_q[n]) / (1 - *tau_high), - exp(xi_high[l])) - 1;
                      dummy_log_det /= exp(xi_high[l]);
                    } else{
                      dummy_log_det = log((1 - tau_q[n]) /(1 - *tau_low));
                    }
                    post_q[i2 + num_keep * (n + *N_tau_q * k)] = 1.0 * theta_keep[u + GJH * (0 + *M * p)];
                    for (m = 0; m < M_1; m++) {
                      post_q[i2 + num_keep * (n + *N_tau_q * k)] += (I_high[m + 1] + M_high[m + 1] * dummy_log_det) * theta_keep[u + GJH * ((m + 1) + *M * p)];
                    }
                  } else{
                    for (m = 0; m < *M; m++) {
                      post_q[i2 + num_keep * (n + *N_tau_q * k)] += B_q[n + *N_tau_q * m] * theta_keep[u + GJH * (m + *M * p)];
                    }
                  }
                }
              }
            }
          }
        }
        if (*cop > 0) {
          for (g = 0; g < *G; g++) {
            for (n = 0; n < n_g[g]; n++) {
              for (h = 0; h < HQ; h++) {
                k = g + *G * n + *G * *N * h;
                GAMMA[i2 + num_keep * k] = gamma[k];
              }
              for (j = 0; j < JH; j++) {
                k = g + *G * j + *G * JH * n;
                ETA[i2 + num_keep * k] = eta[k];
              }
            }
            ALPHA[i2 + num_keep * g] = alpha[g];
            for (h = 0; h < H2; h++) {
              LAMBDA[i2 + num_keep * (g + *G * h)] = lambda[(g + *G * h)];
            }
            for (h = 0; h < HQ; h++) {
              DELTA[i2 + num_keep * (g + *G * h)] = delta[(g + *G * h)];
            }
          }
        }
        for (d = 0; d < *D; d++) {
          for (m = 0; m < *M; m++) {
            for (p = 0; p < *P; p++) {
              prior_index = d + *D * m + *D * *M * p;
              RHO[i2 + num_keep * prior_index] = rho[prior_index];
              for (h = 0; h < *H; h++) {
                MU[i2 + num_keep * (prior_index + DMP * h)] = mu[prior_index + DMP * h];
                for (v = 0; v < *H; v++) {
                  k = prior_index + DMP * (h + *H * v);
                  SIGMA[i2 + num_keep * k] = sig_h[k];
                }
                for (g = 0; g < *G; g++) {
                  k = g + *G * d + *G * *D  * h + *G * *D * *H * (m + *M * p);
                  THETASTAR[i2 + num_keep * k] = thetastar_keep[k];
                }
              }
            }
          }
        }
      }
      /*update model information criteria parameters irrespective of thinning*/
        for(h = 0; h < *H; h++) {
          meanshape[h] += shape[h];
          for (g = 0; g < *G; g++) {
            l = g + *G * h;
            meanxilow[l] += xi_low[l];
            meanxihigh[l] += xi_high[l];
          }
        }
      for (k = 0; k < dim_theta; k++) {
        meantheta[k] += theta[k];
      }
      if (*cop > 0) {
        for (g = 0; g < *G; g++) {
          for (n = 0; n < n_g[g]; n++) {
            for (h = 0; h < HQ; h++) {
              k = g + *G * n + *G * *N * h;
              meangamma[k] += gamma[k];
            }
            for (j = 0; j < JH; j++) {
              k = g + *G * j + *G * JH * n;
              meaneta[k] += eta[k];
            }
          }
          meanalpha[g] += alpha[g];
          for (h = 0; h < H2; h++) {
            meanlambda[g + *G * h] += lambda[g + *G * h];
          }
          for (h = 0; h < HQ; h++) {
            meandelta[g + *G * h] += delta[g + *G * h];
          }
        }
      }
      
      k = (i  + 1 - burn);
      if (k > 1) {
        for (n = 0; n < dim_Y; n++) {
          Y_var[n] = ((k - 2) * Y_var[n] / (k - 1) + pow(y[n] - Y_mean[n],2) / k) ;
          U_var[n] = ((k - 2) * U_var[n] / (k - 1) + pow(qnorm_U[n] - U_mean[n],2) / k) ;
        }
      }
      for (n = 0; n < dim_Y; n++) {
        Y_mean[n] = (Y_mean[n] * (k - 1) + y[n]) / k;
        qnorm_U[n] = qnorm(U[n], 0, 1, 1, 0);
        U_mean[n] = (U_mean[n] * (k - 1) +     qnorm_U[n]) / k;
      }
      if (*cop == 0) {
        for (n = 0; n < dim_Y; n++) {
          CPO[n]        += exp(-log_like[n]);
          lppd_vec[n]   += exp(log_like[n]);
          if(k > 1) {
            log_dens_var[n] = ((k - 2) * log_dens_var[n] / (k - 1) + pow(log_like[n] - log_dens_mean[n], 2) / k);
          }
          log_dens_mean[n] = (log_dens_mean[n] * (k - 1) + log_like[n]) / k;
          epostlogpy    += log_like[n];
        }
      }
      if (*cop > 0) {
        for (g = 0; g < *G; g++) {
          for (n = 0; n < *N; n++) {
            v = g + *G * n;
            /*dummy_log_det is log density of kth observation*/
              this_dens = copula_loglike[v];
              for (j = 0; j < *J; j++) {
                for (h = 0; h < *H; h++) {
                  u = g + *G * j + GJ * (h + *H * n);
                  this_dens += log_like[u];
                }
              }
              CPO[v]      += exp(-this_dens);
              lppd_vec[v] += exp(this_dens);
              epostlogpy += this_dens;
              if(k > 1) {
                log_dens_var[v] = ((k - 2) * log_dens_var[v] / (k - 1) + pow(this_dens - log_dens_mean[v], 2) / k);
              }
              log_dens_mean[v] = (log_dens_mean[v] * (k - 1) + this_dens) / k;
          }
        }
      }
    }
}
PutRNGstate();

if (*cop > 0) {
  for (g = 0; g < *G; g++) {
    for (h = 0; h < *H; h++) {
      for (n = 0; n < *N; n++) {
        for (j = 0; j < *J; j++) {
          W[g + *G * j + *G * *J * h + *G * *J * *H * n] = W2[g + *G * (j * *H + h + JH * n)];
        }
      }
    }
    copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]);
    for (n = 0; n < *N; n++) {
      can_copula_loglike[g + *G * n] = copula_loglike[g + *G * n];
    }
    can_cop_ll_sum[g] = cop_ll_sum[g];
  }
}

epostlogpy /= (double) (keep);

for (n = 0; n < dim_CPO; n++) {
  CPO[n] = (double) (keep) / CPO[n];
  LPML += log(CPO[n]);
  lppd_vec[n] /= (double) (keep);
  lppd_vec[n] = log(lppd_vec[n]);
  lppd_sum += lppd_vec[n];
  pwaic2 += log_dens_var[n];
}
elppd = lppd_sum - pwaic2;
WAIC = -2 * elppd;

/*calculate means of parameters for deviance*/
  for (h = 0; h < *H; h++) {
    meanshape[h] /= (double) (keep);
  }
for (k = 0; k < dim_theta; k++) {
  meantheta[k] /= (double) (keep);
}
for (k = 0; k < dim_xi; k++) {
  meanxilow[k] /= (double) (keep);
  meanxihigh[k] /= (double) (keep);
}
for (k = 0; k < dim_gamma; k++) {
  meangamma[k] /= (double) (keep);
}
for (k = 0; k < dim_eta; k++) {
  meaneta[k] /= (double) (keep);
}
for (g = 0; g < *G; g++) {
  meanalpha[g] /= (double) (keep);
}
for (k = 0; k < dim_lambda; k++) {
  meanlambda[k] /= (double) (keep);
}
for (k = 0; k < dim_delta; k++) {
  meandelta[k] /= (double) (keep);
}
/* need to calculate W_cov again at mean values */
  if (*cop > 1) {
    for (g = 0; g < *G; g++) {
      for (u = 0; u < *J; u++) {
        for (v = 0; v < *J; v++) {
          k = u + *J * v;
          cov_J[k] = pow(meanalpha[g], dist_J[u + *J * v]);
        }
      }
      CholMatrix(J, cov_J, prec_J, &dummy_log_det);
      for (u = 0; u < *J; u++) {
        for (v = 0; v < *J; v++) {
          k = u + *J * v;
          can_AR_chol[g + *G * k] = prec_J[k];
        }
      }
      InvertSymmetricMatrix(J, cov_J, prec_J, &can_alpha_log_det[g]);
      for (u = 0; u < *J; u++) {
        for (v = 0; v < *J; v++) {
          k = u + *J * v;
          can_AR[g + *G * k]      = cov_J[k];
          can_AR_prec[g + *G * k] = prec_J[k];
        }
      }
      for (h = 0; h < H2; h++) {
        cov_H[h] = meanlambda[g + *G * h];
      }
      CholMatrix(H, cov_H, prec_H, &dummy_log_det);
      for (h = 0; h < H2; h++) {
        can_lambda_chol[g + *G * h] = prec_H[u + *H * v];
      }
      InvertSymmetricMatrix(H, cov_H, prec_H, &can_lambda_log_det[g]);
      Kronecker(J, J, H, H, cov_J, cov_H, cov_JH);
      Kronecker(J, J, H, H, prec_J, prec_H, prec_JH);
      for (u = 0; u < JH; u++) {
        for (v = 0; v < JH; v++) {
          k = u + JH * v;
          can_eta_cov[g + *G * k] = cov_JH[k];
          can_eta_prec[g + *G * k] = prec_JH[k];
          can_W_cov[g + *G * k] = can_eta_cov[g + *G * k];
        }
      }
      for (u = 0; u < JH; u++) {
        can_W_cov[g + *G * (u + JH * u)] ++;
        cov_JH[u + JH * u] ++;
        prec_JH[u + JH * u] ++;
      }
      make_D(G, J, H, Q_int, n_g, N, Z, meandelta, can_W_cov, &g, status, can_DD_inv, can_DD);
    }
  }
/*calculate loglikelihood at posterior mean*/
  for (g = 0; g < *G; g++) {
    for (h = 0; h < *H; h++) {
      k = g + *G * h;
      tail_index = g + *G * h;
      make_b(q_ptr[h], M, kappa, &meanshape[h], B);
      for (j = 0; j < *J; j++) {
        u = g + *G * j + *G * *J * h;
        clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                  kappa, X, y, y_low, y_high, meantheta, &meanshape[h],
                  q_ptr[h], p_ptr[h], log_d_ptr[h],
                  AAAA, B, bin, bin_low, bin_high, &ll_sum[u], can_log_like, can_U, U_low, U_high,
                  &pareto_low, &pareto_high, tau_low, tau_high, &meanxilow[tail_index], &meanxihigh[tail_index],
                  w, M_low, M_high, I_low, I_high,
                  &g, &j, &h);
        llpm += ll_sum[u];
      }
    }
  }
if (*cop > 0) {
  for (g = 0; g < *G; g++) {
    for (h = 0; h < *H; h++) {
      for (n = 0; n < *N; n++) {
        for (j = 0; j < *J; j++) {
          q_int = g + *G * j + *G * *J * h + GJH * n;
          can_W[q_int] = qnorm(can_U[q_int], 0, 1, 1, 0);
          can_W2[g + *G * (j * *H + h + JH * n)] = can_W[q_int];
        }
      }
    }
    copula_ll(G, J, H, Q_int, n_g, N, &g, can_W2, Z, meaneta, can_DD, can_DD_inv, meangamma, can_copula_loglike, &can_cop_ll_sum[g]);
    llpm +=  can_cop_ll_sum[g];
  }
}
pdic = 2.0 * (llpm - epostlogpy);
DIC = -2.0 * llpm + 2.0 * pdic;

MIC[0] = LPML;
MIC[1] = DIC;
MIC[2] = WAIC;

if (verbose) {
  for (k = 0; k < DMP; k++) {
    Rprintf("log_prior final = %f\n", lp[k]);
  }
}

for (g = 0; g < *G; g++) {
  for (h = 0; h < *H; h++) {
    k = g + *G * h;
    tail_index = g + *G * h;
    make_b(q_ptr[h], M, kappa, &shape[h], B);
    for (j = 0; j < *J; j++) {
      u = g + *G * j + *G * *J * h;
      clike_ptr(G, J, H, M, P, N, &n_gjh[u], status,
                kappa, X, y, y_low, y_high, theta, &shape[h],
                q_ptr[h], p_ptr[h], log_d_ptr[h],
                AAAA, B, bin, bin_low, bin_high, &ll_sum[u], log_like, U, U_low, U_high,
                &pareto_low, &pareto_high, tau_low, tau_high, &xi_low[tail_index], &xi_high[tail_index],
                w, M_low, M_high, I_low, I_high,
                &g, &j, &h);
      if (verbose) {
        Rprintf("loglike final = %f\n", ll_sum[u]);
      }
    }
  }
}
if (*cop > 0) {
  for (g = 0; g < *G; g++) {
    copula_ll(G, J, H, Q_int, n_g, N, &g, W2, Z, eta, DD, DD_inv, gamma, copula_loglike, &cop_ll_sum[g]); if (verbose) {Rprintf("final cop_ll_sum = %f\n", cop_ll_sum[g]);}
  }
}
for (g = 0; g < *G; g++) {
  for (j = 0; j < *J; j++) {
    for (h = 0; h < *H; h++) {
      for (i = 0; i < *N; i++) {
        k = g + *G * j + *G * *J * h + GJH * i;
        Y_mean[k] = Y_mean[k] * sd_H[h] + mn_H[h];
        Y_var[k]  = Y_var[k]  * sd_H[h] * sd_H[h];
      }
    }
  }
}
for (i = 0; i < dim_Y; i++) {
  Y_array[3 + 7 * i]   = Y_mean[i];
  Y_array[4 + 7 * i]   = Y_var[i];
  Y_array[5 + 7 * i]   = U_mean[i];
  Y_array[6 + 7 * i]   = U_var[i];
}

/***Free memory */
  /*data and regression variables */
  Free(n_g);
Free(status_2);

Free(ll_sum);
Free(can_ll_sum);

Free(log_like);
Free(can_log_like);
Free(CPO);
Free(lppd_vec);
Free(log_dens_mean);
Free(log_dens_var);

Free(y);
Free(y_low);
Free(y_high);
Free(Y_mean);
Free(Y_var);
Free(U_mean);
Free(U_var);

Free(U);
Free(can_U);
Free(U_low);
Free(U_high);
Free(qnorm_U);
Free(bin);
Free(bin_low);
Free(bin_high);

Free(cop_ll_sum);
Free(can_cop_ll_sum);
Free(DD);
Free(DD_inv);
Free(can_DD);
Free(can_DD_inv);

Free(missing_index);
Free(dead_index);
Free(dead_sum);
Free(start_index);

Free(theta);
Free(cantheta);
Free(canthetastar);
Free(tuning_theta);
Free(acc_thetastar);
Free(att_thetastar);
Free(meantheta);
Free(X_D);
Free(dummy_seq);

Free(mn_H);
Free(sd_H);

Free(meanthetastar);
Free(varthetastar);
Free(corrthetastar);
Free(z);

/*basis matrix variables*/
  Free(B);
Free(can_B);
Free(B_q);

Free(AAAA);
/* prior variables */
  Free(mu);
Free(sig_h);
Free(rho);
Free(sig_a);
Free(sig_b);

Free(omega);
Free(canomega);
Free(prec_g);
Free(can_prec_g);
Free(dummy_prec_g);
Free(prec_h);
Free(dummy_sig_h);
Free(dummy_prec_h);

Free(lp);
Free(rhoquad);
Free(canrhoquad_array);
Free(sumalpha);
Free(log_det_g);
Free(log_det_h);
Free(mu_array);
Free(dummy_thetastar_gh);
Free(g1_vec);
Free(g2_vec);


Free(col_sums);
Free(dummy_mean_matrix);
Free(dummy_mean_mu);
Free(mean_mu);
Free(chol_mu);
Free(dummy_chol_mu);
Free(z_mu);
Free(sai_h);
Free(W_h);

/*shape variables */
  Free(shape);
Free(shape_mean);
Free(shape_var);
Free(acc_shape);
Free(tuning_shape);
Free(shapeprec);
Free(meanshape);

/*copula variables*/
  Free(W);
Free(can_W);
Free(W2);
Free(can_W2);
Free(copula_loglike);
Free(can_copula_loglike);


Free(gamma);
Free(omega_gamma);
Free(Omega_gamma);
Free(Omega_gamma_inv);
Free(Omega_gamma_inv_chol);

Free(dummy_vec_HQ);
Free(z_HQ);
Free(Z2_sum);
Free(meangamma);

Free(eta);
Free(omega_eta);
Free(Omega_eta_inv);
Free(Omega_eta_inv_chol);

Free(dummy_vec_JH);
Free(z_JH);

Free(eta_cov);
Free(eta_prec);
Free(W_cov);

Free(can_eta_cov);
Free(can_eta_prec);
Free(can_W_cov);

Free(cov_JH);
Free(prec_JH);

Free(meaneta);

Free(alpha);
Free(can_alpha);
Free(AR);
Free(AR_prec);
Free(AR_chol);
Free(alpha_log_det);

Free(can_AR);
Free(can_AR_prec);
Free(can_AR_chol);
Free(can_alpha_log_det);

Free(cov_J);
Free(prec_J);
Free(meanalpha);
Free(dummy_mat_J2);
Free(alpha_det_matrix);

Free(lambda);
Free(lambda_prec);
Free(lambda_chol);
Free(lambda_log_det);

Free(can_lambda);
Free(can_lambda_prec);
Free(can_lambda_chol);
Free(can_lambda_log_det);

Free(tuning_nu);
Free(acc_lambda);
Free(cov_H);
Free(prec_H);
Free(meanlambda);

Free(delta);
Free(delta_prec);
Free(can_delta);
Free(can_delta_prec);

Free(tuning_delta);
Free(acc_delta);
Free(zeta_A);
Free(zeta_B);

Free(meandelta);

/* tail variables */
  Free(xi_low);
Free(xi_high);

Free(can_xi_low);
Free(can_xi_high);

Free(xi_low_mu);
Free(xi_high_mu);

Free(xi_low_mu_array);
Free(xi_high_mu_array);
Free(xi_low_sai_0);
Free(xi_high_sai_0);

Free(xi_low_prec_g);
Free(xi_low_sig_h);
Free(xi_low_prec_h);
Free(xi_low_omega);
Free(can_xi_low_omega);

Free(xi_high_prec_g);
Free(xi_high_sig_h);
Free(xi_high_prec_h);
Free(xi_high_omega);
Free(can_xi_high_omega);

Free(meanxilow);
Free(tuning_xi_low);
Free(acc_xi_low);

Free(meanxihigh);
Free(tuning_xi_high);
Free(acc_xi_high);

Free(w);
Free(I_low);
Free(I_high);
Free(M_low);
Free(M_high);
      }
