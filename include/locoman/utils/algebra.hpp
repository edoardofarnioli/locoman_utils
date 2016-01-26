#ifndef _ALGEBRA_HPP_
#define _ALGEBRA_HPP_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>
#include <yarp/math/SVD.h>
#include <unistd.h>


#include <locoman/utils/screws.hpp>
#include <locoman/utils/kinematics.hpp>
#include <locoman/utils/kinetostatics.hpp>
#include <locoman/utils/locoman_utils.hpp>
//#include <locoman/utils/algebra.hpp>

using namespace yarp::math;

namespace locoman {
        namespace utils {
                 //----------------------------------------------------------------------------
     /**
     * @brief  RoundMatrix computes the round of M to k decimal places 
     * @param  M yarp matrix to round
     * @param  k number of decimal places 
     * @return rounded yarp matrix 
     */
    yarp::sig::Matrix RoundMatrix( const yarp::sig::Matrix& M, const int k) {
          yarp::sig::Matrix Round_M = M ;
  for( int i=0 ; i<Round_M.cols(); i++) {
      for( int j=0 ; j<Round_M.rows(); j++){
      Round_M[j][i]  = Round_M[j][i]  * std::pow(10,k) ;
      Round_M[j][i]  = floor( Round_M[j][i] +0.5) ;
      Round_M[j][i]  = Round_M[j][i]* std::pow(10,-k) ;
      } ;
  } ;
return Round_M ;
    }
    
    
    //------------------------------------------------------------------------------------
     /**
     * @brief  Pinv_trunc_SVD computes the pseudoinverse of A via the truncated SVD method 
     * @param  A is the matrix to pseudo-inverse
     * @param  k is the maximum ratio admitted between the max and min singular values
     * @return Pinv_trunc_SVD is the pseudo-inverse of A
     */
      yarp::sig::Matrix Pinv_trunc_SVD( const yarp::sig::Matrix& A ,
                                        const double k = 1E-4 
                                      ) {
  int r_A = A.rows() ;
  int c_A = A.cols() ;
  int c_U = c_A;
  if(c_U>r_A ) c_U = r_A ;
  yarp::sig::Matrix U( r_A, c_U )  ;
  yarp::sig::Matrix V( c_A, c_U )  ;
  yarp::sig::Vector S( c_U ) ;  
  yarp::sig::Matrix S_1( c_U , c_U ) ;  
  S_1.zero();
  
  yarp::math::SVD(A, U, S, V );
  
  double k_norm = S(0)*k ;
        for(int i = 0 ;  i < c_U ; i++)
    {
      if( S(i)<k_norm) 
      {
        S_1(i,i)=0.0; 
      }
      else
      {
        S_1(i,i) = 1/S(i) ;
      }
    }  
    yarp::sig::Matrix pinv_A = V * S_1 *U.transposed() ;
  return pinv_A ;
                                    }   // 
//------------------------------------------------------------------------------------
     /**
     * @brief  Pinv_Regularized computes the Levemberg Regularized pseudo-inverse of A 
     * @param  A is the matrix to pseudo-inverse
     * @param  k is the regularization factor
     * @return Pinv_Regularized is the pseudo-inverse of A
     */ 
      yarp::sig::Matrix Pinv_Regularized( const yarp::sig::Matrix& A ,
                                          const double k  
                                        ) {
                                            int r_A = A.rows() ;
                                            int c_A = A.cols() ;
                                            if(c_A< r_A ){    // tall matrix => Defective manipulator
                                            yarp::sig::Matrix At_A = A.transposed()*A ;
                                            yarp::sig::Matrix Eye_temp_l = yarp::math::eye(At_A.rows(), At_A.cols() ) ;
                                            yarp::sig::Matrix pinv_A_l = locoman::utils::Pinv_trunc_SVD(At_A + k*Eye_temp_l , 1E-7 ) * A.transposed() ;
                                            return pinv_A_l ;
                                            }
                                            else {    // fat (or square) matrix => Redundant (or square) manipulator
                                            yarp::sig::Matrix A_At = A*A.transposed() ;
                                            yarp::sig::Matrix Eye_temp_r = yarp::math::eye( A_At.rows(), A_At.cols() ) ;
                                            yarp::sig::Matrix pinv_A_r = A.transposed()* locoman::utils::Pinv_trunc_SVD(A_At + k*Eye_temp_r , 1E-7 )  ;
                                            return pinv_A_r ;
                                            }
                                        }     
                                        
     /**
     * @brief  Pinv_Marq computes the Levemberg-Marquard Regularized pseudo-inverse of A 
     * @param  A is the matrix to pseudo-inverse
     * @param  k is the regularization factor
     * @return Pinv_Marq is the pseudo-inverse of A
     */ 
      yarp::sig::Matrix Pinv_Marq( const yarp::sig::Matrix& A ,
                                          const double k  
                                        ) {
                                    int r_A = A.rows() ;
                                    int c_A = A.cols() ;
                                    if(c_A< r_A ){    // tall matrix => Defective manipulator    
                                    yarp::sig::Matrix At_A = A.transposed()*A ;
                                    int dim_l = At_A.rows() ;
                                    yarp::sig::Matrix Marq_temp_l(At_A.rows(), At_A.cols() ) ;
                                    Marq_temp_l.eye(); 
                                    for(int i=0; i < dim_l ; i++){
                                        Marq_temp_l[i][i] = At_A[i][i] ;   
                                    }   
                                    yarp::sig::Matrix pinv_A_l = Pinv_trunc_SVD( At_A + k*Marq_temp_l,  1E-7   ) *A.transposed() ;
                                    return pinv_A_l ;      
                                    }
                                    else {    // fat (or square) matrix => Redundant (or square) manipulator
                                    yarp::sig::Matrix A_At = A*A.transposed() ;
                                    int dim_r = A_At.rows() ; 
                                    yarp::sig::Matrix Marq_temp_r(A_At.rows(), A_At.cols() ) ;
                                    Marq_temp_r.eye(); 
                                    for(int i=0; i < dim_r ; i++){
                                        Marq_temp_r[i][i] = A_At[i][i] ;   
                                    }   
                                    yarp::sig::Matrix pinv_A_r = A.transposed()* Pinv_trunc_SVD(A_At + k*Marq_temp_r , 1E-7 )  ;
                                    return pinv_A_r ;
                                    }
                                        }    
//------------------------------------------------------------------------------------
     /**
     * @brief  x_Pinv_Iter computes the variable x: Ax=b via the Landweber iteration method
     * @param  A is the matrix to pseudo-inverse
     * @param  b is the vector of known terms
     * @param  n is the maximum number of steps to be performed (less that the minimum dimension of A)
     * @return x_Pinv_Iter is the solution vetor
     */
      yarp::sig::Vector x_Pinv_Iter( const yarp::sig::Matrix& A , 
                                   const yarp::sig::Vector& b , 
                                   double n 
                                   ) {
                                    int r_A = A.rows() ;
                                    int c_A = A.cols() ;
                                    if(n>r_A ) n = r_A ;
                                    if(n>c_A ) n = c_A ;
                                    yarp::sig::Matrix At_A = A.transposed()*A ;
                                    yarp::sig::Matrix Eye_temp = yarp::math::eye(At_A.rows(), At_A.cols() ) ;
                                    yarp::sig::Vector At_b = A.transposed()*b ;
                                    yarp::sig::Vector x_k(c_A , 0.0) ;
                                    yarp::sig::Vector x_k1(c_A , 0.0) ;

                                    for(int i = 0 ;  i < n ; i++)
                                    {
                                    x_k1 = x_k + A.transposed()*b - A.transposed()*A*x_k ;
                                    //  x_k1 = (Eye_temp - At_A) * x_k + At_b ;
                                    x_k = x_k1 ;
                                    } ;
                                    return x_k1 ;                             
          
                                }   //                           
     /**
     * @brief  orth_SVD computes a basis for the span of A via the truncated SVD method 
     * @param  A is the matrix of which a basis is needed
     * @param  k is the maximum ratio admitted between the max and min singular values
     * @return orth_SVD is a basis for the span of A
     */
      yarp::sig::Matrix orth_SVD( const yarp::sig::Matrix& A ,
                                        const double k = 1E-4 
                                      ) {
                                            int r_A = A.rows() ;
                                            int c_A = A.cols() ;
                                            int c_U = c_A;
                                            if(c_U>r_A ) c_U = r_A ;
                                            yarp::sig::Matrix U( r_A, c_U )  ;
                                            yarp::sig::Matrix V( c_A, c_U )  ;
                                            yarp::sig::Vector S( c_U ) ;  
                                            yarp::sig::Matrix S_1( c_U , c_U ) ;  
                                            S_1.zero();
                                            yarp::math::SVD(A, U, S, V );
                                            int cont = 0;
                                            double k_norm = S(0)*k ;
                                                    for(int i = 0 ;  i < c_U ; i++)
                                                {
                                                if( S(i)<k_norm) 
                                                {
                                                    S_1(i,i)=0.0; 
                                                }
                                                else
                                                {//span_A.resize();
                                                //bool a = span_A.setCol(i, U.getCol(i) ) ;//S_1(i,i) = 1/S(i) ;
                                                cont += 1 ;
                                                }
                                                }  
                                                yarp::sig::Matrix span_A = U.submatrix(0, U.rows()-1, 0, (cont-1) );
                                            return span_A ;
                                    }   //                             
     /**
     * @brief  null_SVD computes a basis for the nullspace of A via the truncated SVD method 
     * @param  A is the matrix of which a basis for the nullspace is needed
     * @param  k is the maximum ratio admitted between the max and min singular values
     * @return null_SVD is a basis for the nullspace of A
     */
      yarp::sig::Matrix null_SVD( const yarp::sig::Matrix& A ,
                                        const double k = 1E-4 
                                      ) {
                                            int r_A = A.rows() ;
                                            int c_A = A.cols() ;
                                            int c_U = c_A;
                                            if(c_U>r_A ) c_U = r_A ;
                                            yarp::sig::Matrix U( r_A, c_U )  ;
                                            yarp::sig::Matrix V( c_A, c_U )  ;
                                            yarp::sig::Vector S( c_U ) ;  
                                            yarp::sig::Matrix S_1( c_U , c_U ) ;  
                                            S_1.zero();
                                            yarp::math::SVD(A, U, S, V );
                                            double k_norm = S(0)*k ;
                                                    for(int i = 0 ;  i < c_U ; i++)
                                                {
                                                if( S(i)<k_norm) 
                                                {
                                                    S_1(i,i)=0.0; 
                                                }
                                                else
                                                {
                                                    S_1(i,i) = S(i) ;
                                                }
                                                }  
                                            yarp::sig::Matrix A_filter = U * S_1 * V.transposed() ;
                                            yarp::sig::Matrix null_project_A = yarp::math::nullspaceProjection(A_filter) ;
                                            yarp::sig::Matrix null_A = orth_SVD(null_project_A,k) ;
                                            return null_A ;                                
                                    }  
           
     /**
     * @brief  filter_SVD computes a "filtered" version of A via the truncated SVD  
     * @param  A is the matrix of which the filtered version is needed
     * @param  k is the maximum ratio admitted between the max and min singular values
     * @return filter_SVD is the filtered version of A
     */
      yarp::sig::Matrix filter_SVD( const yarp::sig::Matrix& A ,
                                        const double k = 1E-4 
                                      ) {
                                            int r_A = A.rows() ;
                                            int c_A = A.cols() ;
                                            int c_U = c_A;
                                            if(c_U>r_A ) c_U = r_A ;
                                            yarp::sig::Matrix U( r_A, c_U )  ;
                                            yarp::sig::Matrix V( c_A, c_U )  ;
                                            yarp::sig::Vector S( c_U ) ;  
                                            yarp::sig::Matrix S_1( c_U , c_U ) ;  
                                            S_1.zero();
                                            yarp::math::SVD(A, U, S, V );
                                            double k_norm = S(0)*k ;
                                                    for(int i = 0 ;  i < c_U ; i++)
                                                {
                                                if( S(i)<k_norm) 
                                                {
                                                    S_1(i,i)=0.0; 
                                                }
                                                else
                                                {
                                                    S_1(i,i) = S(i) ;
                                                }
                                                }  
                                            yarp::sig::Matrix A_filter = U * S_1 * V.transposed() ;
                                            return A_filter ;                                
                                    }
        }
}
#endif