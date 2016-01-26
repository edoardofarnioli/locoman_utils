#ifndef _KINETOSTATICS_HPP_
#define _KINETOSTATICS_HPP_

#include <cstdlib>
#include <iostream>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include <locoman/utils/screws.hpp>
#include <locoman/utils/kinematics.hpp>
// #include <locoman/utils/kinetostatics.hpp>
#include <locoman/utils/locoman_utils.hpp>
#include <locoman/utils/algebra.hpp>

namespace locoman {
        namespace utils {
            
            // pre declaration
            yarp::sig::Matrix D_Jacob_spa_i( const yarp::sig::Matrix& J_s, const int i );
            yarp::sig::Matrix Adjoint( const yarp::sig::Matrix& T_ab);
            yarp::sig::Matrix ad_lie( const yarp::sig::Vector& Xi);
                     //----------------------------------------------------------------------------
     /**
     * @brief  Q_ci compute the derivative of the spatial Jacobian 
     * @param  J_spa_i is a 6xc yarp matrix describing a spatial Jacobian
     * @param  T_a_ci is the homogeneous transformation between the floating base and the contact frame 
     * @param  f_ci contact force vector (3 force elements)
     * @return qxq yarp matrix 
     **/
    yarp::sig::Matrix Q_ci( const yarp::sig::Matrix& J_spa_i, 
                            const yarp::sig::Matrix& T_a_ci , 
                            const yarp::sig::Vector& f_ci) {
                                  yarp::sig::Matrix Q_c_i(J_spa_i.cols(), J_spa_i.cols());
                                    yarp::sig::Vector Q_ci_col_i( J_spa_i.cols() )  ; 
                                    yarp::sig::Matrix B(6,3) ;
                                    B[0][0] = 1 ;
                                    B[1][1] = 1 ;
                                    B[2][2] = 1 ;
                                        for ( int i = 0  ; i<(J_spa_i.cols()) ; i++ )  // i<(J_waist_l_c1_spa_0.cols()-1)
                                        {
                                        Q_ci_col_i =  ( locoman::utils::D_Jacob_spa_i(J_spa_i, (i+1)) ).transposed()*
                                                        ( locoman::utils::Adjoint(
                                                            locoman::utils::iHomogeneous(T_a_ci))).transposed()*B *f_ci  +
                                                        J_spa_i.transposed()* 
                                                        ( locoman::utils::Adjoint( locoman::utils::iHomogeneous(T_a_ci) )*
                                                        locoman::utils::ad_lie( -1.0*J_spa_i.getCol(i))).transposed()*  B *f_ci;
                                        Q_c_i.setCol(i, Q_ci_col_i ) ;       
                                        }   
                                    return Q_c_i ;
                            }
     /**
     * @brief  Q_ci_wrench compute the derivative of the spatial Jacobian 
     * @param  J_spa_i is a 6xc yarp matrix describing a spatial Jacobian
     * @param  T_a_ci is the homogeneous transformation between the floating base and the contact frame 
     * @param  w_ci contact force vector (6 elements: forces and couples)
     * @return qxq yarp matrix 
     **/
     yarp::sig::Matrix Q_ci_wrench( const yarp::sig::Matrix& J_spa_i, 
                                    const yarp::sig::Matrix& T_a_ci , 
                                    const yarp::sig::Vector& w_ci) {
                                    
                                    yarp::sig::Matrix Q_c_i(J_spa_i.cols(), J_spa_i.cols());
                                    yarp::sig::Vector Q_ci_col_i( J_spa_i.cols() )  ; 
//                                     yarp::sig::Matrix B(6,3) ;
//                                     B[0][0] = 1 ;
//                                     B[1][1] = 1 ;
//                                     B[2][2] = 1 ;
                                        for ( int i = 0  ; i<(J_spa_i.cols()) ; i++ )  // i<(J_waist_l_c1_spa_0.cols()-1)
                                        {
                                        Q_ci_col_i =  ( locoman::utils::D_Jacob_spa_i(J_spa_i, (i+1)) ).transposed()*
                                                        ( locoman::utils::Adjoint(
                                                            locoman::utils::iHomogeneous(T_a_ci))).transposed() *w_ci  +
                                                        J_spa_i.transposed()* 
                                                        ( locoman::utils::Adjoint( locoman::utils::iHomogeneous(T_a_ci) )*
                                                        locoman::utils::ad_lie( -1.0*J_spa_i.getCol(i))).transposed() *w_ci;
                                        Q_c_i.setCol(i, Q_ci_col_i ) ;       
                                        }   
                                    return Q_c_i ;
                            }
    
    
         //------------------------------------------------------------------------------------
     /**
     * @brief  FLMM_ext computes the FLMM for a compliant humanoid robot 
     * @param  J_c is a contacts x joints yarp matrix describing the body Jacobian of the humanoid
     * @param  S_c is a 6 x contacts yarp matrix describing the body Stance matrix of the humanoid
     * @param  Q_j is a joints x joints yarp matrix about the derivative of the Jacobian 
     * @param  Q_s is a contacts x joints yarp matrix about the derivative of the Jacobian 
     * @param  U_j is a joints x 6 yarp matrix about the derivative of the Jacobian 
     * @param  U_s is a 6 x 6 yarp matrix about the derivative of the Jacobian 
     * @param  K_c is a contacts x contacts yarp matrix describing the contact stiffness matrix
     * @param  K_q is a joints x joints yarp matrix describing the joint stiffness matrix
     * @return FLMM_ext is the Fundamental Loco-Manipulation Matrix
     */
    yarp::sig::Matrix FLMM_ext( const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_j,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_j,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c,
                                const yarp::sig::Matrix& K_q
                              ) {
                                    int size_fc  = J_c.rows();
                                    int size_q  = J_c.cols() ;
                                    int r_FLMM = size_fc + 6 + 2*size_q;
                                    int c_FLMM = size_fc + 6 + 3*size_q;
                                    yarp::sig::Matrix FLMM(r_FLMM, c_FLMM) ;    
                                    yarp::sig::Matrix Eye_fc(size_fc, size_fc) ;
                                    yarp::sig::Matrix Eye_q(size_q, size_q) ;
                                    Eye_fc.eye() ;
                                    Eye_q.eye() ;
                                    yarp::sig::Matrix Eye_tau = Eye_q ;
                                    yarp::sig::Matrix Zeros_fc_q(size_fc, size_q) ;
                                    yarp::sig::Matrix Zeros_q_fc(size_q, size_fc) ;
                                    yarp::sig::Matrix Zeros_q_q(size_q, size_q) ;
                                    yarp::sig::Matrix Zeros_6_q( 6 , size_q) ;
                                    yarp::sig::Matrix Zeros_q_6(size_q, 6 ) ;
                                    Zeros_fc_q.zero();
                                    Zeros_q_fc.zero();
                                    Zeros_q_q.zero();
                                    Zeros_6_q.zero();    
                                    Zeros_q_6.zero();  
                                    // Setting the first block-row of the FLMM
                                    FLMM.setSubmatrix( Eye_fc                           , 0 , 0                   ) ;
                                    FLMM.setSubmatrix( Zeros_fc_q                       , 0 , size_fc             ) ;
                                    FLMM.setSubmatrix( -1.0 * K_c*S_c.transposed()       , 0 , size_fc+size_q      ) ;
                                    FLMM.setSubmatrix( -1.0 * K_c*J_c                    , 0 , size_fc+size_q+6    ) ;
                                    FLMM.setSubmatrix( Zeros_fc_q                       , 0 , size_fc+2*size_q+6  ) ;
                                    // Setting the second block-row of the FLMM
                                    FLMM.setSubmatrix( -1.0*J_c.transposed()        , size_fc , 0                  ) ;
                                    FLMM.setSubmatrix( Eye_tau                      , size_fc , size_fc            ) ;
                                    FLMM.setSubmatrix( -1.0 * U_j                   , size_fc , size_fc+size_q     ) ;
                                    FLMM.setSubmatrix( -1.0 * Q_j                   , size_fc , size_fc+size_q+6   ) ;
                                    FLMM.setSubmatrix( Zeros_q_q                    , size_fc , size_fc+2*size_q+6 ) ;
                                    // Setting the third block-row of the FLMM
                                    FLMM.setSubmatrix( -1.0*S_c      , size_fc +size_q , 0                  ) ;
                                    FLMM.setSubmatrix( Zeros_6_q     , size_fc +size_q , size_fc            ) ;
                                    FLMM.setSubmatrix( -1.0 * U_s    , size_fc +size_q , size_fc+size_q     ) ;
                                    FLMM.setSubmatrix( -1.0 * Q_s    , size_fc +size_q , size_fc+size_q+6   ) ;
                                    FLMM.setSubmatrix( Zeros_6_q     , size_fc +size_q , size_fc+2*size_q+6 ) ;

                                    // Setting the fourth block-row of the FLMM
                                    FLMM.setSubmatrix( Zeros_q_fc  , size_fc +size_q +6  , 0                  ) ;
                                    FLMM.setSubmatrix( Eye_tau     , size_fc +size_q +6  , size_fc            ) ;
                                    FLMM.setSubmatrix( Zeros_q_6   , size_fc +size_q +6  , size_fc+size_q     ) ;
                                    FLMM.setSubmatrix( K_q          , size_fc +size_q +6  , size_fc+size_q+6   ) ;
                                    FLMM.setSubmatrix( -1.0*K_q     , size_fc +size_q +6  , size_fc+2*size_q+6 ) ;
                                    return FLMM ;                          
        
                            }
                              
                                 //------------------------------------------------------------------------------------
     /**
     * @brief  Rf_ext computes the joints-forces map for a compliant humanoid robot 
     * @param  J_c is a contacts x joints yarp matrix describing the body Jacobian of the humanoid
     * @param  S_c is a 6 x contacts yarp matrix describing the body Stance matrix of the humanoid
     * @param  Q_j is a joints x joints yarp matrix about the derivative of the Jacobian 
     * @param  Q_s is a contacts x joints yarp matrix about the derivative of the Jacobian 
     * @param  U_j is a joints x 6 yarp matrix about the derivative of the Jacobian 
     * @param  U_s is a 6 x 6 yarp matrix about the derivative of the Jacobian 
     * @param  K_c is a contacts x contacts yarp matrix describing the contact stiffness matrix
     * @param  K_q is a joints x joints yarp matrix describing the joint stiffness matrix
     * @return Rf_ext computes the joint-forces map
     */
    yarp::sig::Matrix Rf_ext(   const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_j,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_j,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c,
                                const yarp::sig::Matrix& K_q
                              ) {
                                    yarp::sig::Matrix Q_j_1 = -1.0*Q_j - 1.0* J_c.transposed()*K_c*J_c  ;
                                    yarp::sig::Matrix U_j_1 = -1.0*U_j -1.0*J_c.transposed() *K_c*S_c.transposed() ;    
                                    yarp::sig::Matrix Q_s_1 =  -1.0*Q_s-1.0*S_c*K_c*J_c  ;
                                    yarp::sig::Matrix U_s_1 = -1.0*U_s-1.0*S_c*K_c*S_c.transposed() ;    
                                    yarp::sig::Matrix L = yarp::math::luinv(U_s_1)* Q_s_1 ;
                                    yarp::sig::Matrix M = Q_j_1-U_j_1*L ;    
                                    yarp::sig::Matrix H = K_q-M ;
                                    yarp::sig::Matrix F = -1.0*yarp::math::luinv(H)*K_q ;    
                                    yarp::sig::Matrix E = -1.0*K_c* S_c.transposed()* L *F ;
                                    yarp::sig::Matrix R_f_1 = E+K_c*J_c*F  ;
                                    return R_f_1 ;
                            }
    //------------------------------------------------------------------------------------
     /**
     * @brief  FLMM_redu computes a basic version of the FLMM for a rigid humanoid robot 
     * @param  J_c is a contacts x joints yarp matrix describing the body Jacobian of the humanoid
     * @param  S_c is a 6 x contacts yarp matrix describing the body Stance matrix of the humanoid
     * @param  Q_s is a contacts x joints yarp matrix about the derivative of the Jacobian 
     * @param  U_s is a 6 x 6 yarp matrix about the derivative of the Jacobian 
     * @param  K_c is a contacts x contacts yarp matrix describing the contact stiffness matrix
     * @return FLMM_ext is the Fundamental Loco-Manipulation Matrix
     */
    yarp::sig::Matrix FLMM_redu( const yarp::sig::Matrix& J_c ,
                                 const yarp::sig::Matrix& S_c ,
                                 const yarp::sig::Matrix& Q_s,
                                 const yarp::sig::Matrix& U_s,
                                 const yarp::sig::Matrix& K_c
                              ) {
                                    int size_fc  = J_c.rows();
                                    int size_q  = J_c.cols() ;
                                    int r_FLMM = size_fc + 6  ;
                                    int c_FLMM = size_fc + 6 + size_q;
                                    yarp::sig::Matrix FLMM(r_FLMM, c_FLMM) ;    
                                    yarp::sig::Matrix Eye_fc(size_fc, size_fc) ;
                                    Eye_fc.eye() ;
                                    // Setting the first block-row of the FLMM
                                    FLMM.setSubmatrix( Eye_fc                      , 0 , 0                   ) ;
                                    FLMM.setSubmatrix( -1.0 * K_c*S_c.transposed() , 0 , size_fc             ) ;
                                    FLMM.setSubmatrix( -1.0 * K_c*J_c              , 0 , size_fc+6      ) ;
                                    // Setting the second block-row of the FLMM
                                    FLMM.setSubmatrix( -1.0 * S_c    , size_fc , 0                  ) ;
                                    FLMM.setSubmatrix( -1.0 * U_s    , size_fc , size_fc            ) ;
                                    FLMM.setSubmatrix( -1.0 * Q_s    , size_fc , size_fc+6     ) ;    
                                    return FLMM ;
                            }                            
     //------------------------------------------------------------------------------------
     /**
     * @brief  Rf_redu computes the joints-forces map for a rigid humanoid robot 
     * @param  J_c is a contacts x joints yarp matrix describing the body Jacobian of the humanoid
     * @param  S_c is a 6 x contacts yarp matrix describing the body Stance matrix of the humanoid
     * @param  Q_s is a contacts x joints yarp matrix about the derivative of the Jacobian 
     * @param  U_s is a 6 x 6 yarp matrix about the derivative of the Jacobian 
     * @param  K_c is a contacts x contacts yarp matrix describing the contact stiffness matrix
     * @return Rf_redu is the Fundamental Loco-Manipulation Matrix
     */
    yarp::sig::Matrix Rf_redu(  const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c
                              ) {
                                    yarp::sig::Matrix Q_s_1 =  Q_s +  S_c*K_c*J_c  ;
                                    yarp::sig::Matrix U_s_1 =  U_s + S_c*K_c*S_c.transposed() ;    
                                    yarp::sig::Matrix L = yarp::math::luinv(U_s_1)* Q_s_1 ; //Pinv_trunc_SVD(U_s_1, 1E4)* Q_s_1; //
                                    yarp::sig::Matrix R_f_1 = - 1.0*K_c*J_c + K_c*S_c.transposed()* L  ; 
                                    return R_f_1 ;
                            }                             
                              
       //------------------------------------------------------------------------------------
     /**
     * @brief  Ru_redu computes the joints-movements map for a a rigid humanoid robot 
     * @param  J_c is a contacts x joints yarp matrix describing the body Jacobian of the humanoid
     * @param  S_c is a 6 x contacts yarp matrix describing the body Stance matrix of the humanoid
     * @param  Q_s is a contacts x joints yarp matrix about the derivative of the Jacobian 
     * @param  U_s is a 6 x 6 yarp matrix about the derivative of the Jacobian 
     * @param  K_c is a contacts x contacts yarp matrix describing the contact stiffness matrix
     * @return Ru_redu is the Fundamental Loco-Manipulation Matrix
     */
    yarp::sig::Matrix Ru_redu(  const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c
                              ) {
                                    yarp::sig::Matrix Q_s_1 =  -1.0*Q_s-1.0*S_c*K_c*J_c  ;
                                    yarp::sig::Matrix U_s_1 = -1.0*U_s-1.0*S_c*K_c*S_c.transposed() ;    
                                    yarp::sig::Matrix L = yarp::math::luinv(U_s_1)* Q_s_1 ;
                                    return L ;
                            }                                             
                              
                              
        }
}
#endif