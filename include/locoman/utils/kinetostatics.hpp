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
	    yarp::sig::Matrix getRot(const yarp::sig::Matrix& T_ab) ;
	    yarp::sig::Vector WB_Cartesian_Tasks( 
                            const yarp::sig::Matrix& T_l_hand_des,
                            const yarp::sig::Matrix& T_r_hand_des,
                            const yarp::sig::Matrix& T_l1_foot_des ,
                            const yarp::sig::Matrix& T_r1_foot_des ,
                            const yarp::sig::Vector& CoM_waist_cmd ,
                            const yarp::sig::Matrix& J_l_hand_body ,
                            const yarp::sig::Matrix& J_r_hand_body ,
                            const yarp::sig::Matrix& J_l1_foot_body ,
                            const yarp::sig::Matrix& J_r1_foot_body ,
                            const yarp::sig::Matrix& J_waist_CoM 
                                        ) ;
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
                            
                            
         
    //RobotUtils& robot 
                            
   /**
     * @brief  Q_mg compute the derivative of the gravitational term
     * @param  q_actual current robot configuration configuration 
     * @param  mg_vector contact force vector (3 force elements)
     * @param  T_aw_w_0 is the homogeneous transformation between {AW} and {W} (Auxiliary wolrd and wolrd)
     * @param  robot is the robot model
     * @return (q+6)x(q+6) yarp matrix 
     **/
    yarp::sig::Matrix Q_mg( const yarp::sig::Vector& q_actual , 
                            const yarp::sig::Vector& mg_vector,
			    const yarp::sig::Matrix& T_aw_w_0 , 
			    RobotUtils& robot) {
				  robot.idynutils.updateiDyn3Model( q_actual, true );
                                  yarp::sig::Matrix Q_mg_matrix( q_actual.length()+6 , q_actual.length()+6 );
				  Q_mg_matrix.zero();
				  yarp::sig::Vector q_h = q_actual ;
                                  double h_incremental = std::pow(std::numeric_limits<double>::epsilon(), 1.0/2.0);  
				  yarp::sig::Matrix J_com_w(6, q_actual.length() + 6 ) ;
				  yarp::sig::Matrix J_com_w_redu(6, q_actual.length()+6);
				  yarp::sig::Matrix J_com_aw(6, q_actual.length()+6);
				  yarp::sig::Matrix J_com_w_q_i(6, q_actual.length() + 6 ) ;
				  yarp::sig::Matrix J_com_w_redu_q_i(6, q_actual.length()+6);
				  yarp::sig::Matrix J_com_aw_q_i(6, q_actual.length()+6);				  
				  
				  yarp::sig::Vector tau_mg_0( q_actual.length()+6, 0.0) ; //= J_com_aw.transposed()*mg_vector ;
				  yarp::sig::Vector tau_mg_q_i(q_actual.length()+6, 0.0 ) ;// = tau_mg_0 ;
				  yarp::sig::Vector tau_mg_differential_i( q_actual.length()+6, 0.0 )  ;
				  
				  yarp::sig::Matrix Rot_aw_w_0 = locoman::utils::getRot(T_aw_w_0)  ;
  
				  // Calcolare i jacobiani nella configurazione attuale
				  robot.idynutils.iDyn3_model.getCOMJacobian(J_com_w) ;
				  J_com_w_redu = J_com_w.submatrix(0,2 , 0 , J_com_w.cols()-1 ) ;  
				  J_com_aw     = Rot_aw_w_0 *J_com_w_redu;   
				  tau_mg_0 = J_com_aw.transposed()*mg_vector ;
  
                                  for ( int i = 0  ; i<q_actual.length()  ; i++ )  // i<(J_waist_l_c1_spa_0.cols()-1)
                                        {
					    tau_mg_q_i.zero();
					    q_h = q_actual ;
					    q_h[i] += h_incremental ;
					    robot.idynutils.updateiDyn3Model( q_h, true );
					    robot.idynutils.iDyn3_model.getCOMJacobian(J_com_w_q_i) ;
					    J_com_w_redu_q_i = J_com_w_q_i.submatrix(0,2 , 0 , J_com_w_q_i.cols()-1 ) ;  
					    J_com_aw_q_i = Rot_aw_w_0 *J_com_w_redu_q_i ; 
					    tau_mg_q_i = J_com_aw_q_i.transposed() * mg_vector ;
					    tau_mg_differential_i = (tau_mg_q_i - tau_mg_0)/h_incremental ;
					    Q_mg_matrix.setCol( i, tau_mg_differential_i ) ;
   
//                                         Q_ci_col_i =  ( locoman::utils::D_Jacob_spa_i(J_spa_i, (i+1)) ).transposed()*
//                                                         ( locoman::utils::Adjoint(
//                                                             locoman::utils::iHomogeneous(T_a_ci))).transposed()*B *f_ci  +
//                                                         J_spa_i.transposed()* 
//                                                         ( locoman::utils::Adjoint( locoman::utils::iHomogeneous(T_a_ci) )*
//                                                         locoman::utils::ad_lie( -1.0*J_spa_i.getCol(i))).transposed()*  B *f_ci;
//                                         Q_c_i.setCol(i, Q_ci_col_i ) ;       
                                        }   
                                    return Q_mg_matrix ;
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
                            
 yarp::sig::Vector rg_hand_dw_touch_simple( RobotUtils& robot ,
		                  yarp::sig::Vector& q_current, 
				  const yarp::sig::Matrix& Big_J, 
				  const int r_ankle_index,
				  const int r_hand_index,
				  double f_limit = 5.0,
				  double step_limit = 0.1
                              ) {
      		      // simple it means that we are on flat terrain
                     // absolute z-axis = z axis of the feet
			  yarp::sig::Vector d_q_move = 0.0*q_current ;
			  yarp::sig::Vector ft_r_wrist(6,0.0) ;
			  robot.senseftSensor("r_arm_ft", ft_r_wrist) ;  
			  f_limit = -1.0*abs(f_limit);
			  
			  if(ft_r_wrist[2]>f_limit){ 
			  yarp::sig::Matrix T_w_r_ankle_0 = robot.idynutils.iDyn3_model.getPosition(r_ankle_index) ;
			  yarp::sig::Matrix T_w_r_hand_0  = robot.idynutils.iDyn3_model.getPosition(r_hand_index ) ;   

			  yarp::sig::Matrix T_rfoot_rhand = locoman::utils::iHomogeneous(T_w_r_hand_0)*T_w_r_ankle_0;
			  yarp::sig::Matrix T_rfoot_target =T_rfoot_rhand ;
			  T_rfoot_target[2][3] -=0.1 ;
			  yarp::sig::Matrix T_rhand_target = locoman::utils::iHomogeneous(T_rfoot_rhand)*T_rfoot_target ;
 
			  yarp::sig::Matrix Eye_4(4,4) ;
			  Eye_4.eye();
			  yarp::sig::Vector zero_3(3,0.0) ;

// 			  yarp::sig::Matrix T_rg_dw(4,4) ;
// 			  T_rg_dw.eye() ;
// 			  T_rg_dw[2][3] = -0.1 ;
	      
			  d_q_move = locoman::utils::WB_Cartesian_Tasks( 
									  Eye_4,             // T_l_hand_des,
									  T_rhand_target,             // T_r_hand_des,
									  Eye_4,             // T_l1_foot_des ,
									  Eye_4,             // T_r1_foot_des ,
									  zero_3,            //CoM_err ,
									  Big_J.submatrix(0,5,0,Big_J.cols()-1) ,
									  Big_J.submatrix(6,11,0,Big_J.cols()-1) ,
									  Big_J.submatrix(12,17,0,Big_J.cols()-1) ,
									  Big_J.submatrix(18,23,0,Big_J.cols()-1) ,
									  Big_J.submatrix(24,Big_J.rows()-1,0,Big_J.cols()-1) 
										      ) ;
			   if(norm(d_q_move) > step_limit) { d_q_move =  step_limit*d_q_move/ norm(d_q_move) ; } 
			
		      }
		      else  {d_q_move = 0.0*q_current ; }
		      
		      return d_q_move ;
				  
                            }     
  
  
  yarp::sig::Vector lf_hand_dw_touch_simple( RobotUtils& robot ,
		                  yarp::sig::Vector& q_current, 
				  const yarp::sig::Matrix& Big_J, 
				  const int l_ankle_index,
				  const int l_hand_index,
				  double f_limit = 5.0,
				  double step_limit = 0.1
                              ) {
      		      // simple it means that we are on flat terrain
                     // absolute z-axis = z axis of the feet
			  yarp::sig::Vector d_q_move = 0.0*q_current ;
			  yarp::sig::Vector ft_l_wrist(6,0.0) ;
			  robot.senseftSensor("l_arm_ft", ft_l_wrist) ;  
			  f_limit = -1.0*abs(f_limit);
			  
			  if(ft_l_wrist[2]>f_limit){ 
			  yarp::sig::Matrix T_w_l_ankle_0 = robot.idynutils.iDyn3_model.getPosition(l_ankle_index) ;
			  yarp::sig::Matrix T_w_l_hand_0  = robot.idynutils.iDyn3_model.getPosition(l_hand_index ) ;   

			  yarp::sig::Matrix T_lfoot_rhand = locoman::utils::iHomogeneous(T_w_l_hand_0)*T_w_l_ankle_0;
			  yarp::sig::Matrix T_lfoot_target =T_lfoot_rhand ;
			  T_lfoot_target[2][3] -=0.1 ;
			  yarp::sig::Matrix T_lhand_target = locoman::utils::iHomogeneous(T_lfoot_rhand)*T_lfoot_target ;
 
			  yarp::sig::Matrix Eye_4(4,4) ;
			  Eye_4.eye();
			  yarp::sig::Vector zero_3(3,0.0) ;

	      
			  d_q_move = locoman::utils::WB_Cartesian_Tasks( 
									  T_lhand_target,             // T_l_hand_des,
									  Eye_4,             // T_r_hand_des,
									  Eye_4,             // T_l1_foot_des ,
									  Eye_4,             // T_r1_foot_des ,
									  zero_3,            //CoM_err ,
									  Big_J.submatrix(0,5,0,Big_J.cols()-1) ,
									  Big_J.submatrix(6,11,0,Big_J.cols()-1) ,
									  Big_J.submatrix(12,17,0,Big_J.cols()-1) ,
									  Big_J.submatrix(18,23,0,Big_J.cols()-1) ,
									  Big_J.submatrix(24,Big_J.rows()-1,0,Big_J.cols()-1) 
										      ) ;
			   if(norm(d_q_move) > step_limit) { d_q_move =  step_limit*d_q_move/ norm(d_q_move) ; } 
			
		      }
		      else  {d_q_move = 0.0*q_current ; }
		      
		      return d_q_move ;
				  
                            }  
     
     yarp::sig::Vector hands_dw_touch_simple( RobotUtils& robot ,
					    yarp::sig::Vector& q_current, 
					    const yarp::sig::Matrix& Big_J, 
					    const int l_ankle_index,
					    const int r_ankle_index,
					    const int l_hand_index,
					    const int r_hand_index,
					    double f_limit = 5.0,
					    double step_limit = 0.1
                              ) {
      		      // simple it means that we are on flat terrain
                     // absolute z-axis = z axis of the feet
			  yarp::sig::Vector d_q_move = 0.0*q_current ;
			  
			  f_limit = -1.0*abs(f_limit);
			  
			  yarp::sig::Vector ft_l_wrist(6,0.0) ;
 			  yarp::sig::Vector ft_r_wrist(6,0.0) ;
			  robot.senseftSensor("l_arm_ft", ft_l_wrist) ;  
			  robot.senseftSensor("r_arm_ft", ft_r_wrist) ;  
			  
			  yarp::sig::Matrix T_lhand_target(4,4) ;
			  T_lhand_target.eye() ;
// 			  
			  if(ft_l_wrist[2]>f_limit){ 
			  yarp::sig::Matrix T_w_l_ankle_0 = robot.idynutils.iDyn3_model.getPosition(l_ankle_index) ;
			  yarp::sig::Matrix T_w_l_hand_0  = robot.idynutils.iDyn3_model.getPosition(l_hand_index ) ;   

			  yarp::sig::Matrix T_lfoot_rhand = locoman::utils::iHomogeneous(T_w_l_hand_0)*T_w_l_ankle_0;
			  yarp::sig::Matrix T_lfoot_target =T_lfoot_rhand ;
			  T_lfoot_target[2][3] -=0.1 ;
			  T_lhand_target = locoman::utils::iHomogeneous(T_lfoot_rhand)*T_lfoot_target ;
			  }
			  
			  yarp::sig::Matrix T_rhand_target(4,4) ;
			  T_rhand_target.eye() ;
			  if(ft_r_wrist[2]>f_limit){ 
			  yarp::sig::Matrix T_w_r_ankle_0 = robot.idynutils.iDyn3_model.getPosition(r_ankle_index) ;
			  yarp::sig::Matrix T_w_r_hand_0  = robot.idynutils.iDyn3_model.getPosition(r_hand_index ) ;   

			  yarp::sig::Matrix T_rfoot_rhand = locoman::utils::iHomogeneous(T_w_r_hand_0)*T_w_r_ankle_0;
			  yarp::sig::Matrix T_rfoot_target =T_rfoot_rhand ;
			  T_rfoot_target[2][3] -=0.1 ;
			  T_rhand_target = locoman::utils::iHomogeneous(T_rfoot_rhand)*T_rfoot_target ;
			  }
			  
			  yarp::sig::Matrix Eye_4(4,4) ;
			  Eye_4.eye();
			  yarp::sig::Vector zero_3(3,0.0) ;
 		      if( !(ft_l_wrist[2]<f_limit&& ft_r_wrist[2]<f_limit)){ 
			  d_q_move = locoman::utils::WB_Cartesian_Tasks( 
									  T_lhand_target,             // T_l_hand_des,
									  T_rhand_target,             // T_r_hand_des,
									  Eye_4,             // T_l1_foot_des ,
									  Eye_4,             // T_r1_foot_des ,
									  zero_3,            //CoM_err ,
									  Big_J.submatrix(0,5,0,Big_J.cols()-1) ,
									  Big_J.submatrix(6,11,0,Big_J.cols()-1) ,
									  Big_J.submatrix(12,17,0,Big_J.cols()-1) ,
									  Big_J.submatrix(18,23,0,Big_J.cols()-1) ,
									  Big_J.submatrix(24,Big_J.rows()-1,0,Big_J.cols()-1) 
										      ) ;
			   if(norm(d_q_move) > step_limit) { d_q_move =  step_limit*d_q_move/ norm(d_q_move) ; } 
			}

 		      else 
			  {d_q_move = 0.0*q_current ; }
		      
		      return d_q_move ;
				  
                            }  
     
     
           
                                       
                              
                              
        }
}
#endif