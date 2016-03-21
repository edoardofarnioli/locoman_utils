#ifndef _KINEMATICS_HPP_
#define _KINEMATICS_HPP_

#include <cstdlib>
#include <iostream>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include <locoman/utils/screws.hpp>
// #include <locoman/utils/kinematics.hpp>
#include <locoman/utils/kinetostatics.hpp>
#include <locoman/utils/locoman_utils.hpp>
#include <locoman/utils/algebra.hpp>

namespace locoman {
        namespace utils {
            
            // pre declaration
            yarp::sig::Vector getTrasl( const yarp::sig::Matrix& T_ab);
            yarp::sig::Matrix getRot( const yarp::sig::Matrix& T_ab);
            
                 //-----------------------------------------------------------------
     /**
     * @brief  SkewToVect transforms a cross product matrix into the original vector
     * @param  Skew is a 3x3 yarp matrix describing an cross product
     * @return 3x1 yarp vector
     */
    yarp::sig::Vector SkewToVect( const yarp::sig::Matrix& Skew) {
                                            yarp::sig::Vector Vect(3) ;
                                            Vect[0] = Skew[2][1] ;
                                            Vect[1] = Skew[0][2] ;
                                            Vect[2] = Skew[1][0] ;
                                            return Vect ;
                                                }
    
       /**
     * @brief  Rot2Quat computes the quaterion components given the rotation matrix
     * @param  Rot is 3 x 3 rotation matrix
     * @return a vector of 4 elements; the first element is the scalar part of the quaternion
     */
    yarp::sig::Vector Rot2Quat( const yarp::sig::Matrix& Rot) {
                                        double tol = 0.001 ;
                                        double r_11 = Rot[0][0] ;
                                        double r_12 = Rot[0][1] ;
                                        double r_13 = Rot[0][2] ;
                                        double r_21 = Rot[1][0] ;
                                        double r_22 = Rot[1][1] ;
                                        double r_23 = Rot[1][2] ;
                                        double r_31 = Rot[2][0] ;
                                        double r_32 = Rot[2][1] ;
                                        double r_33 = Rot[2][2] ;
                                        yarp::sig::Vector quat(4,0.0) ;
                                        quat[0] = (1.0/2.0)*sqrt(r_11+ r_22+ r_33+ 1) ; 
                                        double a1 = yarp::math::sign(r_32-r_23) ;
                                        double a2 = yarp::math::sign(r_13-r_31) ;
                                        double a3 = yarp::math::sign(r_21-r_12) ;
                                        if(a1==0){a1 = 1 ;};
                                        if(a2==0){a2 = 1 ;};
                                        if(a3==0){a3 = 1 ;}; 
                                        //--
                                        double to_sq_1 = r_11- r_22- r_33 + 1 ;
                                        double to_sq_2 = r_22- r_33- r_11 + 1 ;
                                        double to_sq_3 = r_33- r_11- r_22 + 1;
                                        double sq_1 ;
                                        double sq_2 ;
                                        double sq_3 ;
                                        if(to_sq_1>tol){sq_1 = sqrt(to_sq_1) ; }
                                        else {sq_1 = 0 ;}
                                        if(to_sq_2>tol){sq_2 = sqrt(to_sq_2) ; }
                                        else {sq_2 = 0 ;}
                                        if(to_sq_3>tol){sq_3 = sqrt(to_sq_3) ; }
                                        else {sq_3 = 0 ;}
                                        quat[1] = (1.0/2.0)*a1*sq_1 ; // sqrt(r_11- r_22- r_33 + 1) ; 
                                        quat[2] = (1.0/2.0)*a2*sq_2 ; // sqrt(r_22- r_33- r_11 + 1) ; 
                                        quat[3] = (1.0/2.0)*a3*sq_3 ; // sqrt(r_33- r_11- r_22 + 1) ; 
                                        //  std::cout << "  a1 = " <<  std::endl << a1  << std::endl;   
                                        return quat ;
                                            }          
     
     /**
     * @brief  Orient_Error computes the orientation error based on the quaternion representation of the rotation matrices
     * @param  Rot_des is 3 x 3 rotation matrix representing the desired orientation
     * @param  Rot_cur is 3 x 3 rotation matrix representing the current orientation
     * @return a vector of 3 elements representing the orientation error
     */
    yarp::sig::Vector Orient_Error( const yarp::sig::Matrix& Rot_des,
                                    const yarp::sig::Matrix& Rot_cur) {
                                        yarp::sig::Vector q_des = Rot2Quat(Rot_des) ;
                                        yarp::sig::Vector q_cur = Rot2Quat(Rot_cur) ;
                                        double eta_des = q_des[0] ;
                                        yarp::sig::Vector eps_des = q_des.subVector(1,q_des.length()-1);
                                        double eta_cur = q_cur[0] ;
                                        yarp::sig::Vector eps_cur = q_cur.subVector(1,q_cur.length()-1);
                                        yarp::sig::Vector e_o = eta_cur*eps_des- eta_des*eps_cur - 1.0* yarp::math::crossProductMatrix(eps_des)*eps_cur  ;
                                        return e_o ;
                                    }         

     /**
     * @brief  Inv_quaternion computes the inverse of a given quaternion
     * @param  quat is 4 elements vector describing a quaternion
     * @return a vector of 4 elements representing the inverse of the input quaternion
     */
    yarp::sig::Vector Inv_quaternion( const yarp::sig::Vector& quat ) {
                                            yarp::sig::Vector e_inv = -1.0*quat.subVector(1,3) ;
                                            yarp::sig::Vector quat_inv(4,0.0) ;  
                                            quat_inv[0] = quat[0] ;
                                            quat_inv.setSubvector(1, e_inv) ;
                                            return quat_inv; 
                                                }  

     /**
     * @brief  Computes the product of two quaternions
     * @param  quat_1 and quat_2 are 4 elements vectors describing quaternions
     * @return a vector of 4 elements representing the product of the two input quaternions
     */
    yarp::sig::Vector quat_Product( const yarp::sig::Vector& quat_1 , 
                                    const yarp::sig::Vector& quat_2) {
                                          double eta_1 = quat_1[0] ;
                                            double eta_2 = quat_2[0];
                                            yarp::sig::Vector eps_1 = quat_1.subVector(1,3) ;
                                            yarp::sig::Vector eps_2 = quat_2.subVector(1,3) ;
                                            yarp::sig::Vector quat_prod(4,0.0) ;  
                                            quat_prod[0] = eta_1*eta_2 -1.0*yarp::math::dot(eps_1,eps_2) ;
                                            yarp::sig::Vector eps_prod =  eta_1*eps_2 + eta_2*eps_1 + yarp::math::crossProductMatrix(eps_1)*eps_2 ;
                                            quat_prod.setSubvector(1, eps_prod) ;
                                            return quat_prod ; 
                                    }  
    
    /**
     * @brief  Computes the rotation matrix about x axis
     * @param  phi_x rotation angle
     * @return a 3x3 yarp Matrix
     */
    yarp::sig::Matrix Rot_x( const double phi_x ) {    
                                    yarp::sig::Matrix Rotx(3,3) ;
                                    Rotx.eye() ;
                                    Rotx[1][1] = cos(phi_x) ;
                                    Rotx[2][1] = sin(phi_x) ;
                                    Rotx[1][2] =-sin(phi_x) ;
                                    Rotx[2][2] = cos(phi_x) ;
                                    return Rotx ;
                                            }     
        
    /**
     * @brief  Computes the rotation matrix about y axis
     * @param  theta_y rotation angle
     * @return a 3x3 yarp Matrix
     */   
    yarp::sig::Matrix Rot_y( const double theta_y ) {
                                    yarp::sig::Matrix Roty(3,3) ;
                                    Roty.eye() ;
                                    Roty[0][0] = cos(theta_y) ;
                                    Roty[0][2] = sin(theta_y) ;
                                    Roty[2][0] =-sin(theta_y) ;
                                    Roty[2][2] = cos(theta_y) ;
                                    return Roty ;
                                        }
    
     /**
     * @brief  Computes the rotation matrix about z axis
     * @param  psi_z rotation angle
     * @return a 3x3 yarp Matrix
     */
    yarp::sig::Matrix Rot_z( const double psi_z) {
                                    yarp::sig::Matrix Rotz(3,3) ;
                                    Rotz.eye() ;
                                    Rotz[0][0] = cos(psi_z) ;
                                    Rotz[1][0] = sin(psi_z) ;
                                    Rotz[0][1] =-sin(psi_z) ;
                                    Rotz[1][1] = cos(psi_z) ;
                                    return Rotz ;
                                        }      
     
    
         /**
     * @brief  WB_Cartesian_Tasks computes the whole-body displacement for achieving the 
     *         desired configurations of the End Effectors and of the CoM (all together)
     * @param  T_l1_foot_des desired pose of the left foot (first contact) in local frame
     * @param  T_r1_foot_des desired pose of the right foot (first contact) in local frame
     * @param  T_l_hand_des  desired pose of the left hand in local frame
     * @param  T_r_hand_des  desired pose of the right hand in local frame
     * @param  CoM_waist_cmd desired position of the CoM w.r.t the wrist frame
     * 
     * @param  J_l1_foot_body is a (6 x (joints+6)) body Jacobian of the left  foot (first contact) 
     * @param  J_r1_foot_body is a (6 x (joints+6)) body Jacobian of the right foot (first contact) 
     * @param  J_l_hand_body  is a (6 x (joints+6)) body Jacobian of the left  hand
     * @param  J_r_hand_body  is a (6 x (joints+6)) body Jacobian of the right hand 
     * @param  J_waist_CoM    is a (3 x (joints+6)) Jacobian of the CoM expressed in {Waist} frame
     * @return the desired delta_q vector 
     */
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
                                        ) {
                                        yarp::sig::Matrix Eye_3(3,3); Eye_3.eye() ;
                                        // filtering the error
                                        double  max_trasl_err  = 0.001  ;  // maximum displacement allowed in a loop
                                        double  max_orient_err = 0.01   ;  // maximum rotation allowed in a loop
                                        double  min_trasl_err  = 0.00002 ;  // minimum displacement, otherwise is approximated to zero
                                        double  min_orient_err = 0.001  ;  // minimum rotation, otherwise is approximated to zero
                                        
                                        yarp::sig::Vector trasl_l_hand_err   = locoman::utils::getTrasl(T_l_hand_des) ;
                                        yarp::sig::Vector orient_l_hand_err  = Orient_Error( locoman::utils::getRot(T_l_hand_des ) , Eye_3 )  ;
                                        yarp::sig::Vector trasl_r_hand_err   = locoman::utils::getTrasl(T_r_hand_des) ;
                                        yarp::sig::Vector orient_r_hand_err  = Orient_Error( locoman::utils::getRot(T_r_hand_des ) , Eye_3 )  ;
                                        yarp::sig::Vector trasl_l1_foot_err  = locoman::utils::getTrasl(T_l1_foot_des);
                                        yarp::sig::Vector orient_l1_foot_err = Orient_Error( locoman::utils::getRot(T_l1_foot_des), Eye_3 )   ;
                                        yarp::sig::Vector trasl_r1_foot_err  = locoman::utils::getTrasl(T_r1_foot_des);
                                        yarp::sig::Vector orient_r1_foot_err = Orient_Error( locoman::utils::getRot(T_r1_foot_des ) , Eye_3 ) ;
                                        yarp::sig::Vector CoM_err = CoM_waist_cmd  ;

                                        if(yarp::math::norm(trasl_l_hand_err)  > max_trasl_err ) {trasl_l_hand_err   = max_trasl_err  * trasl_l_hand_err/(yarp::math::norm(trasl_l_hand_err)    ) ;}  
                                        if(yarp::math::norm(orient_l_hand_err) > max_orient_err) {orient_l_hand_err  = max_orient_err * orient_l_hand_err/(yarp::math::norm(orient_l_hand_err)  ) ;}  
                                        if(yarp::math::norm(trasl_r_hand_err)  > max_trasl_err ) {trasl_r_hand_err   = max_trasl_err  * trasl_r_hand_err/(yarp::math::norm(trasl_r_hand_err)    ) ;}
                                        if(yarp::math::norm(orient_r_hand_err) > max_orient_err) {orient_r_hand_err  = max_orient_err * orient_r_hand_err/(yarp::math::norm(orient_r_hand_err)  ) ;}
                                        if(yarp::math::norm(trasl_l1_foot_err) > max_trasl_err ) {trasl_l1_foot_err  = max_trasl_err  * trasl_l1_foot_err/(yarp::math::norm(trasl_l1_foot_err)  ) ;}
                                        if(yarp::math::norm(orient_l1_foot_err)> max_orient_err) {orient_l1_foot_err = max_orient_err * orient_l1_foot_err/(yarp::math::norm(orient_l1_foot_err)) ;}
                                        if(yarp::math::norm(trasl_r1_foot_err) > max_trasl_err ) {trasl_r1_foot_err  = max_trasl_err  * trasl_r1_foot_err/(yarp::math::norm(trasl_r1_foot_err)  ) ;}
                                        if(yarp::math::norm(orient_r1_foot_err)> max_orient_err) {orient_r1_foot_err = 0.0*orient_r1_foot_err ;}
                                        if(yarp::math::norm(CoM_err) > max_trasl_err) {CoM_err = max_trasl_err*CoM_err ;}
                                        
                                        if(yarp::math::norm(trasl_l_hand_err)<min_trasl_err) {trasl_l_hand_err = 0.0*trasl_l_hand_err ;}
                                        if(yarp::math::norm(orient_l_hand_err)<min_orient_err) {orient_l_hand_err = 0.0*orient_l_hand_err ;}  
                                        if(yarp::math::norm(trasl_r_hand_err)<min_trasl_err) {trasl_r_hand_err = 0.0*trasl_r_hand_err ;}
                                        if(yarp::math::norm(orient_r_hand_err)<min_orient_err) {orient_r_hand_err = 0.0*orient_r_hand_err ;}
                                        if(yarp::math::norm(trasl_l1_foot_err)<min_trasl_err) {trasl_l1_foot_err = 0.0*trasl_l1_foot_err ;}
                                        if(yarp::math::norm(orient_l1_foot_err)<min_orient_err) {orient_l1_foot_err = 0.0*orient_l1_foot_err ;}
                                        if(yarp::math::norm(trasl_r1_foot_err)<min_trasl_err) {trasl_r1_foot_err = 0.0*trasl_r1_foot_err ;}
                                        if(yarp::math::norm(orient_r1_foot_err)<min_orient_err) {orient_r1_foot_err = 0.0*orient_r1_foot_err ;}
                                        if(yarp::math::norm(CoM_err)<min_trasl_err) {CoM_err = 0.0*CoM_err ;}
                                        
                                        // Error Vector
                                        yarp::sig::Vector d_C(27, 0.0) ;
                                        d_C.setSubvector(0 , trasl_l_hand_err    ) ;
                                        d_C.setSubvector(3 , orient_l_hand_err   ) ; 
                                        d_C.setSubvector(6 , trasl_r_hand_err    ) ;
                                        d_C.setSubvector(9 , orient_r_hand_err   ) ;
                                        d_C.setSubvector(12, trasl_l1_foot_err   ) ;
                                        d_C.setSubvector(15, orient_l1_foot_err  ) ;
                                        d_C.setSubvector(18, trasl_r1_foot_err   ) ;
                                        d_C.setSubvector(21, orient_r1_foot_err  ) ;
                                        d_C.setSubvector(24, CoM_err             ) ;
                                        
                                        //   std::cout << " d_C = " <<  std::endl << d_C.toString() << std::endl; 

                                        // Whole Jacobian 
                                        yarp::sig::Matrix Whole_Jac(d_C.length(), J_l_hand_body.cols() ) ; 
                                        Whole_Jac.zero(); 
                                        Whole_Jac.setSubmatrix( J_l_hand_body , 0  , 0 ) ;
                                        Whole_Jac.setSubmatrix( J_r_hand_body , 6  , 0 ) ;
                                        Whole_Jac.setSubmatrix( J_l1_foot_body, 12 , 0 ) ;
                                        Whole_Jac.setSubmatrix( J_r1_foot_body, 18 , 0 ) ;
                                        Whole_Jac.setSubmatrix( J_waist_CoM   , 24 , 0 ) ;
                                        
                                        //pseudoinverse
                                        yarp::sig::Vector d_q_long = locoman::utils::Pinv_trunc_SVD( Whole_Jac , 1E-10 ) * d_C ;

                                        //selections
                                        yarp::sig::Vector d_q_move = d_q_long.subVector(6, d_q_long.length()-1) ;
                                        
                                        return d_q_move ; 
                                        } 

        }
}

#endif