#ifndef _LOCOMAN_UTILS_HPP_
#define _LOCOMAN_UTILS_HPP_

#include <cstdlib>
#include <iostream>

#include <cstdio>
#include <ctime>

#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include <idynutils/RobotUtils.h>

#include <locoman/utils/screws.hpp>
#include <locoman/utils/kinematics.hpp>
#include <locoman/utils/kinetostatics.hpp>
#include <locoman/utils/algebra.hpp>

#define RIGHT_ARM_JOINT_NUM 7
#define LEFT_ARM_JOINT_NUM 7
#define TORSO_JOINT_NUM 3
#define RIGHT_LEG_JOINT_NUM 6
#define LEFT_LEG_JOINT_NUM 6
#define HEAD_JOINT_NUM 2

namespace locoman {
        namespace utils {
            
            const unsigned int& getNumberOfKinematicJoints(const RobotUtils& robot) {
                return robot.left_arm.getNumberOfJoints()  +
                       robot.right_arm.getNumberOfJoints() +
                       robot.left_leg.getNumberOfJoints()  +
                       robot.left_leg.getNumberOfJoints()  +
                       robot.torso.getNumberOfJoints()  +
                       robot.head.getNumberOfJoints()  ;
            }
            
            const unsigned int& getNumberOfActuatedJoints(const RobotUtils& robot) {
                return robot.left_arm.getNumberOfJoints()  +
                       robot.right_arm.getNumberOfJoints() +
                       robot.left_leg.getNumberOfJoints()  +
                       robot.left_leg.getNumberOfJoints()  +
                       robot.torso.getNumberOfJoints()  +
                       robot.head.getNumberOfJoints()  +
                       robot.left_hand.getNumberOfJoints()  +
                       robot.right_hand.getNumberOfJoints()  ;
            }
            
            
            
            // pre declaration 
            yarp::sig::Matrix Homogeneous( const yarp::sig::Matrix& R_ab,  const yarp::sig::Vector& d_ab  );
            yarp::sig::Matrix iHomogeneous( const yarp::sig::Matrix& T_ab);
            yarp::sig::Matrix Adjoint_MT( const yarp::sig::Matrix& T_ab);
                 
	    yarp::sig::Vector sense_position_no_hands(RobotUtils& robot) {
                //yarp::sig::Vector wb_input_q = robot.sensePositionRefFeedback();
                yarp::sig::Vector wb_input_q = robot.sensePosition() ;
                yarp::sig::Vector input_q_no_hands(getNumberOfKinematicJoints(robot));

                for(int i=0;i<input_q_no_hands.size();i++) {
                    input_q_no_hands[i]=wb_input_q[i]; //hands not considered
                }
                return input_q_no_hands;
            }
            
            yarp::sig::Vector sense_torque_no_hands(RobotUtils& robot) {
                    yarp::sig::Vector wb_input_tau = robot.senseTorque();
                    yarp::sig::Vector input_tau_no_hands(getNumberOfKinematicJoints(robot));

                    for(int i=0;i<input_tau_no_hands.size();i++) {
                        input_tau_no_hands[i]=wb_input_tau[i]; //hands not considered
                    }
                    return input_tau_no_hands;
            }
     
            
            //------------------------------------------------------------
      /**
     * @brief getKq returns the joint stiffness matrix (hardcoded for now)
     * @return a yarp Matrix of dimension =  robot.getNumberOfJoints() x robot.getNumberOfJoints()
     */
    yarp::sig::Matrix getKq_Walkman(RobotUtils& robot ) {
                            yarp::sig::Vector Kq_vec_right_arm(  RIGHT_ARM_JOINT_NUM  ) ;
                            yarp::sig::Vector Kq_vec_left_arm(   LEFT_ARM_JOINT_NUM   ) ;  
                            yarp::sig::Vector Kq_vec_torso(      TORSO_JOINT_NUM      ) ;
                            yarp::sig::Vector Kq_vec_right_leg(  RIGHT_LEG_JOINT_NUM  ) ;
                            yarp::sig::Vector Kq_vec_left_leg(   LEFT_LEG_JOINT_NUM   ) ;
                            yarp::sig::Vector Kq_vec_head(       HEAD_JOINT_NUM       ) ;
                            
                            // RIGHT ARM    
                            Kq_vec_right_arm[0] = 1500.0 ; //1000.0 ;
                            Kq_vec_right_arm[1] = 1500.0 ; //1000.0 ;
                            Kq_vec_right_arm[2] = 1500.0 ; //600.0 ;
                            Kq_vec_right_arm[3] = 1500.0 ; //1000.0 ;
                            Kq_vec_right_arm[4] = 800.0 ; //100.0 ;
                            Kq_vec_right_arm[5] = 800.0 ; //100.0 ;
                            Kq_vec_right_arm[6] = 800.0 ; //10.0 ;
                            // LEFT ARM   
                            Kq_vec_left_arm[0] = 1500.0 ; //1000.0 ;
                            Kq_vec_left_arm[1] = 1500.0 ; //1000.0 ;
                            Kq_vec_left_arm[2] = 1500.0 ; //600.0 ;
                            Kq_vec_left_arm[3] = 1500.0 ; //1000.0 ;
                            Kq_vec_left_arm[4] = 800.0 ; //100.0 ;
                            Kq_vec_left_arm[5] = 800.0 ; //100.0 ;
                            Kq_vec_left_arm[6] = 800.0 ; //10.0 ;
                            //TORSO
                            Kq_vec_torso[0] = 4000.0 ; // 1000.0 ;
                            Kq_vec_torso[1] = 5500.0 ; // 1000.0 ;
                            Kq_vec_torso[2] = 4000.0 ; // 1000.0 ;
                            // RIGHT LEG
                            Kq_vec_right_leg[0] = 6500.0 ; //3000.0 ;
                            Kq_vec_right_leg[1] = 1063.0 ; //5000.0 ;
                            Kq_vec_right_leg[2] = 6500.0 ; //3000.0 ;
                            Kq_vec_right_leg[3] = 6500.0 ; //3000.0 ;
                            Kq_vec_right_leg[4] = 6500.0 ; //4000.0 ;
                            Kq_vec_right_leg[5] = 6500.0 ; //3000.0 ;
                            // LEFT LEG
                            Kq_vec_left_leg[0] = 6500.0 ; //3000.0 ;
                            Kq_vec_left_leg[1] = 1063.0 ; //5000.0 ;
                            Kq_vec_left_leg[2] = 6500.0 ; //3000.0 ;
                            Kq_vec_left_leg[3] = 6500.0 ; //3000.0 ;
                            Kq_vec_left_leg[4] = 6500.0 ; //4000.0 ;
                            Kq_vec_left_leg[5] = 6500.0 ; //3000.0 ;
                            // HEAD
                            Kq_vec_head[0] = 100.0 ;
                            Kq_vec_head[1] = 500.0 ;
                            //
                            yarp::sig::Vector Kq_vec( getNumberOfActuatedJoints(robot) )  ;     
                            robot.fromRobotToIdyn( Kq_vec_right_arm ,
                                                Kq_vec_left_arm  ,
                                                Kq_vec_torso  ,
                                                Kq_vec_right_leg ,
                                                Kq_vec_left_leg  ,
                                                Kq_vec_head,
                                                Kq_vec ); 
                            yarp::sig::Matrix Kq_matrix(  getNumberOfKinematicJoints(robot), getNumberOfKinematicJoints(robot) )  ; 
                            Kq_matrix.diagonal(  Kq_vec ) ;
                            return Kq_matrix ;
                            }        
    /**
     * @brief getKq returns the joint stiffness matrix (hardcoded for now)
     * @return a yarp Matrix of dimension =  robot.getNumberOfJoints() x robot.getNumberOfJoints()
     */
    yarp::sig::Matrix getKq(RobotUtils& robot ) {
                            yarp::sig::Vector Kq_vec_right_arm(  RIGHT_ARM_JOINT_NUM ) ;
                            yarp::sig::Vector Kq_vec_left_arm(   LEFT_ARM_JOINT_NUM ) ;  
                            yarp::sig::Vector Kq_vec_torso(      TORSO_JOINT_NUM ) ;
                            yarp::sig::Vector Kq_vec_right_leg(  RIGHT_LEG_JOINT_NUM ) ;
                            yarp::sig::Vector Kq_vec_left_leg(   LEFT_LEG_JOINT_NUM ) ;
                            // RIGHT ARM    
                            Kq_vec_right_arm[0] = 1000.0 ;
                            Kq_vec_right_arm[1] = 1000.0 ;
                            Kq_vec_right_arm[2] = 600.0 ;
                            Kq_vec_right_arm[3] = 1000.0 ;
                            Kq_vec_right_arm[4] = 100.0 ;
                            Kq_vec_right_arm[5] = 100.0 ;
                            Kq_vec_right_arm[6] = 10.0 ;
                            // LEFT ARM   
                            Kq_vec_left_arm[0] = 1000.0 ;
                            Kq_vec_left_arm[1] = 1000.0 ;
                            Kq_vec_left_arm[2] = 600.0 ;
                            Kq_vec_left_arm[3] = 1000.0 ;
                            Kq_vec_left_arm[4] = 100.0 ;
                            Kq_vec_left_arm[5] = 100.0 ;
                            Kq_vec_left_arm[6] = 10.0 ;
                            //TORSO
                            Kq_vec_torso[0] = 1000.0 ;
                            Kq_vec_torso[1] = 1000.0 ;
                            Kq_vec_torso[2] = 1000.0 ;
                            // RIGHT LEG
                            Kq_vec_right_leg[0] = 3000.0 ;
                            Kq_vec_right_leg[1] = 5000.0 ;
                            Kq_vec_right_leg[2] = 3000.0 ;
                            Kq_vec_right_leg[3] = 3000.0 ;
                            Kq_vec_right_leg[4] = 4000.0 ;
                            Kq_vec_right_leg[5] = 3000.0 ;
                            // LEFT LEG
                            Kq_vec_left_leg[0] = 3000.0 ;
                            Kq_vec_left_leg[1] = 5000.0 ;
                            Kq_vec_left_leg[2] = 3000.0 ;
                            Kq_vec_left_leg[3] = 3000.0 ;
                            Kq_vec_left_leg[4] = 4000.0 ;
                            Kq_vec_left_leg[5] = 3000.0 ;
                            //
                            yarp::sig::Vector Kq_vec( getNumberOfKinematicJoints(robot)  )  ;     
                            robot.fromRobotToIdyn( Kq_vec_right_arm ,
                                                Kq_vec_left_arm  ,
                                                Kq_vec_torso  ,
                                                Kq_vec_right_leg ,
                                                Kq_vec_left_leg  ,
                                                Kq_vec ); 
                            yarp::sig::Matrix Kq_matrix(  getNumberOfKinematicJoints(robot), getNumberOfKinematicJoints(robot) )  ; 
                            Kq_matrix.diagonal(  Kq_vec ) ;
                            return Kq_matrix ;
                            }

    
     /**
     * @brief senseMotorPosition we use this method to obtain the motor position
     * @return a yarp vector of dimension equal to robot.getNumberOfJoints() 
     */
    yarp::sig::Vector senseMotorPosition(RobotUtils& robot, bool flag_robot ) {
                                        yarp::sig::Vector q_motor( getNumberOfKinematicJoints(robot)  )  ; 
                                        if(flag_robot){
                                            yarp::sig::Vector q_link = locoman::utils::sense_position_no_hands(robot) ;
                                            return q_link ; 
                                        }
                                        else{
                                            yarp::sig::Vector q_link = locoman::utils::sense_position_no_hands(robot) ;
                                            yarp::sig::Vector tau    = locoman::utils::sense_torque_no_hands(robot) ;
                                            // 
                                            yarp::sig::Matrix Kq_matrix;
                                            if(robot.idynutils.getRobotName() == "coman") {
                                                Kq_matrix = getKq(robot) ;
                                            }
                                            else {
                                                Kq_matrix = getKq_Walkman(robot) ;
                                            }
                                            
                                           // std::cout << Kq_matrix.toString() << std::endl;
                                            yarp::sig::Matrix Cq_matrix = yarp::math::luinv(Kq_matrix) ;
                                            q_motor = Cq_matrix*tau  + q_link ;
                                            return q_motor ;
					    }
                                        }
        
    //-----------------------------------------------------------
    //-----------------------------------------------------------------
     /**
     * @brief  fConToSens maps contact forces on the sensor frame, is equivalent to a grasp matrix
     * @param  sens_index, ... c4_index are the indexes of the sensor and of the 4 contact points
     * @return is a 6x6 yarp matrix, able to map contact forces in sensor wrenches
     */
    yarp::sig::Matrix fConToSens( const int sens_index,
                                  const int c1_index,
                                  const int c2_index,
                                  const int c3_index,
                                  const int c4_index,
                                  iDynUtils& model 
                                 // bool flag_robot
                                ) 
    {
    
      yarp::sig::Matrix map_fConToSens( 6 , 12 ) ;
                                        yarp::sig::Matrix B_select( 6 , 3 ) ;
                                        yarp::sig::Matrix T_w_sensor( 4 , 4 ) ;
                                        yarp::sig::Matrix T_w_c1( 4 , 4 ) ;
                                        yarp::sig::Matrix T_w_c2( 4 , 4 ) ;
                                        yarp::sig::Matrix T_w_c3( 4 , 4 ) ;
                                        yarp::sig::Matrix T_w_c4( 4 , 4 ) ;
                                        yarp::sig::Matrix Eye_3( 3 , 3 )  ;
                                        Eye_3.eye() ;
                                        // 
                                        B_select.setSubmatrix( Eye_3 , 0 , 0 ) ;
                                        /*if(flag_robot)
                                        {
                                        yarp::sig::Matrix T_w_ankle = model.iDyn3_model.getPosition(sens_index) ;
                                        yarp::sig::Vector d_sens_ankle (3, 0.0) ;
                                        d_sens_ankle[2] = 0.07 ; //   // measured 7 cm
                                        yarp::sig::Matrix T_sensor_ankle = locoman::utils::Homogeneous(Eye_3, d_sens_ankle) ;
                                        //
                                        yarp::sig::Matrix T_ankle_w = locoman::utils::iHomogeneous(T_w_ankle) ;
                                        yarp::sig::Matrix T_sensor_w = T_sensor_ankle * T_ankle_w ;
                                        T_w_sensor =  locoman::utils::iHomogeneous(T_sensor_w) ;
                                        //  
                                        }*/ 
                                        // else
                                        // {
                                        T_w_sensor = model.iDyn3_model.getPosition(sens_index) ;     
                                        //}
                                        //
                                        T_w_c1  = model.iDyn3_model.getPosition(c1_index)  ;
                                        T_w_c2  = model.iDyn3_model.getPosition(c2_index)  ;    
                                        T_w_c3  = model.iDyn3_model.getPosition(c3_index)  ;    
                                        T_w_c4  = model.iDyn3_model.getPosition(c4_index)  ;    
                                        
                                        yarp::sig::Matrix Ad_1 = locoman::utils::Adjoint_MT( locoman::utils::iHomogeneous(T_w_sensor)*T_w_c1  ) *B_select  ;
                                        yarp::sig::Matrix Ad_2 = locoman::utils::Adjoint_MT( locoman::utils::iHomogeneous(T_w_sensor)*T_w_c2  ) *B_select  ;
                                        yarp::sig::Matrix Ad_3 = locoman::utils::Adjoint_MT( locoman::utils::iHomogeneous(T_w_sensor)*T_w_c3  ) *B_select  ;
                                        yarp::sig::Matrix Ad_4 = locoman::utils::Adjoint_MT( locoman::utils::iHomogeneous(T_w_sensor)*T_w_c4  ) *B_select  ;

                                        map_fConToSens.setSubmatrix( Ad_1 , 0 , 0 ) ;
                                        map_fConToSens.setSubmatrix( Ad_2 , 0 , 3 ) ;
                                        map_fConToSens.setSubmatrix( Ad_3 , 0 , 6 ) ;
                                        map_fConToSens.setSubmatrix( Ad_4 , 0 , 9 ) ;

                                        return map_fConToSens ;
                                }

     /**
     * @brief The function computes and returns the pose of the frame auxiliary world frame {AW} with 
     * @return respect to the world
     */
     yarp::sig::Matrix AW_world_posture( iDynUtils& model , 
                                         RobotUtils& robot  ) {
                                                    yarp::sig::Vector zero_3(3, 0.0) ;
                                                    int imu_link_index = model.iDyn3_model.getLinkIndex("imu_link") ; 
                                                    yarp::sig::Matrix T_w_imu_0 = model.iDyn3_model.getPosition( imu_link_index ) ;    
                                                    yarp::sig::Matrix T_imu_w_0 = locoman::utils::iHomogeneous(T_w_imu_0) ; 
                                                    RobotUtils::IMUPtr IMU_ptr = robot.getIMU()  ;
                                                    yarp::sig::Vector IMU_sense = IMU_ptr->sense(); ;
                                                    yarp::sig::Vector IMU_sense_lin_acc(3) ; 
                                                            
                                                    IMU_sense_lin_acc[0] = IMU_sense[3] ;    
                                                    IMU_sense_lin_acc[1] = IMU_sense[4] ;    
                                                    IMU_sense_lin_acc[2] = IMU_sense[5] ;    
                                                        
                                                    double norm_imu = norm(IMU_sense_lin_acc)     ;
                                                    //
                                                    // Defining z-axis for the Auxiliary World => frame {AW} used as world in loop of the run function
                                                    yarp::sig::Vector z_imu_aw = -1.0* IMU_sense_lin_acc/norm_imu ;   // Z axis vertical, going up

                                                    double norm_z = norm(z_imu_aw)  ;
                                                    
                                                    yarp::sig::Matrix z_imu_aw_matr(3,1);
                                                    z_imu_aw_matr[0][0] = z_imu_aw[0] ;
                                                    z_imu_aw_matr[1][0] = z_imu_aw[1] ;
                                                    z_imu_aw_matr[2][0] = z_imu_aw[2] ;   

                                                    int waist_index = model.iDyn3_model.getLinkIndex("Waist") ;
                                                        
                                                    yarp::sig::Matrix T_w_waist_0 = model.iDyn3_model.getPosition( waist_index )  ;   
                                                    //yarp::sig::Matrix T_waist_w_0 = locoman::utils::iHomogeneous(T_w_waist_0) ;
                                                    yarp::sig::Matrix T_imu_waist_0 = T_imu_w_0* T_w_waist_0  ;
                                                        
                                                    yarp::sig::Vector x_imu_waist(3) ;
                                                    x_imu_waist[0] =  T_imu_waist_0[0][0];
                                                    x_imu_waist[1] =  T_imu_waist_0[1][0];
                                                    x_imu_waist[2] =  T_imu_waist_0[2][0];

                                                    yarp::sig::Matrix Null_z_tr =  nullspaceProjection(z_imu_aw_matr.transposed()) ;
                                                    yarp::sig::Vector x_imu_aw = Null_z_tr*x_imu_waist  ;

                                                    yarp::sig::Vector y_imu_aw(3) ;
                                                    yarp::sig::Matrix R_imu_aw_0(3,3) ;    
                                                        
                                                    if (norm(x_imu_aw)>0.01)
                                                    {
                                                    x_imu_aw = x_imu_aw/norm(x_imu_aw) ;   
                                                    R_imu_aw_0[1][0] = x_imu_aw[1] ;
                                                    R_imu_aw_0[2][0] = x_imu_aw[2] ;

                                                    R_imu_aw_0[0][1] = y_imu_aw[0] ;
                                                    R_imu_aw_0[1][1] = y_imu_aw[1] ;
                                                    R_imu_aw_0[2][1] = y_imu_aw[2] ;

                                                    R_imu_aw_0[0][2] = z_imu_aw[0] ;
                                                    R_imu_aw_0[1][2] = z_imu_aw[1] ;
                                                    R_imu_aw_0[2][2] = z_imu_aw[2] ; 
                                                    }
                                                    else   {   
                                                    x_imu_aw[0] = Null_z_tr[0][0] ; 
                                                    x_imu_aw[1] = Null_z_tr[1][0] ; 
                                                    x_imu_aw[2] = Null_z_tr[2][0] ;        
                                                    x_imu_aw = x_imu_aw/norm(x_imu_aw) ;     
                                                    }  
                                                    y_imu_aw = cross(z_imu_aw, x_imu_aw  );   
                                                    
                                                    R_imu_aw_0[0][0] = x_imu_aw[0] ;
                                                    R_imu_aw_0[1][0] = x_imu_aw[1] ;
                                                    R_imu_aw_0[2][0] = x_imu_aw[2] ;

                                                    R_imu_aw_0[0][1] = y_imu_aw[0] ;
                                                    R_imu_aw_0[1][1] = y_imu_aw[1] ;
                                                    R_imu_aw_0[2][1] = y_imu_aw[2] ;

                                                    R_imu_aw_0[0][2] = z_imu_aw[0] ;
                                                    R_imu_aw_0[1][2] = z_imu_aw[1] ;
                                                    R_imu_aw_0[2][2] = z_imu_aw[2] ;
                                                        
                                                    // Origin of {AW} coincident with the origin of {IMU}
                                                    yarp::sig::Matrix T_imu_aw_0 = locoman::utils::Homogeneous( R_imu_aw_0, zero_3 ) ;
                                                    // yarp::sig::Matrix T_aw_imu_0 = locoman::utils::iHomogeneous(T_imu_aw_0) ;    
                                                    yarp::sig::Matrix T_w_aw = T_w_imu_0 * T_imu_aw_0 ;
                                                    // yarp::sig::Matrix T_aw_w_0 = locoman::utils::iHomogeneous(T_w_aw_0) ;  
                                                    return T_w_aw ;  
    }  
        
                          
                                
//    /*  /**
//      * @brief  sigma_frict computes the metrics measuring the goodness of the contact force with respect to friction limits 
//      * @param  fc is the contact force. The normal is assumed to be n = [0 0 1]^T 
//      * @param  mu is the friction coefficient
//      * @return the value of sigma_frict
//      */
//       double sigma_frict( const yarp::sig::Vector& fc ,
//                           const double mu 
//                                       ) {
//                                         yarp::sig::Vector fc_aux = fc ;
//                                         if(fc_aux(2)<0){
//                                             fc_aux= -1.0*fc; 
//                                         } ;
//                                         yarp::sig::Vector normal(3);
//                                         normal(0) = 0 ;
//                                         normal(1) = 0 ;
//                                         normal(2) = 1 ;
//                                         //
//                                         double alpha_frict = 1/(sqrt( 1+ pow(mu,2))) ;
//                                         double sigma_frict = alpha_frict*yarp::math::norm(fc_aux) - yarp::math::dot(fc_aux,normal) ;
//                                         return sigma_frict ;
//                                                                             } 
//         /**                                      
//         * @brief  sigma_min computes the distance (along a certain metric) with respect to minimum force allowed 
//         * @param  fc_i    is the contact force.  
//         * @param  f_min_i is the minimum module of the force allowed
//         * @return the value of sigma_min
//         */
//             double sigma_min( const yarp::sig::Vector& fc_i,
//                                 const double f_min_i 
//                                             ) {
//                                         yarp::sig::Vector fc_aux = fc_i ;
//                                         if(fc_aux(2)<0){
//                                             fc_aux= -1.0*fc_i; 
//                                         } ;
//                                         yarp::sig::Vector normal(3) ;
//                                         normal(0) = 0 ;
//                                         normal(1) = 0 ;
//                                         normal(2) = 1 ;
//                                         //
//                                         double sigma_min =  f_min_i - yarp::math::dot(fc_aux,normal) ;
//                                         return sigma_min;
//                                     } 
                                      
//                                            /**
//      * @brief  sigma_max computes the distance (along a certain metric) with respect to maximum force allowed 
//      * @param  fc_i is the contact force.  
//      * @param  f_max is the maximum module of the force allowed
//      * @return the value of sigma_min
//      */
//       double sigma_max( const yarp::sig::Vector& fc_i,
//                         const double f_max 
//                                       ) {
//                                             yarp::sig::Vector fc_aux = fc_i ;
//                                             if(fc_aux(2)<0){
//                                                 fc_aux= -1.0*fc_i; 
//                                             } ;
//                                             yarp::sig::Vector normal(3) ;
//                                             normal(0) = 0 ;
//                                             normal(1) = 0 ;
//                                             normal(2) = 1 ;
//                                             //
//                                             double sigma_max = - f_max +   yarp::math::dot(fc_aux,normal) ;
//                                             return sigma_max ;
//                                     }*/ 
        /**                                      
        * @brief  sigma_i computes the sigma values for friction, minimum and maximum force constraints
        * @param  fc_i is the contact force vector to the i-th contact point (3 elements, only linear components of the force are admitted)
	* @param  n_i  is the direction of the normal to the contact in local coordinates at the i-th contact point (3 elements)
	* @param  mu_i is the friction coefficient at the i-th contact point
        * @param  f_min_i is the minimum module of the normal force allowed
	* @param  f_max_i is the maximum module of the normal force allowed
        * @return sigma_i: 3-elemnts vector containing the sigma value for friction, minimum and maximum force constrints
        */
        yarp::sig::Vector sigma_i( const yarp::sig::Vector& fc_i, 
				   const yarp::sig::Vector& n_i, 
				   const double mu_i, 
                                   const double f_min_i,
				   const double f_max_i 
                                            ) {
// 	  		          yarp::sig::Matrix f_c_matrix_row(1,3) ;
// 				  f_c_matrix_row[0][0] = fc_i(0) ;
// 				  f_c_matrix_row[0][1] = fc_i(1) ;
// 				  f_c_matrix_row[0][2] = fc_i(2) ;
	  
				  yarp::sig::Vector alpha(3,0.0) ;
				  yarp::sig::Vector beta(3,0.0) ;
				  yarp::sig::Vector gamma(3,0.0) ;
				  alpha(0) = 1.0/(std::sqrt(1.0 + std::pow(mu_i,2.0) )) ;
				  beta(0) = -1.0 ;
				  beta(1) = -1.0 ;
				  beta(2) =  1.0 ;
				  gamma(1) =  f_min_i ;
				  gamma(2) = -f_max_i ;
				  yarp::sig::Vector sigma_i_vect(3,0.0) ; 				  
				  for(int i = 0; i <sigma_i_vect.length() ; i++){
				    sigma_i_vect(i) = alpha(i)*norm(fc_i) + beta(i)
				      *yarp::math::dot(fc_i,n_i)  + gamma(i) ;
				  }
				  return sigma_i_vect ;
                                    }                                
        /**                              
        * @brief  sigma_tot computes the distance (along a certain metric) with respect to minimum force allowed 
        * @param  fc is the vector collecting all the contact forces (3 elements each)
	* @param  normals is the vector collecting all normal vector to the contact point (3 elements each)
	* @param  mu is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @param  f_min is the vector collecting all the friction coefficient of each contact point (1 element each)
	* @param  f_max is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @return sigma_tot: vector collecting all the sigma_ij (=> computed for each contact and each constraint)
        */
        yarp::sig::Vector sigma_tot(  const yarp::sig::Vector& fc, 
				      const yarp::sig::Vector& normals, 
				      const yarp::sig::Vector& mu, 
				      const yarp::sig::Vector& f_min,
				      const yarp::sig::Vector& f_max 
                                            ) {
				  int n_forces = mu.length() ; // 1 mu for each contact point
				  
				  yarp::sig::Vector start_i(n_forces,0.0) ;
				  yarp::sig::Vector end_i(n_forces,0.0) ;
				  //start_i(0) = 0 ;
				  for(int i=1 ; i<n_forces;i++){
				    start_i(i) = start_i(i-1) + 3.0 ;
				  }
	  			  for(int i=0 ; i<n_forces;i++){
				    end_i(i) = start_i(i) + 2.0 ;
				  }
				  yarp::sig::Vector sigma_vect(fc.length(),0.0); // 3 elements of force each contact, 3 sigma each contact
				  for(int i=0 ; i<n_forces;i++){
				  yarp::sig::Vector fc_loop = fc.subVector(start_i(i), end_i(i));
				  yarp::sig::Vector n_loop = normals.subVector(start_i(i), end_i(i)) ;
				  double mu_loop = mu(i);
				  double fmin_loop = f_min(i);
				  double fmax_loop = f_max(i);
				    
				 std::cout << " fc_loop =  "<< std::endl << fc_loop.toString() << std::endl  ; 
                                 std::cout << " n_loop =  "<< std::endl << n_loop.toString() << std::endl  ; 
				 std::cout << " mu_loop =  "<< std::endl << mu_loop << std::endl  ; 
				 std::cout << " fmin_loop =  "<< std::endl << fmin_loop << std::endl  ; 
				 std::cout << " fmax_loop =  "<< std::endl << fmax_loop << std::endl  ; 

				 
				    yarp::sig::Vector sigma_vect_i = sigma_i( fc_loop, 
								     n_loop,
								     mu_loop, fmin_loop,fmax_loop )  ;
				  sigma_vect.setSubvector(start_i(i) , sigma_vect_i  )  ;
				  }				  
				  return sigma_vect ;
                                    }                                         
	
	/**                                      
        * @brief  sigma_min computes the distance (along a certain metric) with respect to minimum force allowed 
        * @param  fc_i is the contact force vector to the i-th contact point (3 elements, only linear components of the force are admitted)
	* @param  n_i  is the direction of the normal to the contact in local coordinates at the i-th contact point (3 elements)
	* @param  mu_i is the friction coefficient at the i-th contact point
        * @param  f_min_i is the minimum module of the normal force allowed
	* @param  f_max_i is the maximum module of the normal force allowed
        * @return sigma_i: 3-elemnts vector containing the sigma value for friction, minimum and maximum force constrints
        */
        double V_i( const yarp::sig::Vector& fc_i, 
		    const yarp::sig::Vector& n_i, 
		    const double mu_i, 
		    const double f_min_i,
		    const double f_max_i 
			      ) {
		    yarp::sig::Vector sigma_contact = sigma_i( fc_i,  n_i,  mu_i, f_min_i, f_max_i  ) ;
		    double epsilon = std::pow(std::numeric_limits<double>::epsilon(), 1.0/8.0);
		    double a = (3.0/2.0)*(1.0/(std::pow(epsilon,4.0)));
		    double b = 4.0*(1.0/(std::pow(epsilon,3.0))) ;
		    double c = 3.0/(std::pow(epsilon,2.0) ) ;
		    yarp::sig::Vector V_i_vect(3,0.0) ;
		    for(int i = 0; i < sigma_contact.length(); i++){
		      if(sigma_contact(i)<epsilon){
		      V_i_vect(i) = 1.0/(2.0*std::pow(sigma_contact(i),2.0)) ;
		      }
		      else{
		      V_i_vect(i) = a*std::pow(sigma_contact(i),2) + b*sigma_contact(i) + c ;
		      } 
		    }
		    yarp::sig::Vector ones_1_3(3,1.0) ;
// 		    ones_1_3[0][0] = 1.0 ;
// 		    ones_1_3[0][1] = 1.0 ;
// 		    ones_1_3[0][2] = 1.0 ;
		    double V_sum_i = yarp::math::dot(ones_1_3, V_i_vect ); 
		    return V_sum_i;
		      }                                         
       
        /**                              
        * @brief  V_tot computes the value of the function V for the actual contact force distribution
        * @param  fc is the vector collecting all the contact forces (3 elements each)
	* @param  normals is the vector collecting all normal vector to the contact point (3 elements each)
	* @param  mu is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @param  f_min is the vector collecting all the friction coefficient of each contact point (1 element each)
	* @param  f_max is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @return V_final: the value of the function V
        */
        double V_tot( const yarp::sig::Vector& fc, 
		      const yarp::sig::Vector& normals, 
		      const yarp::sig::Vector& mu, 
		      const yarp::sig::Vector& f_min,
		      const yarp::sig::Vector& f_max 
			    ) {
		      int n_forces = mu.length() ; // 1 mu for each contact point
		      yarp::sig::Vector start_i(n_forces, 0.0) ;
		      yarp::sig::Vector end_i(n_forces,  0.0) ;
		      //start_i(0) = 0 ;
		      for(int i=1 ; i<n_forces;i++){
			start_i(i) = start_i(i-1) + 3 ;
		      }
		      for(int i=0 ; i<n_forces;i++){
			end_i(i) = start_i(i) + 2;
		      }
		      //
		      double V_final  = 0.0 ; 
		      for(int i=0 ; i<n_forces;i++){
			 V_final += V_i(fc.subVector(start_i(i), end_i(i)) , 
					normals.subVector(start_i(i), end_i(i)),
					mu(i), f_min(i),f_max(i) )  ;
			 }				  
		      return V_final ; 
		      }    
       
        /**                                      
        * @brief  D_V_i computes the derivative of V_i wrt y
        * @param  fc_i is the contact force vector to the i-th contact point (3 elements, only linear components of the force are admitted)
	* @param  n_i  is the direction of the normal to the contact in local coordinates at the i-th contact point (3 elements)
	* @param  mu_i is the friction coefficient at the i-th contact point
        * @param  f_min_i is the minimum module of the normal force allowed
	* @param  f_max_i is the maximum module of the normal force allowed
	* @param  E_i is the portion of the bassis of the controllable contact forces relative to the the i-th contact force
        * @return D_V_i_vect: the derivative of the V_i
        */
        yarp::sig::Vector D_V_i( const yarp::sig::Vector& fc_i, 
		    const yarp::sig::Vector& n_i, 
		    const double mu_i, 
		    const double f_min_i,
		    const double f_max_i,
		    const yarp::sig::Matrix& E_i 
			      ) {
		    yarp::sig::Vector sigma_contact = sigma_i( fc_i,  n_i,  mu_i, f_min_i, f_max_i  ) ;
		    yarp::sig::Vector alpha(3,0.0) ;
		    yarp::sig::Vector beta(3,0.0) ;
		    //yarp::sig::Vector gamma(3,0.0) ;
		    alpha(0) = 1.0/(std::sqrt(1.0+ std::pow(mu_i,2.0))) ;
		    beta(0) = -1.0 ;
		    beta(1) = -1.0 ;
		    beta(2) =  1.0 ;
		    //gamma(1) =  f_min_i ;
		    //gamma(2) = -f_max_i ;
		    double epsilon = std::pow(std::numeric_limits<double>::epsilon(), (1.0/8.0) ) ;
		    double a = (3.0/2.0)*(1.0/(std::pow(epsilon,4.0)));
		    double b =  4.0*(1.0/(std::pow(epsilon,3.0))) ;
		    //double c = 3.0/(std::pow(epsilon,2.0) ;
		    
// 		    yarp::sig::Matrix f_c_matrix_row(1,3) ;
// 		    f_c_matrix_row[0][0] = fc_i(0) ;
// 		    f_c_matrix_row[0][1] = fc_i(1) ;
// 		    f_c_matrix_row[0][2] = fc_i(2) ;
// 		    yarp::sig::Matrix f_c_matrix_col = f_c_matrix_row.transposed() ;
// 		    
		    int n_y_elements = E_i.cols() ;
		    yarp::sig::Vector D_V_i_vect( n_y_elements ,0.0) ;
		    
		    for(int i = 0; i < sigma_contact.length()-1; i++){
		      if(sigma_contact(i)<epsilon){
			  D_V_i_vect += -1.0 * ( pow(sigma_contact(i),-3.0))*(
					  alpha(i)* pow( yarp::math::dot(fc_i,fc_i) , -1.0/2.0)* E_i.transposed()*fc_i +
					  beta(i)* E_i.transposed()*n_i 
			  );
			  //1.0/(2.0*pow(sigma_contact(i),2.0)) ;
		      }
		      else{
			  D_V_i_vect += (2.0*a*sigma_contact(i) + b)*(
			    alpha(i)* pow( yarp::math::dot(fc_i,fc_i) , -1.0/2.0)* E_i.transposed()*fc_i +
					  beta(i)* E_i.transposed()*n_i 
			  )   ;
		      } 
		    }
		    return D_V_i_vect ;
		      }  
       
        /**                              
        * @brief  D_V_tot computes the value of the derivative of V for the actual contact force distribution
        * @param  fc is the vector collecting all the contact forces (3 elements each)
	* @param  normals is the vector collecting all normal vector to the contact point (3 elements each)
	* @param  mu is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @param  f_min is the vector collecting all the friction coefficient of each contact point (1 element each)
	* @param  f_max is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @return D_V_final: the value of the derivative of V wrt the y coefficients
        */
        yarp::sig::Vector D_V_tot( const yarp::sig::Vector& fc, 
		      const yarp::sig::Vector& normals, 
		      const yarp::sig::Vector& mu, 
		      const yarp::sig::Vector& f_min,
		      const yarp::sig::Vector& f_max,
		      const yarp::sig::Matrix& E 
			    ) {
		      int n_forces = mu.length() ; // 1 mu for each contact point
		      yarp::sig::Vector start_i(n_forces,0.0) ;
		      yarp::sig::Vector end_i(n_forces,0.0) ;
		      //start_i(0) = 0 ;
		      for(int i=1 ; i<n_forces;i++){
			start_i(i) = start_i(i-1) + 3 ;
		      }
		      for(int i=0 ; i<n_forces;i++){
			end_i(i) = start_i(i) + 2;
		      }      
		      int n_y_elements = E.cols() ;
		      yarp::sig::Vector D_V_final( n_y_elements ,0.0) ;
		      for(int i=0 ; i<n_forces;i++){
			 D_V_final += D_V_i(fc.subVector(start_i(i), end_i(i)) , 
					normals.subVector(start_i(i), end_i(i)),
					mu(i), f_min(i),f_max(i),
					E.submatrix(start_i(i), end_i(i) , 0,n_y_elements-1 )) ;
			 }
		       return D_V_final ; 
		      }           
       
       
         /**                                      
        * @brief  H_V_i computes the Hessian matrix of V_i 
        * @param  fc_i is the contact force vector to the i-th contact point (3 elements, only linear components of the force are admitted)
	* @param  n_i  is the direction of the normal to the contact in local coordinates at the i-th contact point (3 elements)
	* @param  mu_i is the friction coefficient at the i-th contact point
        * @param  f_min_i is the minimum module of the normal force allowed
	* @param  f_max_i is the maximum module of the normal force allowed
	* @param  E_i is the portion of the bassis of the controllable contact forces relative to the the i-th contact force
        * @return H_V_i_matr: the Hessian of V_i
        */
        yarp::sig::Matrix H_V_i( const yarp::sig::Vector& fc_i, 
		    const yarp::sig::Vector& n_i, 
		    const double mu_i, 
		    const double f_min_i,
		    const double f_max_i,
		    const yarp::sig::Matrix& E_i 
			      ) {
		    yarp::sig::Vector sigma_contact = sigma_i( fc_i,  n_i,  mu_i, f_min_i, f_max_i  ) ;
		    yarp::sig::Vector alpha(3,0.0) ;
		    yarp::sig::Vector beta(3,0.0) ;
		    //yarp::sig::Vector gamma(3,0.0) ;
		    alpha(0) = 1.0/(std::sqrt(1.0+ std::pow(mu_i,2.0))) ;
		    beta(0) = -1.0 ;
		    beta(1) = -1.0 ;
		    beta(2) =  1.0 ;
		    //gamma(1) =  f_min_i ;
		    //gamma(2) = -f_max_i ;
		    double epsilon = std::pow(std::numeric_limits<double>::epsilon(), (1.0/8.0));
		    double a = (3.0/2.0)*(1.0/(std::pow(epsilon,4.0)));
		    double b = 4.0*(1.0/(std::pow(epsilon,3.0))) ;
		    //double c = 3.0/(std::pow(epsilon,2.0) ;
		    		    
		    int n_y_elements = E_i.cols() ;
		    yarp::sig::Matrix H_V_i_matrix( n_y_elements ,n_y_elements) ;
		    H_V_i_matrix.zero(); 
		    //
		    for(int i = 0; i < sigma_contact.length()-1; i++){
		      yarp::sig::Vector d_sig_y = (alpha(i)* pow( yarp::math::dot(fc_i,fc_i) , -1.0/2.0) * E_i.transposed()*fc_i +
					  beta(i)* E_i.transposed()*n_i );
		      yarp::sig::Matrix eye_rows_Ei(E_i.rows(), E_i.rows()  ) ;
		      eye_rows_Ei.eye() ;
		      yarp::sig::Matrix d_sig_y_2 = alpha(i)* pow( yarp::math::dot(fc_i,fc_i) , -1.0/2.0)*E_i.transposed() *
					  (eye_rows_Ei -
					  (yarp::math::outerProduct(fc_i,fc_i)/yarp::math::dot(fc_i,fc_i))  
					  )*E_i;
		      
		      if(sigma_contact(i)<epsilon){
			  H_V_i_matrix += -1.0*( pow(sigma_contact(i),-3.0))*d_sig_y_2 
						+ 3.0*( pow(sigma_contact(i),-4.0))
						  *yarp::math::outerProduct(d_sig_y,d_sig_y);
		      }
		      else{
			  H_V_i_matrix += (2.0*a*sigma_contact(i) + b)*d_sig_y_2+
					    2.0*a*yarp::math::outerProduct(d_sig_y,d_sig_y);
			  /*(
			    alpha(i)* pow( dot(fc_i,fc_i) , -1/2)* E_i.transposed*fc_i +
					  beta(i)* E_i.transposed*n_i 
			  )   ;*/
		      } 
		    }
		    return H_V_i_matrix ;
		    }  
       
	/**                              
        * @brief  H_V_tot computes the Hessian of V for the actual contact force distribution
        * @param  fc is the vector collecting all the contact forces (3 elements each)
	* @param  normals is the vector collecting all normal vector to the contact point (3 elements each)
	* @param  mu is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @param  f_min is the vector collecting all the friction coefficient of each contact point (1 element each)
	* @param  f_max is the vector collecting all the friction coefficient of each contact point (1 element each)
        * @return H_V_final: the Hessian of V at the actual contact force distribution
        */
        yarp::sig::Matrix H_V_tot( const yarp::sig::Vector& fc, 
		      const yarp::sig::Vector& normals, 
		      const yarp::sig::Vector& mu, 
		      const yarp::sig::Vector& f_min,
		      const yarp::sig::Vector& f_max,
		      const yarp::sig::Matrix& E 
			    ) {
		      int n_forces = mu.length() ; // 1 mu for each contact point
		      yarp::sig::Vector start_i(n_forces,0.0) ;
		      yarp::sig::Vector end_i(n_forces,0.0) ;
		      for(int i=1 ; i<n_forces;i++){
			start_i(i) = start_i(i-1) + 3 ;
		      }
		      for(int i=0 ; i<n_forces;i++){
			end_i(i) = start_i(i) + 2;
		      }
		      int n_y_elements = E.cols() ;
		      yarp::sig::Matrix H_V_final( n_y_elements , n_y_elements ) ;
		      for(int i=0 ; i<n_forces;i++){
			 H_V_final += H_V_i( fc.subVector(start_i(i), end_i(i)) , 
					normals.subVector(start_i(i), end_i(i)) ,
					mu(i), f_min(i),f_max(i),
					E.submatrix(start_i(i), end_i(i) , 0,n_y_elements-1 )) ;
			 }
		       return H_V_final ; 
		      }              
              
                                      
                                      
//      /**
//      * @brief  V_ij computes the distance (along a certain metric) with respect to the contact limits
//      * @param  sigma is one exit of the functions sigma_frict, sigma_min, sigma_max
//      * @param  toll is the admitted tolerance with respect to sigma limit value (= 0)
//      * @return the value of V_ij
//      */
//       double V_ij( const double sigma, 
//                    const double toll  = 1E-7
//                  ) {
//                                     double V_ij ;
//                                     if(sigma<toll){
//                                     V_ij = 1/(2*pow(sigma,2)) ;
//                                     }
//                                     else{
//                                     double a = 3/(2*pow(toll,4)) ;
//                                     double b = 4/(  pow(toll,3)) ;
//                                     double c = 3/(  pow(toll,2)) ;
//                                     V_ij = a*pow(sigma,2) + b*sigma + c ;
//                                     }
//                                     return V_ij ;
                                                                       // } 
       /**
        * @brief  provide the initial configuration
        * @return  yarp vector 
        */
        yarp::sig::Vector q_init(  RobotUtils& robot,
                                   yarp::sig::Vector& right_arm_configuration,
                                   yarp::sig::Vector& left_arm_configuration, 
                                   yarp::sig::Vector& torso_configuration, 
                                   yarp::sig::Vector& right_leg_configuration,
                                   yarp::sig::Vector& left_leg_configuration 
                                ) {
                                        yarp::sig::Vector q_motor_init(getNumberOfKinematicJoints(robot) ) ;                

                                        robot.fromRobotToIdyn(  right_arm_configuration ,
                                                                left_arm_configuration  ,
                                                                torso_configuration     ,
                                                                right_leg_configuration ,
                                                                left_leg_configuration  ,
                                                                q_motor_init     );      
                                        
                                        return q_motor_init ; 
            
                                    } 
     
     /**
     * @brief  provide the desired contact force distribution => 100% on the right
     * @return  0 if ok 
     */
     int FC_DES_right( yarp::sig::Vector& FC_DES, double mg  ) {
                        double part = -10.0/10.0 ;  // - => moving on the right; + => moving on the left
                        // On the left foot
                        FC_DES[2] = - (mg/8.0 + part*(mg/8.0) ) ;
                        FC_DES[5] = - (mg/8.0 + part*(mg/8.0) )   ;
                        FC_DES[8] = - (mg/8.0 + part*(mg/8.0) )   ;
                        FC_DES[11] = - (mg/8.0 + part*(mg/8.0) )   ;
                        // On the right foot
                        FC_DES[14] = - (mg/8.0 - part*(mg/8.0) )  ;
                        FC_DES[17] = - (mg/8.0 - part*(mg/8.0) )   ;
                        FC_DES[20] = - (mg/8.0 - part*(mg/8.0) )   ;
                        FC_DES[23] = - (mg/8.0 - part*(mg/8.0) )   ; 
                        //     
                        return 0 ; // FC_DES ;
                        }  
                        
        /**
        * @brief  provide the desired contact force distribution => 100% on the left
        * @return  0 if ok
        */
        int FC_DES_left( yarp::sig::Vector& FC_DES, double mg   ) {
                        double part = 10.0/10.0 ;  // - => moving on the right; + => moving on the left
                        // On the left foot
                        FC_DES[2] = - (mg/8.0 + part*(mg/8.0) ) ;
                        FC_DES[5] = - (mg/8.0 + part*(mg/8.0) )   ;
                        FC_DES[8] = - (mg/8.0 + part*(mg/8.0) )   ;
                        FC_DES[11] = - (mg/8.0 + part*(mg/8.0) )   ;
                        // On the right foot
                        FC_DES[14] = - (mg/8.0 - part*(mg/8.0) )  ;
                        FC_DES[17] = - (mg/8.0 - part*(mg/8.0) )   ;
                        FC_DES[20] = - (mg/8.0 - part*(mg/8.0) )   ;
                        FC_DES[23] = - (mg/8.0 - part*(mg/8.0) )   ; 
                        //     
                        return 0 ; // FC_DES ;
    }   

     /**
     * @brief  provide the desired contact force distribution => 50% on the left/right
     * @return  0 if ok
     */
     int FC_DES_center( yarp::sig::Vector& FC_DES, double mg ) {
                            double part = 0.0/10.0 ;  // - => moving on the right; + => moving on the left
                            // On the left foot
                            FC_DES[2] = - (mg/8.0 + part*(mg/8.0) ) ;
                            FC_DES[5] = - (mg/8.0 + part*(mg/8.0) )   ;
                            FC_DES[8] = - (mg/8.0 + part*(mg/8.0) )   ;
                            FC_DES[11] = - (mg/8.0 + part*(mg/8.0) )   ;
                            // On the right foot
                            FC_DES[14] = - (mg/8.0 - part*(mg/8.0) )  ;
                            FC_DES[17] = - (mg/8.0 - part*(mg/8.0) )   ;
                            FC_DES[20] = - (mg/8.0 - part*(mg/8.0) )   ;
                            FC_DES[23] = - (mg/8.0 - part*(mg/8.0) )   ; 
                            //     
                            return 0 ; // FC_DES ;
    }  
     
    /**
    * @brief  provide the desired contact force distribution => 0% on the left/right
    * @param  front_part fraction of the weight desired on the front part of the feet, it has to be <1.0
    * @return  0 if ok
    */
    int FC_DES_front( yarp::sig::Vector& FC_DES, double mg , double front_part) {
			  double rear_part = (1.0-front_part) ; 
			  // On the left foot
			  FC_DES[2] =  - ( front_part *(mg/8.0) )   ;
			  FC_DES[5] =  - ( front_part *(mg/8.0) )   ;
			  FC_DES[8] =  - ( rear_part  *(mg/8.0) )   ;
			  FC_DES[11] = - ( rear_part  *(mg/8.0) )   ;
			  // On the right foot
			  FC_DES[14] = - ( front_part *(mg/8.0) )   ;
			  FC_DES[17] = - ( front_part *(mg/8.0) )   ;
			  FC_DES[20] = - ( rear_part  *(mg/8.0) )   ;
			  FC_DES[23] = - ( rear_part  *(mg/8.0) )   ; 
			  //     
			  return 0 ; // FC_DES ;
    }   
     
     /**
     * @brief  provide the desired contact force distribution => 50% on the left/right Hand
     * @return  0 if ok
     */
     int FC_DES_center_hands( yarp::sig::Vector& FC_DES_hands, double mg_hands ) {
                            // On the left hand
                            FC_DES_hands[1]  = - (mg_hands/8.0 ) ; // Left Hand
                            FC_DES_hands[4]  = - (mg_hands/8.0 ) ;
                            FC_DES_hands[7]  = - (mg_hands/8.0 ) ;
                            FC_DES_hands[10] = - (mg_hands/8.0 ) ;
                            // On the right hand
                            FC_DES_hands[13] =  (mg_hands/8.0 ) ;
                            FC_DES_hands[16] =  (mg_hands/8.0 ) ;
                            FC_DES_hands[19] =  (mg_hands/8.0 ) ;
                            FC_DES_hands[22] =  (mg_hands/8.0 ) ; 
                            //     
                            return 0 ; // FC_DES_hands ;
    }   
     
     
     
     /**
     * @brief  easy way for rotating the right shoulder 
     * @param alpha rotation angle [rad], angluar step at each loop
     * @return the uptated joint vector configuration 
     */
     yarp::sig::Vector moving_right_arm( const double alpha,
                                         RobotUtils& robot,
                                         bool flag_robot, 
					 yarp::sig::Vector q_current
    ) {
                                            yarp::sig::Vector q_ref_ToMove = q_current ; 

					    yarp::sig::Vector q_ref_ToMove_right_arm(RIGHT_ARM_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_left_arm(LEFT_ARM_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_torso(TORSO_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_right_leg(RIGHT_LEG_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_left_leg(LEFT_LEG_JOINT_NUM) ;
                                        
                                            //q_ref_ToMove = q_current ; //locoman::utils::senseMotorPosition(robot, flag_robot)  ;

                                            robot.fromIdynToRobot(  q_ref_ToMove,
                                                                    q_ref_ToMove_right_arm,
                                                                    q_ref_ToMove_left_arm,
                                                                    q_ref_ToMove_torso,
                                                                    q_ref_ToMove_right_leg,
                                                                    q_ref_ToMove_left_leg  ) ; 

                                            q_ref_ToMove_right_arm[0] += alpha ;  
                       
                                            robot.fromRobotToIdyn( q_ref_ToMove_right_arm ,
                                                                q_ref_ToMove_left_arm  ,
                                                                q_ref_ToMove_torso     ,
                                                                q_ref_ToMove_right_leg ,
                                                                q_ref_ToMove_left_leg  ,
                                                                q_ref_ToMove           );   

                                            return q_ref_ToMove ; 
                                              }  

     /**
     * @brief linear function from the values (0,err_min) to (1,err_max)
     * @param err the point on which the filtering function has to be computed
     * @param err_min minimum error value (positive)
     * @param err_max maximum error value (greather than err_min)
     * @return the filtering value
     */
     double alpha_filter( double err, double err_min, double err_max ) {
                                    double tol = 0.00001 ;
                                    if(err_min <tol){ err_min = tol;}
                                    if(err_max<=(err_min+tol)){err_max = err_min + tol; }
                                    double m = 1/(err_max-err_min );
                                    double q = -m*err_min ;
                                    double alpha = m*err + q ;
                                    if(alpha<tol )    { alpha = 0.0 ; } 
                                    if(alpha>(1.0-tol)) { alpha = 1.0 ; }   
                                    return alpha ;
                                        }
                                        
     /**
     * @brief Joint_Trajectory 
     * @param robot is the robot model
     * @param flag_robot (bool) 1 if you are using the code on the real robot, 0 if you using the simulator
     * @param q_motor_current current joint configuration (motor side)
     * @param q_des desired joint configuration
     * @param steps (int) number of steps of trajectory
     * @return a yarp vector of dimension equal to robot.getNumberOfJoints() 
     */
    bool Joint_Trajectory(RobotUtils& robot, 
                          const bool flag_robot,
                          const yarp::sig::Vector& q_motor_current,
                          const yarp::sig::Vector& q_des,
                          const int steps = 300.0 , const bool verbose = 0) {           
/*              yarp::sig::Vector q_motor_current = locoman::utils::senseMotorPosition(robot, flag_robot) ; // this function uses manually imposed joint stiffness values             
                yarp::sig::Vector q_des(robot.getNumberOfJoints() ) ;     */               
                //if(1-flag_robot){steps = 100.0 ; } //faster on the simulator
                
                double steps_aux = steps ; // This is helpful to avoid casting to int the subsequent division 
                yarp::sig::Vector d_q_des = (q_des - q_motor_current); //
                for(int i =0; i <steps_aux+1; i++){
                    robot.moveNoHead( q_motor_current+(i/steps_aux)*  d_q_des) ; // robot.move(q_des) ;
                    if(verbose){std::cout << "step = " << i  << "/" << steps_aux <<std::endl ; }
                    usleep(30*1000) ;
                }
        return 0 ;
    }
    
    /**
     * @brief Tic function for tic-toc matlab style
     * @return the tic time
     */
    double Tic( ) {
        double tt_tic = std::clock() ;
        return tt_tic ;   
            }   
          /**
     * @brief 
     * @return the tic-toc time (seconds)
     */
    double Toc( double tt_tic ) {
        double tt_toc =  ( std::clock() - tt_tic ) / ((double) CLOCKS_PER_SEC )  ;
        return tt_toc ;              
            }       

     
   }   
}
#endif
