#ifndef _LOCOMAN_UTILS_HPP_
#define _LOCOMAN_UTILS_HPP_

#include <cstdlib>
#include <iostream>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include <idynutils/RobotUtils.h>

#include <locoman/utils/screws.hpp>
#include <locoman/utils/kinematics.hpp>
#include <locoman/utils/kinetostatics.hpp>
// #include <locoman/utils/locoman_utils.hpp>
#include <locoman/utils/algebra.hpp>

#define RIGHT_ARM_JOINT_NUM 7
#define LEFT_ARM_JOINT_NUM 7
#define TORSO_JOINT_NUM 3
#define RIGHT_LEG_JOINT_NUM 6
#define LEFT_LEG_JOINT_NUM 6
#define HEAD_JOINT_NUM 2

namespace locoman {
        namespace utils {
            // pre declaration 
            yarp::sig::Matrix Homogeneous( const yarp::sig::Matrix& R_ab,  const yarp::sig::Vector& d_ab  );
            yarp::sig::Matrix iHomogeneous( const yarp::sig::Matrix& T_ab);
            yarp::sig::Matrix Adjoint_MT( const yarp::sig::Matrix& T_ab);
            
            
            yarp::sig::Vector sense_position_no_hands(RobotUtils& robot) {
                yarp::sig::Vector wb_input_q = robot.sensePositionRefFeedback();
                yarp::sig::Vector input_q_no_hands(robot.getNumberOfKinematicJoints());

                for(int i=0;i<input_q_no_hands.size();i++) {
                    input_q_no_hands[i]=wb_input_q[i]; //hands not considered
                }
                return input_q_no_hands;
            }
            
            yarp::sig::Vector sense_torque_no_hands(RobotUtils& robot) {
                    yarp::sig::Vector wb_input_tau = robot.senseTorque();
                    yarp::sig::Vector input_tau_no_hands(robot.getNumberOfKinematicJoints());

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
                            yarp::sig::Vector Kq_vec( robot.getNumberOfActuatedJoints()  )  ;     
                            robot.fromRobotToIdyn31( Kq_vec_right_arm ,
                                                Kq_vec_left_arm  ,
                                                Kq_vec_torso  ,
                                                Kq_vec_right_leg ,
                                                Kq_vec_left_leg  ,
                                                Kq_vec_head,
                                                Kq_vec ); 
                            yarp::sig::Matrix Kq_matrix(  robot.getNumberOfKinematicJoints(), robot.getNumberOfKinematicJoints() )  ; 
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
                            yarp::sig::Vector Kq_vec( robot.getNumberOfKinematicJoints()  )  ;     
                            robot.fromRobotToIdyn29( Kq_vec_right_arm ,
                                                Kq_vec_left_arm  ,
                                                Kq_vec_torso  ,
                                                Kq_vec_right_leg ,
                                                Kq_vec_left_leg  ,
                                                Kq_vec ); 
                            yarp::sig::Matrix Kq_matrix(  robot.getNumberOfKinematicJoints(), robot.getNumberOfKinematicJoints() )  ; 
                            Kq_matrix.diagonal(  Kq_vec ) ;
                            return Kq_matrix ;
                            }

    
     /**
     * @brief senseMotorPosition we use this method to obtain the motor position
     * @return a yarp vector of dimension equal to robot.getNumberOfJoints() 
     */
    yarp::sig::Vector senseMotorPosition(RobotUtils& robot, bool flag_robot ) {
                                        yarp::sig::Vector q_motor( robot.getNumberOfKinematicJoints()  )  ; 
                                        if(flag_robot){
                                            yarp::sig::Vector q_link = sense_position_no_hands(robot) ;
                                            return q_link ; 
                                        }
                                        else{
                                            yarp::sig::Vector q_link = sense_position_no_hands(robot) ;
                                            yarp::sig::Vector tau    = sense_torque_no_hands(robot) ;
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
                                ) {
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
                                                    yarp::sig::Vector z_imu_aw =  IMU_sense_lin_acc/norm_imu ; 
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
                                
                                
     /**
     * @brief  sigma_frict computes the metrics measuring the goodness of the contact force with respect to friction limits 
     * @param  fc is the contact force. The normal is assumed to be n = [0 0 1]^T 
     * @param  mu is the friction coefficient
     * @return the value of sigma_frict
     */
      double sigma_frict( const yarp::sig::Vector& fc ,
                          const double mu 
                                      ) {
                                        yarp::sig::Vector fc_aux = fc ;
                                        if(fc_aux(2)<0){
                                            fc_aux= -1.0*fc; 
                                        } ;
                                        yarp::sig::Vector normal(3);
                                        normal(0) = 0 ;
                                        normal(1) = 0 ;
                                        normal(2) = 1 ;
                                        //
                                        double alpha_frict = 1/(sqrt( 1+ pow(mu,2))) ;
                                        double sigma_frict = alpha_frict*yarp::math::norm(fc_aux) - yarp::math::dot(fc_aux,normal) ;
                                        return sigma_frict ;
                                                                            } 
        /**                                      
        * @brief  sigma_min computes the distance (along a certain metric) with respect to minimum force allowed 
        * @param  fc is the contact force.  
        * @param  f_min is the minimum module of the force allowed
        * @return the value of sigma_min
        */
            double sigma_min( const yarp::sig::Vector& fc,
                                const double f_min 
                                            ) {
                                        yarp::sig::Vector fc_aux = fc ;
                                        if(fc_aux(2)<0){
                                            fc_aux= -1.0*fc; 
                                        } ;
                                        yarp::sig::Vector normal(3) ;
                                        normal(0) = 0 ;
                                        normal(1) = 0 ;
                                        normal(2) = 1 ;
                                        //
                                        double sigma_min =  f_min - yarp::math::dot(fc_aux,normal) ;
                                        return sigma_min;
                                    } 
                                      
                                           /**
     * @brief  sigma_max computes the distance (along a certain metric) with respect to maximum force allowed 
     * @param  fc is the contact force.  
     * @param  f_max is the maximum module of the force allowed
     * @return the value of sigma_min
     */
      double sigma_max( const yarp::sig::Vector& fc,
                        const double f_max 
                                      ) {
                                            yarp::sig::Vector fc_aux = fc ;
                                            if(fc_aux(2)<0){
                                                fc_aux= -1.0*fc; 
                                            } ;
                                            yarp::sig::Vector normal(3) ;
                                            normal(0) = 0 ;
                                            normal(1) = 0 ;
                                            normal(2) = 1 ;
                                            //
                                            double sigma_max = - f_max +   yarp::math::dot(fc_aux,normal) ;
                                            return sigma_max ;
                                    } 
                                      
                                           /**
     * @brief  V_ij computes the distance (along a certain metric) with respect to the contact limits
     * @param  sigma is one exit of the functions sigma_frict, sigma_min, sigma_max
     * @param  toll is the admitted tolerance with respect to sigma limit value (= 0)
     * @return the value of V_ij
     */
      double V_ij( const double sigma, 
                   const double toll  = 1E-7
                 ) {
                                    double V_ij ;
                                    if(sigma<toll){
                                    V_ij = 1/(2*pow(sigma,2)) ;
                                    }
                                    else{
                                    double a = 3/(2*pow(toll,4)) ;
                                    double b = 4/(  pow(toll,3)) ;
                                    double c = 3/(  pow(toll,2)) ;
                                    V_ij = a*pow(sigma,2) + b*sigma + c ;
                                    }
                                    return V_ij ;
                                                                        } 
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
                                        yarp::sig::Vector q_motor_init(robot.getNumberOfKinematicJoints() ) ;                

                                        robot.fromRobotToIdyn29(  right_arm_configuration ,
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
     * @brief  easy way for rotating the right shoulder 
     * @param alpha rotation angle [rad], angluar step at each loop
     * @return the uptated joint vector configuration 
     */
     yarp::sig::Vector moving_right_arm( const double alpha,
                                         RobotUtils& robot,
                                         bool flag_robot
    ) {
                                            yarp::sig::Vector q_ref_ToMove(RIGHT_ARM_JOINT_NUM) ; 
                                            yarp::sig::Vector q_ref_ToMove_right_arm(RIGHT_ARM_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_left_arm(LEFT_ARM_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_torso(TORSO_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_right_leg(RIGHT_LEG_JOINT_NUM) ;
                                            yarp::sig::Vector q_ref_ToMove_left_leg(LEFT_LEG_JOINT_NUM) ;
                                        
                                            q_ref_ToMove =  senseMotorPosition(robot, flag_robot)  ;

                                            robot.fromRobotToIdyn29(  q_ref_ToMove,
                                                                    q_ref_ToMove_right_arm,
                                                                    q_ref_ToMove_left_arm,
                                                                    q_ref_ToMove_torso,
                                                                    q_ref_ToMove_right_leg,
                                                                    q_ref_ToMove_left_leg  ) ; 

                                            q_ref_ToMove_right_arm[0] += alpha ;  
                                            
                                            robot.fromRobotToIdyn29( q_ref_ToMove_right_arm ,
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
                          const int steps = 300.0 ) {
              
/*              yarp::sig::Vector q_motor_current = locoman::utils::senseMotorPosition(robot, flag_robot) ; // this function uses manually imposed joint stiffness values             
                yarp::sig::Vector q_des(robot.getNumberOfJoints() ) ;     */               
                //if(1-flag_robot){steps = 100.0 ; } //faster on the simulator
                
                double steps_aux = steps ; // This is helpful to avoid casting to int the subsequent division 
                yarp::sig::Vector d_q_des = (q_des - q_motor_current); //
                for(int i = 1; i <steps_aux+1; i++){
                    robot.move29( q_motor_current+(i/steps_aux)*  d_q_des) ; // robot.move(q_des) ;
                   
                    usleep(30*1000) ;
                }
        return 0 ;
    }                                        
     
     
     /**
     * @brief Tic fucntion for tic-toc matlab style
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
        double tt_toc =  ( std::clock() - tt_tic ) / ( CLOCKS_PER_SEC )  ;
        return tt_toc ;              
            }       
        }
}
#endif
