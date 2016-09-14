#ifndef _LOCO_DECLARATIONS_H_
#define _LOCO_DECLARATIONS_H_

#include <yarp/sig/all.h>
#include <idynutils/RobotUtils.h>

namespace locoman {
        namespace utils {
        yarp::sig::Matrix RoundMatrix( const yarp::sig::Matrix& M, const int k);
        yarp::sig::Matrix Pinv_trunc_SVD( const yarp::sig::Matrix& A, const double k);
        yarp::sig::Matrix Pinv_Regularized( const yarp::sig::Matrix& A, const double k );
        yarp::sig::Matrix Pinv_Marq( const yarp::sig::Matrix& A, const double k );
        yarp::sig::Vector x_Pinv_Iter( const yarp::sig::Matrix& A ,const yarp::sig::Vector& b ,double n);
        yarp::sig::Matrix orth_SVD( const yarp::sig::Matrix& A ,const double k);
        yarp::sig::Matrix null_SVD( const yarp::sig::Matrix& A, const double k);
        yarp::sig::Matrix filter_SVD( const yarp::sig::Matrix& A, const double k);

        yarp::sig::Matrix getRot( const yarp::sig::Matrix& T_ab) ;
        yarp::sig::Vector getTrasl( const yarp::sig::Matrix& T_ab) ;
        yarp::sig::Matrix Homogeneous( const yarp::sig::Matrix& R_ab,  const yarp::sig::Vector& d_ab  ) ;
        yarp::sig::Matrix iHomogeneous( const yarp::sig::Matrix& T_ab) ;
        yarp::sig::Matrix Adjoint( const yarp::sig::Matrix& T_ab) ;
        yarp::sig::Matrix Adjoint_MT( const yarp::sig::Matrix& T_ab) ;
        yarp::sig::Matrix xi_hat( const yarp::sig::Vector& xi) ;
        yarp::sig::Matrix exp_omega_theta( const yarp::sig::Vector& omega, const double theta) ;
        yarp::sig::Matrix twistexp( const yarp::sig::Vector& xi, const double theta) ;
        yarp::sig::Matrix ad_lie( const yarp::sig::Matrix& Xi) ;
        yarp::sig::Matrix ad_lie( const yarp::sig::Vector& Xi);
        yarp::sig::Matrix D_Jacob_spa_i( const yarp::sig::Matrix& J_s, const int i );
        yarp::sig::Matrix AdjToPose( const yarp::sig::Matrix& Adj);

        yarp::sig::Matrix Q_ci( const yarp::sig::Matrix& J_spa_i, 
                            const yarp::sig::Matrix& T_a_ci , 
                            const yarp::sig::Vector& f_ci) ;
        yarp::sig::Matrix Q_ci_wrench( const yarp::sig::Matrix& J_spa_i, 
                                    const yarp::sig::Matrix& T_a_ci , 
                                    const yarp::sig::Vector& w_ci);
        yarp::sig::Matrix FLMM_ext( const yarp::sig::Matrix& J_c ,
                            const yarp::sig::Matrix& S_c ,
                            const yarp::sig::Matrix& Q_j,
                            const yarp::sig::Matrix& Q_s,
                            const yarp::sig::Matrix& U_j,
                            const yarp::sig::Matrix& U_s,
                            const yarp::sig::Matrix& K_c,
                            const yarp::sig::Matrix& K_q
                              ) ;
        yarp::sig::Matrix Rf_ext( const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_j,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_j,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c,
                                const yarp::sig::Matrix& K_q
                                );
        yarp::sig::Matrix FLMM_redu( const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c
                                );
        yarp::sig::Matrix Rf_redu( const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c
                                );
        yarp::sig::Matrix Ru_redu( const yarp::sig::Matrix& J_c ,
                                const yarp::sig::Matrix& S_c ,
                                const yarp::sig::Matrix& Q_s,
                                const yarp::sig::Matrix& U_s,
                                const yarp::sig::Matrix& K_c
                                );

        yarp::sig::Matrix getKq_Walkman(RobotUtils& robot ) ;
        yarp::sig::Matrix getKq(RobotUtils& robot );
        yarp::sig::Vector senseMotorPosition(RobotUtils& robot, bool flag_robot ) ;
        yarp::sig::Matrix fConToSens( const int sens_index,
                                    const int c1_index,
                                    const int c2_index,
                                    const int c3_index,
                                    const int c4_index,
                                    iDynUtils& model 
                                    // bool flag_robot
                                    ) ;
        yarp::sig::Matrix AW_world_posture( iDynUtils model , 
                                            RobotUtils& robot  ) ;
        double sigma_frict( const yarp::sig::Vector& fc ,
                            const double mu 
                                        ) ;
            double sigma_min( const yarp::sig::Vector& fc,
                                    const double f_min 
                                                ) ;
        double sigma_max( const yarp::sig::Vector& fc,
                            const double f_max 
                                        );
        double V_ij( const double sigma, 
                    const double toll
                    ) ;
            yarp::sig::Vector q_init(  RobotUtils& robot,
                                    yarp::sig::Vector right_arm_configuration,
                                    yarp::sig::Vector left_arm_configuration, 
                                    yarp::sig::Vector torso_configuration, 
                                    yarp::sig::Vector right_leg_configuration,
                                    yarp::sig::Vector left_leg_configuration 
                                    ) ;
        int FC_DES_right( yarp::sig::Vector& FC_DES, double mg  );
        int FC_DES_left( yarp::sig::Vector& FC_DES, double mg   ) ;
        int FC_DES_center( yarp::sig::Vector& FC_DES, double mg );
        yarp::sig::Vector moving_right_arm( const double alpha,
                                            RobotUtils& robot,
                                            bool flag_robot
        ) ;
        double alpha_filter( double err, double err_min, double err_max );
        bool Joint_Trajectory(RobotUtils& robot, const bool flag_robot, const yarp::sig::Vector& q_motor_current, const yarp::sig::Vector& q_des, const int steps ) ;

        yarp::sig::Vector SkewToVect( const yarp::sig::Matrix& Skew) ;
        yarp::sig::Vector Rot2Quat( const yarp::sig::Matrix& Rot) ;
        yarp::sig::Vector Orient_Error( const yarp::sig::Matrix& Rot_des,
                                        const yarp::sig::Matrix& Rot_cur) ;
        yarp::sig::Vector Inv_quaternion( const yarp::sig::Vector& quat ) ;
        yarp::sig::Vector quat_Product( const yarp::sig::Vector& quat_1 , 
                                        const yarp::sig::Vector& quat_2) ;
        yarp::sig::Matrix Rot_x( const double phi_x );
        yarp::sig::Matrix Rot_y( const double theta_y );
        yarp::sig::Matrix Rot_z( const double psi_z);
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
                                            );
    }
}
#endif //_LOCO_DECLARATIONS_H_