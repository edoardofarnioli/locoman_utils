#ifndef _SCREWS_HPP_
#define _SCREWS_HPP_

#include <cstdlib>
#include <iostream>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

// #include <locoman/utils/screws.hpp>
#include <locoman/utils/kinematics.hpp>
#include <locoman/utils/kinetostatics.hpp>
#include <locoman/utils/locoman_utils.hpp>
#include <locoman/utils/algebra.hpp>

namespace locoman {
        namespace utils {
                    /**
                    * @brief  getRot extracts the rotation matrix from a homogenous matrix
                    * @param  T_ab is a 4x4 yarp matrix describing a homogenous transformation
                    * @return R_ab is a 3x3 yarp matrix
                    */
                    yarp::sig::Matrix getRot( const yarp::sig::Matrix& T_ab) {
                                                yarp::sig::Matrix R_ab( 3 , 3 ) ;  
                                                R_ab = T_ab.submatrix(0,2,0,2) ;
                                                return R_ab ;
                           }
            //------------------------------------------------------------
    
                    /**
                    * @brief  getTrasl extracts the translation vector from a homogenous matrix
                    * @param  T_ab is a 4x4 yarp matrix describing a homogenous transformation
                    * @return d_ab  is a 3x1 yarp vector 
                    */
                    yarp::sig::Vector getTrasl( const yarp::sig::Matrix& T_ab) {
                            yarp::sig::Matrix d_ab_mat( 3, 1 ) ;  
                            d_ab_mat = T_ab.submatrix(0,2,3, 3) ;   
                            yarp::sig::Vector d_ab( 3  ) ;  
                            d_ab[0] = d_ab_mat[0][0] ;
                            d_ab[1] = d_ab_mat[1][0] ;
                            d_ab[2] = d_ab_mat[2][0] ;
                            return d_ab ;                             
                    }
                        //------------------------------------------------------------
    
                    /**
                    * @brief  homogeneous computes the homogenous transformation composed by:
                    * @param  R_ab is a 3x3 yarp matrix describing the rotational part
                    * @param  d_ab is a 3x1 yarp vector describing the translational part   
                    * @return T_ab is a 4x4 yarp matrix 
                    */
                    yarp::sig::Matrix Homogeneous( const yarp::sig::Matrix& R_ab,  const yarp::sig::Vector& d_ab  ) {
                        
                        yarp::sig::Matrix T_ab( 4 , 4 ) ;
                        yarp::sig::Matrix d_ab_mat( 3 , 1 ) ;
                        T_ab.setSubmatrix(R_ab, 0, 0 ) ;
                        d_ab_mat[0][0] = d_ab[0] ;
                        d_ab_mat[1][0] = d_ab[1] ;
                        d_ab_mat[2][0] = d_ab[2] ;
                        
                        T_ab.setSubmatrix(d_ab_mat, 0, 3 ) ;
                        T_ab[3][3] = 1 ; 
                        return T_ab ; 
                                        }
                   
                       //-----------------------------------------------------------------
     /**
     * @brief  iHomogeneous computes the inverse of an homogenous transformation
     * @param  T_ab is a 4x4 yarp matrix describing an homogenous transformation
     * @return T_ba  is a 4x4 yarp matrix, the inverse of T_ab
     */
    yarp::sig::Matrix iHomogeneous( const yarp::sig::Matrix& T_ab) {
                    
                        yarp::sig::Matrix T_ba( 4 , 4 ) ;  
                        yarp::sig::Matrix R_ab( 3 , 3 ) ; 
                        yarp::sig::Vector d_ab( 3  ) ;  
                        yarp::sig::Matrix R_ba( 3 , 3 ) ; 
                        yarp::sig::Vector d_ba( 3  ) ;    
                        
                        R_ab = locoman::utils::getRot(T_ab) ;
                        d_ab = locoman::utils::getTrasl(T_ab)  ;
                        
                        R_ba = R_ab.transposed() ; 
                        d_ba = -1.0*R_ba*d_ab ; // 
                        
                        T_ba = locoman::utils::Homogeneous(R_ba, d_ba) ;
                        //
                        return T_ba ; 
                        }
                             //-----------------------------------------------------------------
     /**
     * @brief  Adjoint computes the adjoint matrix := Ad_g from a homogenous transformation
     * @param  T_ab is a 4x4 yarp matrix describing an homogenous transformation
     * @return Ad_T_ab is a 6x6 yarp matrix, able to map twists in twists
     */
    yarp::sig::Matrix Adjoint( const yarp::sig::Matrix& T_ab) {    yarp::sig::Matrix Adj( 6 , 6 ) ;
                                yarp::sig::Matrix R_ab( 3 , 3 ) ; 
                                yarp::sig::Vector d_ab( 3  ) ;  
                                yarp::sig::Matrix d_ab_skew( 3 , 3 )   ;
                                
                                R_ab = locoman::utils::getRot(T_ab) ;
                                d_ab = locoman::utils::getTrasl(T_ab)  ;
                                d_ab_skew = yarp::math::crossProductMatrix( d_ab ) ; 
                                
                                Adj.setSubmatrix(R_ab, 0 , 0 ) ;
                                Adj.setSubmatrix(R_ab, 3 , 3 ) ;
                                Adj.setSubmatrix(d_ab_skew*R_ab, 0 , 3 ) ;
                                
                                return Adj ;

    }
    
     //-----------------------------------------------------------------
     /**
     * @brief  Adjoint_MT computes the inverse transpose of the adjoint  matrix := (Ad_(g^{-1}))^{T} from a homogenous transformation
     * @param  T_ab is a 4x4 yarp matrix describing an homogenous transformation
     * @return ( Ad_(T_ab)^{-1} )^{T}  is a 6x6 yarp matrix, able to map wrenches in wrenches
     */
    yarp::sig::Matrix Adjoint_MT( const yarp::sig::Matrix& T_ab) {
                                yarp::sig::Matrix Adj_MT( 6 , 6 ) ;
                                yarp::sig::Matrix R_ab( 3 , 3 ) ; 
                                yarp::sig::Vector d_ab( 3  ) ;  
                                yarp::sig::Matrix d_ab_skew( 3 , 3 )   ;
                                
                                R_ab = locoman::utils::getRot(T_ab) ;
                                d_ab = locoman::utils::getTrasl(T_ab)  ;
                                d_ab_skew = yarp::math::crossProductMatrix( d_ab ) ; 
                                
                                Adj_MT.setSubmatrix(R_ab, 0 , 0 ) ;
                                Adj_MT.setSubmatrix(R_ab, 3 , 3 ) ;
                                Adj_MT.setSubmatrix(d_ab_skew*R_ab, 3 , 0 ) ;
                                
                                return Adj_MT ;
                                }
        
         //-----------------------------------------------------------------
     /**
     * @brief  xi_hat returns the homogeneous form of a tiwst
     * @param  xi is a 6 dimentional yarp vector describing a twits = [v^T, w^T]^T
     * @return 4x4 yarp matrix describing the homogenous form of a the twist xi
     */
    yarp::sig::Matrix xi_hat( const yarp::sig::Vector& xi) {
                                yarp::sig::Vector vel = xi.subVector(0,2)   ;
                                yarp::sig::Vector omega  = xi.subVector(3,5)  ;
                                yarp::sig::Matrix vel_matr(3,1)  ;
                                vel_matr[0][0] = vel[0] ;
                                vel_matr[1][0] = vel[1] ;
                                vel_matr[2][0] = vel[2] ;
                                yarp::sig::Matrix xi_hat_matr(4,4) ;
                                yarp::sig::Matrix omega_hat = yarp::math::crossProductMatrix(omega) ;
                                xi_hat_matr.setSubmatrix( omega_hat, 0, 0 )   ;
                                xi_hat_matr.setSubmatrix( vel_matr, 0, 3 )   ;
                                return xi_hat_matr ;
                        }
         //-----------------------------------------------------------------
     /**
     * @brief  exp_omega_theta returns the rotation matrix provided by the Rodrigues formula
     * @param  omega is the rotation axis
     * @param  theta is the rotation amount
     * @return 4x4 yarp matrix describing the homogenous form of a the twist xi
     */
    yarp::sig::Matrix exp_omega_theta( const yarp::sig::Vector& omega, const double theta) {
                                                yarp::sig::Matrix Eye_3(3,3) ;
                                                Eye_3.eye() ;
                                                yarp::sig::Matrix R(3,3) ;
                                                yarp::sig::Matrix omega_hat = yarp::math::crossProductMatrix(omega) ;
                                                R = Eye_3 + (omega_hat/yarp::math::norm(omega) )*sin(yarp::math::norm(omega)*theta) + 
                                                            omega_hat*omega_hat/
                                                            (yarp::math::norm(omega)*yarp::math::norm(omega))*
                                                            (1 - cos(yarp::math::norm(omega)*theta)) ; 
                                                return R ;
                                        }
    
         //-----------------------------------------------------------------
     /**
     * @brief  exp_xi_theta returns the homogenous transformation associated to the twist
     * @param  xi is the twist = [v^T, w^T]^T 
     * @param  theta is the transformation amount
     * @return 4x4 yarp matrix describing the homogenous transformation
     */
    yarp::sig::Matrix twistexp( const yarp::sig::Vector& xi, const double theta) {
                                            yarp::sig::Matrix Eye_3(3,3) ;
                                            Eye_3.eye() ;
                                            yarp::sig::Vector vel = xi.subVector(0,2)   ;
                                            yarp::sig::Vector omega  = xi.subVector(3,5)  ;
                                            yarp::sig::Matrix omega_matr(3,1)  ;
                                            omega_matr[0][0] = omega[0] ;
                                            omega_matr[1][0] = omega[1] ;
                                            omega_matr[2][0] = omega[2] ;
                                            yarp::sig::Matrix R(3,3) ;
                                            yarp::sig::Vector t(3) ;
                                            if (yarp::math::norm(omega)< 1E-8)
                                            {
                                                R = Eye_3 ;
                                                t = vel*theta ;
                                            }  
                                            else {
                                            R = locoman::utils::exp_omega_theta(omega, theta) ;
                                            // std::cout << "qui_3 " << std::endl ;
                                            t = (Eye_3 - 1.0*R )*(yarp::math::crossProductMatrix(omega)*vel ) + omega_matr*omega_matr.transposed()*vel *theta  ;
                                            }
                                            yarp::sig::Matrix T(4,4) ;
                                            T = locoman::utils::Homogeneous(R, t) ;
                                            return T ;
                                                }
     //----------------------------------------------------------------------------
     /**
     * @brief  ad_lie transforms a tiwst given as a yarp matrix into the Lie adjoint matrix
     * @param  Xi is a 6x1 yarp matrix describing a twist
     * @return 6xc yarp matrix 
     */
    yarp::sig::Matrix ad_lie( const yarp::sig::Matrix& Xi) {
                                    if(!(Xi.rows() == 6)) std::cout << "ERROR DIMENSIONS OF Xi" << std::endl; 
                                    if(!(Xi.cols() == 1)) std::cout << "ERROR DIMENSIONS OF Xi" << std::endl; 
                                    yarp::sig::Vector v(3) ;
                                    yarp::sig::Vector omega(3) ;
                                    v[0] = Xi[0][0];
                                    v[1] = Xi[1][0];
                                    v[2] = Xi[2][0];
                                    omega[0] = Xi[3][0];
                                    omega[1] = Xi[4][0];
                                    omega[2] = Xi[5][0];
                                    yarp::sig::Matrix v_skew = yarp::math::crossProductMatrix(v) ;
                                    yarp::sig::Matrix omega_skew = yarp::math::crossProductMatrix(omega) ;
                                    yarp::sig::Matrix Zero_3_3(3,3) ;
                                    Zero_3_3.zero();
                                    yarp::sig::Matrix AD_LIE(6,6) ;
                                    AD_LIE.setSubmatrix(omega_skew, 0, 0 ) ;
                                    AD_LIE.setSubmatrix(v_skew,     0, 3 ) ;
                                    AD_LIE.setSubmatrix(Zero_3_3 ,  3, 0 ) ;
                                    AD_LIE.setSubmatrix(omega_skew ,3, 3 ) ;
                                    return AD_LIE ;
    }
    
    yarp::sig::Matrix ad_lie( const yarp::sig::Vector& Xi) {
                                if(!(Xi.length() == 6)) std::cout << "ERROR DIMENSIONS OF Xi" << std::endl; 
                                yarp::sig::Vector v(3) ;
                                yarp::sig::Vector omega(3) ;
                                v[0] = Xi[0] ;
                                v[1] = Xi[1] ;
                                v[2] = Xi[2] ;
                                omega[0] = Xi[3] ;
                                omega[1] = Xi[4] ;
                                omega[2] = Xi[5] ;
                                yarp::sig::Matrix v_skew = yarp::math::crossProductMatrix(v) ;
                                yarp::sig::Matrix omega_skew = yarp::math::crossProductMatrix(omega) ;
                                yarp::sig::Matrix Zero_3_3(3,3) ;
                                Zero_3_3.zero();
                                yarp::sig::Matrix AD_LIE(6,6) ;
                                AD_LIE.setSubmatrix(omega_skew, 0, 0 ) ;
                                AD_LIE.setSubmatrix(v_skew,     0, 3 ) ;
                                AD_LIE.setSubmatrix(Zero_3_3 ,  3, 0 ) ;
                                AD_LIE.setSubmatrix(omega_skew ,3, 3 ) ;
                                return AD_LIE ;
    }
         //----------------------------------------------------------------------------
     /**
     * @brief  D_Jacob_spa_i compute the derivative of the spatial Jacobian with respect to the i-th q
     * @param  J_s is a 6xc yarp matrix describing a spatial Jacobian
     * @param  i -th joint () with repect to the derivative is computed 
     * @return 6x6 yarp matrix describing the Lie adjoint matrix
     */
    yarp::sig::Matrix D_Jacob_spa_i( const yarp::sig::Matrix& J_s, const int i ) {
                                                    
                                        int r_J = J_s.rows() ;
                                        int c_J = J_s.cols() ;
                                        yarp::sig::Matrix D_J_i(r_J,c_J) ;
                                        yarp::sig::Vector zero_6(6, 0.0) ;
                                        D_J_i.zero();
                                        if(!(r_J==6)){std::cout << "ERROR DIMENSIONS OF Xi" << std::endl ;  }  
                                        if(i>=c_J ){ return D_J_i;    }
                                        else{  
                                            for ( int k = 0  ; k<i ; k++ )  
                                        {
                                            D_J_i.setCol(k, zero_6) ;
                                        };
                                        for(int k = i  ; k<c_J ; k++  )
                                        {    
                                        yarp::sig::Matrix temp = locoman::utils::ad_lie(J_s.submatrix( 0, 5, i-1, i-1 ))*J_s.submatrix(0,5, k , k )  ;
                                        yarp::sig::Vector temp_2(6) ;
                                        temp_2[0] =  temp[0][0] ;
                                        temp_2[1] =  temp[1][0] ;
                                        temp_2[2] =  temp[2][0] ;
                                        temp_2[3] =  temp[3][0] ;
                                        temp_2[4] =  temp[4][0] ;
                                        temp_2[5] =  temp[5][0] ;
                                        D_J_i.setCol(k,  temp_2) ;
                                        };
                                        return D_J_i;
                                        }
                                            }
    
         //----------------------------------------------------------------------------
     /**
     * @brief  AdjToPose transforms an Adjoint matrix into the original homogenous matrix
     * @param  Adj is a 6x6 yarp matrix describing an Adjoint transformation
     * @return 4x4 yarp matrix describing a homogenous transformation
     */
    yarp::sig::Matrix AdjToPose( const yarp::sig::Matrix& Adj) {
                                        yarp::sig::Matrix T(4,4) ;
                                        T.setSubmatrix( Adj.submatrix(0,2,0,2), 0,0 ); // = Adj.submatrix(0,2,0,2) ;
                                        yarp::sig::Matrix R = locoman::utils::getRot(T)  ;
                                        //std::cout << "R = "  << R.toString()  << std::endl ;
                                        yarp::sig::Matrix Temp = Adj.submatrix(0,2,3,5) ;
                                        //std::cout << "Temp = "  << Temp.toString()  << std::endl ;
                                        yarp::sig::Matrix Temp_2 =  Temp*R.transposed(); 
                                        //std::cout << "Temp_2 = "  <<  Temp_2.toString()  << std::endl ;
                                        yarp::sig::Vector d = SkewToVect(Temp_2) ;
                                        //std::cout << "d = "  << d.toString()  << std::endl ;
                                        yarp::sig::Matrix d_matr(3,1) ;
                                        d_matr[0][0] = d[0] ;
                                        d_matr[1][0] = d[1] ;    
                                        d_matr[2][0] = d[2] ;
                                        //    std::cout << "d_matr = "  << d_matr.toString()  << std::endl ;
                                        T.setSubmatrix(d_matr , 0,3)  ; 
                                        //    std::cout << "T = "  << T.toString()  << std::endl ;
                                        T[3][3] = 1 ;
                                        return T ;
                                        }
            
        }
        
}

#endif
