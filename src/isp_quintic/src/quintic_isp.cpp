
/**
\file   quintic_spline.cpp
\brief  quintic_spline parameterization
 *
\author  Mahmoud Ali
\date    21/6/2019
*/


#include<vector>
#include<iostream>
#include "ros/ros.h"
#include "std_msgs/Float64MultiArray.h"

//#include<normal_toppra_traj_instant_1.h>
#include<normal_toppra_traj_instant_3.h>
#include <Eigen/Dense>
#include <math.h>

#include <python2.7/Python.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace Eigen;
using namespace std;
using Eigen::MatrixXd;

const double sm=180, vm= 130, am=250, jm=1000; // pos, vel, acc, jrk max limits


// finds in which segment time instant tg belons to
int find_seg_number (double tg, std::vector<double> T_vec){
    // segment changes from 0 to n-2 which mean n-1 seg
    int seg = 0;
    if(tg<= T_vec[0]) //less than tstart
        return seg;
    else if(tg>= T_vec[T_vec.size() -2 ] ) //if tg is greater than the starting time of the last segment
        return T_vec.size() - 2;
    else {
        for (int i=0; i< T_vec.size()-2 ; i++) { // between tstart and tend
            if(tg>= T_vec[i] && tg< T_vec[i+1])
                seg = i;
        }
        return seg;
    }
}


// to sample trajectory, find state(pos, vel, acc, jrk) at a given time tg, T_vec is time vector of waypoints, cof is coef matrix of quintic spline
void sample (const std::vector<double>& T_vec, const double& tg, const std::vector<std::vector<double>> &cof, std::vector<double> & state){//, bool pos_plt=true, bool vel_plt=true, bool acc_plt=true, bool jrk_plt=true ){
    std::vector<double> a=cof[0], b=cof[1], c=cof[2], d=cof[3], e=cof[4], f=cof[5];
    double pos=0, vel=0, acc=0, jrk=0;
    int n = T_vec.size(), seg=0; // n:number of points
    double t=0; // to represent tg-tseg
    //check for vector sizes
    //    std::cout<<"\nsize of a= "<< a.size() << " size of T= "<< T_vec.size();
    if( a.size() != T_vec.size()-1 )
        std::cout<<"error: coef and time have different size\n";


    // check for tg inside T_vec or not, if yes find corrresponding segment
    if(tg < T_vec[0]){
        std::cout<<"Warn: request time instant is before start time of trajectory\n";
        pos= a[0];
        vel= b[0];
        acc= 2*c[0];
        jrk= 6*d[0];
    }
    else if(tg > T_vec[n-1] ){
        std::cout<<"Warn: request time instant is after the end time of trajectory\n ";
        t = tg - T_vec[n-1];
        pos=   a[n-2] +   b[n-2]*t +         c[n-2]*pow( t, 2) +        d[n-2]*pow( t, 3) +    e[n-2]*pow( t, 4) +    f[n-2]* pow( t, 5)  ;
        vel=              b[n-2]        +  2*c[n-2]*t          +      3*d[n-2]*pow( t, 2) +  4*e[n-2]*pow( t, 3) +  5*f[n-2]* pow( t, 4)  ;
        acc=                               2*c[n-2]                 + 6*d[n-2]*t          + 12*e[n-2]*pow( t, 2) + 20*f[n-2]* pow( t, 3)  ;
        jrk=                                                          6*d[n-2]                 + 24*e[n-2]*t          + 60*f[n-2]* pow( t, 2)  ;
    }
    else{
        seg= find_seg_number(tg, T_vec);
        t = tg - T_vec[seg];
        //        std::cout<<"      tg= " <<tg<< "   seg_number= "<<  seg<<"   t= "<<t<<std::endl;
        pos=   a[seg] +   b[seg]*t +    c[seg]*pow( t, 2) +   d[seg]*pow( t, 3) +    e[seg]*pow( t, 4) +    f[seg]* pow( t, 5) ;
        vel=              b[seg]        +  2*c[seg]*t          + 3*d[seg]*pow( t, 2) +  4*e[seg]*pow( t, 3) +  5*f[seg]* pow( t, 4) ;
        acc=                               2*c[seg]                 + 6*d[seg]*t          + 12*e[seg]*pow( t, 2) + 20*f[seg]* pow( t, 3) ;
        jrk=                                                          6*d[seg]                 + 24*e[seg]*t          + 60*f[seg]* pow( t, 2) ;
    }
    // assign state to state vector
    state[0] = pos;
    state[1] = vel;
    state[2] = acc;
    state[3] = jrk;

}



// generate C1 matrix which is depends on the time of each segment, idx is segment indx, h vectore contains segments times
void generate_coef_matrix_c1(const int idx, const VectorXd &h, MatrixXd &C1){
    MatrixXd C(5,6);
    C<< 1,  h(idx),  pow(h(idx), 2),    pow(h(idx), 3),    pow(h(idx), 4),     pow(h(idx), 5),
            0,       1,        2*h(idx),  3*pow(h(idx), 2),  4*pow(h(idx), 3),   5*pow(h(idx), 4),
            0,       0,               1,          3*h(idx),  6*pow(h(idx), 2),  10*pow(h(idx), 3),
            0,       0,               0,                 1,          4*h(idx),  10*pow(h(idx), 2),
            0,       0,               0,                 0,                 1,           5*h(idx);
    C1 = C;
}















int main(int argc, char **argv)
{
    ros::init(argc, argv, "quintic_isp");
    ros::NodeHandle nh;

    ROS_INFO("initializing path .... ");
    //============ read trajectory ==========
    trajectory_msgs::JointTrajectory  traj;
    traj = generate_traj();
    int n_jts = traj.joint_names.size();
    int n_pts = traj.points.size();

    // get waypoints from trajectory
    std::vector< std::vector<double> > P_jt_wpt;
    P_jt_wpt.resize( n_jts);
    for(int jt=0; jt<n_jts; jt++){
        for(int pt=0; pt<n_pts; pt++){
            P_jt_wpt[jt].push_back( traj.points[pt].positions[jt] );
        }
    }


    //pos vector
    VectorXd pos(n_pts);
    for (int pt=0; pt< n_pts; pt++) {
        pos(pt) = P_jt_wpt[4][pt];
    }
    //    for (int pt=0; pt< n_pts; pt++) {
    //        pos(pt+n_pts) = P_jt_wpt[5][pt];
    //    }
    // put equations in form of AX=B
    int n = pos.size(); //number of waypoints
    int n_eqs = 6*(n-1); // number of equations/variables
    double vel0=0, veln=0, acc0=0, accn=0; // initial/final  conditions(vel, acc)

    // store stop point at which trajectory reverse direction
    std::vector<int> stp_pts_idx;
    for (int i=1; i<n-1; i++) {
        if( (  (pos(i+1) - pos(i))> 0 && (pos(i) - pos(i-1))> 0  ) || ( (pos(i+1) - pos(i))< 0 && (pos(i) - pos(i-1))< 0 ) )
            continue;
        stp_pts_idx.push_back( i );
    }
    //    ROS_INFO_STREAM("stp_pts_idx"); //check
    //    for (int i=0; i< stp_pts_idx.size(); i++)
    //        ROS_INFO_STREAM(" "<< stp_pts_idx[i] );



    //from pos vector: generate initial time for way points. t: time vector of waypoinys times, h is vector for segments times
    VectorXd t(n), h(n-1); //t0 ==> tn-1  ...   h0 ==> hn-2 (#n of t and #n-1 of h)
    for (int i = 0; i < n-1; ++i) {
        h(i) = abs( (abs(pos(i+1)) - abs(pos(i))) ) /vm;
    }
    //  ROS_INFO_STREAM("h vector: \n"<<h);
    t(0)=0; //starting time of trajectory
    for (int i = 0; i < n-1; ++i) {
        t(i+1) = h(i) + t(i);
    }
    //  ROS_INFO_STREAM("t vector: \n"<<t);


    // define coef matrix and B matrix: AX=B
    VectorXd B(n_eqs), X(n_eqs);
    MatrixXd A =   MatrixXd::Zero(n_eqs, n_eqs);
    MatrixXd C1 =  MatrixXd::Zero(5, 6);
    MatrixXd C2 =  MatrixXd::Zero(5, 6);

    // C1 and C2 store coef of 5(n-2) equations
    for (int i=0; i<5; i++)
        C2(i,i)=-1.0; // C2(0,0)  to   C2(4,4)

    for (int i=0; i<n-1; i++) {
        generate_coef_matrix_c1(i, h, C1);
        ROS_INFO_STREAM(" C1("<<i <<") initialization:\n "<< C1);
    }


    //>>>>>>>>>>>step_1:  3 boundary conditions for first point (pos0, vel0, acc0), 3 rows
    A(0,0) = 1;
    A(1,1) = 1;
    A(2,2) = 1;  // 3 rows are filled so far, so next r=3


    //>>>>>>>>>>> step_2:  pos condition for all waypoints except first and last, (n-2) rows
    int col_idx= 6;
    for (int r = 3; r < 3+(n-2); ++r) {  //start filling from 4th row, #rows= n-2
        A(r,col_idx) = 1;
        col_idx+=6;
    }  // n-2 rows are filled so far, so total rows are n-2+3= n+1


    int count =0; // count no of iteration
    std::vector<int> seg_idx;  //store segments indices that violate jerk
    bool update_all_seg_time= false; // select: true=> updates times of all segments, false: update time for segments that voilate jerk
    //================================================= loop: update time till the jerk become within the limit ================================
    // bif while loop
    bool lmtd_jrk = false; //will be set to true when all jrk values are in the limits
    while (!lmtd_jrk){
        if(count>0){ //no update in the first loop
            if(update_all_seg_time) // update all segment time
                for (int idx = 0; idx < n-1; ++idx)
                    h(idx)=h(idx)+0.01;

            else{  // only update the time of segment that violate jerk limit
                for (int idx = 0; idx < seg_idx.size(); ++idx){
                    ROS_INFO_STREAM("seg_idx["<<idx<< "]:  "<<seg_idx[idx]);
                    if (seg_idx[idx] == n-1)
                        seg_idx[idx] = seg_idx[idx] -1;
                    h(seg_idx[idx])=h(seg_idx[idx])+0.01;
                }
            }
        }
        seg_idx.clear();
        for (int i = 0; i < n-1; ++i) {
            t(i+1) = h(i) + t(i);
        }

        // checks for velocity and acc limits

        //>>>>>>>>>>>step_3: boundary conditions for last point (posn, veln, accn), 3 rows (depends on time so it is inside the loop)
        // ROS_INFO("3 boundary conditions for last poit");
        for (int r = n+1 ; r < 3+(n+1); ++r) {// for 3 row
            for (int c = n_eqs-6 ; c < n_eqs; ++c) {  // last 6 columns
                generate_coef_matrix_c1( n-2, h, C1);
                A(r, c) = C1(r-(n+1), c-(n_eqs-6));     // last h ==> hn-2
            }
        }//3+(n-2)+3 = n+4 rows are filled so far,  so next row to fill is (n+5)th row [its is index= n+4]
        //>>>>>>>>>>>step_4: 5(n-2) coninuity conditions for each waypoints except first and last, 5(n-2) rows (depend on time so it is inside the loop)
        //   ROS_INFO(" 5(n-2) continuity conditions for midle poits");
        int row_idx = n+4; col_idx=0;  //start filling from (n+2)th row [its is index= n+1],
        for (int blk=0; blk< n-2; blk++) { // n-2 blks, each blk is 5 equations with 6 coef
            //       ROS_INFO_STREAM("row_idx: "<< row_idx <<"   col_idx:"<< col_idx);
            for (int r = 0 ; r < 5 ; ++r) {  // 5(n-2) rows
                for (int c = 0 ; c < 6; ++c) {  // last 6 columns
                    //               ROS_INFO_STREAM("r: "<< r <<"   c:"<< c);
                    generate_coef_matrix_c1( blk, h, C1);
                    A(row_idx+r, col_idx+c)   = C1(r,c);     // last h ==> hn-2
                    A(row_idx+r, col_idx+c+6) = C2(r,c);
                }
            }
            col_idx += 6;
            row_idx += 5;
        } //all row s are filled


        //====================================check for stop points:==================================


        //>>>>>>>>>>>step_5: fill B vector B
        B(0) = pos(0);
        B(1) = vel0;
        B(2) = acc0; //  3 rows so far
        for (int i=1; i< n-1; i++) // all pos from pos1 to pos(n-1)= (n-2) rows
            B(i+2)= pos(i); // so far 3+n-2 = n+1

        for (int i=n+4; i< n_eqs; i++) // last 5(n-2) ros
            B(i)= 0; // so far 3+n-2 = n+1

        B(n+1) = pos(n-1);
        B(n+2) = veln;
        B(n+3) = accn;
        ROS_INFO_STREAM("B vector: \n"<<B);


        //find the solution: using matirx inverse
        //   MatrixXd A_inv = A.inverse();
        //   X = A.inverse()*B;  // some times gives nan or inf
        //   ROS_INFO_STREAM("X vector: \n"<<X);

        HouseholderQR<MatrixXd> qr(A);
        X = qr.solve(B); // computes A^-1 * b
        ROS_INFO_STREAM("X vector: \n"<<X );



        // extract spline coefs from the vector  ========================================================
        std::vector<double> a, b, c, d, e, f;
        for (int i = 0; i < n_eqs; i+=6) {
            a.push_back( X(i) );
            b.push_back( X(i+1) );
            c.push_back( X(i+2) );
            d.push_back( X(i+3) );
            e.push_back( X(i+4) );
            f.push_back( X(i+5) );
        }


        // sample trajectory using sample func: create variables,  sampling T_vec, Coef_mtrx
        std::vector<double> T_vec;
        T_vec.resize(n);
        for (int i=0; i< n; i++) {
            T_vec[i] = t(i);
        }
        std::vector<double> POS, VEL, ACC, JRK, T;
        std::vector<std::vector<double>> cof_mtrx;
        cof_mtrx.push_back(a);
        cof_mtrx.push_back(b);
        cof_mtrx.push_back(c);
        cof_mtrx.push_back(d);
        cof_mtrx.push_back(e);
        cof_mtrx.push_back(f);
        std::vector<double> state;
        state.resize(4);
        // states at each intermdiate times
        double frq =125, tg=0;
        while (tg< T_vec[n-1]) {
            sample(T_vec, tg, cof_mtrx, state);
            //    std::cout<< "state:"<< "p= "<< state[0]<<"   v= "<< state[1]<<"   a= "<< state[2]<<"  j= "<< state[3] ;
            T.push_back(tg);
            POS.push_back( state[0]);
            VEL.push_back( state[1]);
            ACC.push_back( state[2]);
            JRK.push_back( state[3]);
            tg += 1/frq;
        }


        // states at waypoints and check limit conditions
        lmtd_jrk = true;
        std::vector<std::vector<double>> wpt_states;
        wpt_states.resize(4);
        for (int pt = 0; pt < n; ++pt) {
            tg= T_vec[pt];
            sample(T_vec, tg, cof_mtrx, state);
            wpt_states[0].push_back( state[0] );
            wpt_states[1].push_back( state[1] );
            wpt_states[2].push_back( state[2] );
            wpt_states[3].push_back( state[3] );
            std::cout<< "jrk["<<pt<<"]= "<< state[3] ;
            if(  state[3] > 1000 || state[3] < -1000 ){
                seg_idx.push_back( pt);
                lmtd_jrk = false;
            }
        }
        // ROS_INFO( "waypoints times:    ");
        // for (int i=0; i< n; i++ )
        //    ROS_INFO_STREAM( "    "<< T_vec[i] );
        // ROS_INFO( "segment peroids:    ");
        // for (int i=0; i< n-1; i++ )
        //    ROS_INFO_STREAM( "    "<< h[i] );

        // plot section =================================
        count++;
        int plot_per_loop = 3; //to plot 9 iteration to fit on sbplts number
        int sbplt_rows= 3, sbplt_cols= 3;
        bool pos_plt=false, vel_plt=true, acc_plt=true, jrk_plt = false;// select states to plot
        pos_plt=false, vel_plt=false, acc_plt=false, jrk_plt = true;
        std::string  sbplt_name;

        //   here choose few iterations, plot them  in subplot
        if(  (count)%plot_per_loop == 0 ){
            plt::figure(1);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop   );
            plt::plot( T, POS);
            plt::plot( T_vec, wpt_states[0],"r*");
            sbplt_name= "pos_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);
            plt::figure(2);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop   );
            plt::plot(  T, VEL);
            plt::plot( T_vec, wpt_states[1],"r*");
            sbplt_name= "vel_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);
            plt::figure(3);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop   );
            plt::plot(  T, ACC);
            plt::plot( T_vec, wpt_states[2],"r*");
            sbplt_name= "acc_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);
            plt::figure(4);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop   );
            plt::plot(  T, JRK);
            plt::plot( T_vec, wpt_states[3],"r*");
            sbplt_name= "jrk_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);
            //    plt::legend(); // Enable legend.

        }
        if(lmtd_jrk){ // to plot last iteration
            plt::figure(1);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop  +1);
            plt::plot( T, POS);
            plt::plot( T_vec, wpt_states[0],"r*");
            sbplt_name= "pos_last_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);

            plt::figure(2);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop    +1);
            plt::plot(  T, VEL);
            plt::plot( T_vec, wpt_states[1],"r*");
            sbplt_name= "vel_last_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);

            plt::figure(3);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop    +1);
            plt::plot(  T, ACC);
            plt::plot( T_vec, wpt_states[2],"r*");
            sbplt_name= "acc_last_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);

            plt::figure(4);
            plt::subplot(  sbplt_rows, sbplt_cols,count/plot_per_loop    +1);
            plt::plot(  T, JRK);
            plt::plot( T_vec, wpt_states[3],"r*");
            sbplt_name= "jrk_last_iter_" + std::to_string(count);
            plt::title(  sbplt_name);
            plt::grid(true);

            //    plt::legend(); // Enable legend.
        }
        ROS_INFO_STREAM("#iteration= "<<count <<"is_lmtd_jrk: "<< lmtd_jrk);
    }

    ROS_INFO_STREAM("#loops= "<<count);
    plt::show();


}

