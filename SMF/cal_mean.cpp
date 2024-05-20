#include "krig.h"
double mean_cal(int num_pts, MatD L, VecD y0)
{
//............ CALCULATE MEAN 
        VecD ones(num_pts,1.0),mean_dummy(num_pts);
        chol_sol(num_pts,L,ones,mean_dummy);
        double mean_denom = std::accumulate(mean_dummy.begin(), mean_dummy.end(), 0.0);
        chol_sol(num_pts,L,y0,mean_dummy);
        double mean_num = std::accumulate(mean_dummy.begin(), mean_dummy.end(), 0.0);
        return mean_num/mean_denom;
}

double sigma2_cal(int num_pts, MatD L, VecD y0, double mean)
{

//............ CALCULATE SIGMA^2
	int i;

        VecD y0_minus_mean(num_pts),var_dummy(num_pts);
	for(i=0;i<num_pts;i++)
		y0_minus_mean[i] = y0[i] - mean;

        chol_sol(num_pts,L,y0_minus_mean,var_dummy);
	double total = 0;

	for(i=0;i<num_pts;i++)
		total = total + y0_minus_mean[i]*var_dummy[i];
	return total/num_pts;
}

double err_cal(model S,VecD x_trial)
{
// Calculate error at trial points
// Calculate rx 

//	std::cout<<"in err cal \n";

	VecI th_index ; 
	for(int i=0;i<S.dim;i++)
		th_index.push_back(i);

	VecD rx ; 
	expcov(S.num_pts,S.dim,S.dim,S.x0,x_trial,th_index,S.theta,rx);



// Using rx, calculate error

	
	double mean0 = mean_cal(S.num_pts,S.L,S.y0);
	double sigma2 = sigma2_cal(S.num_pts,S.L,S.y0,mean0);

	VecD v(S.num_pts) ; 
//	std::cout<<S.num_pts<<"\n";
	chol_sol(S.num_pts,S.L,rx,v);

	double norm2 = 0.0; 

	for(int i=0; i<v.size();i++)
		norm2 +=rx[i]*v[i]; 

//	std::cout<<norm2<<" is the norm \n\n";

	if(norm2 >1) 
		throw std::runtime_error("\nIll conditioning error: Error function gets imaginary values for std-dev\n");

	
	double err = -sigma2*(1.0-norm2); 
	return err; 


}


