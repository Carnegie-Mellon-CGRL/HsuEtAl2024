#include "krig.h"

void set_th_index(VecI &th_index, int dim)
{

//	MODIFY IN ACCORDANCE TO TH_INDEX THAT YOU WANT	
	for(int i=0;i<dim;i++)
		th_index[i] = i;
}

void set_xlog(VecI &xlog_index)
{
	int dim = xlog_index.size();

	for(int i =0;i<dim;i++)
		xlog_index[i]=0;	
}

void set_th_bounds(VecD &thlow,VecD &thupp,int dim_th)
{
//	MODIFY IN ACCORDANCE TO TH_BOUNDS THAT YOU WANT
//  ALSO KEEP TRACK OF MAPPING BETWEEN TH_INDEX TO ACTUAL THETA	
	for(int i=0;i<dim_th;i++)
	{
		thlow[i] = 1e-4;
		thupp[i] = 200;		
	}

}

void find_xmin(MatD x,VecD J, VecD H, VecD &xmin, double &Jmin, double &Hmin)	// Find set of parameters for current (possibly infeasible) incumbent solution
{
	int num_pts = J.size();
	int min_pos = 0;
	Jmin = 1E10;
	for(int i=0;i<num_pts;i++)
	{
		if(Jmin > J[i])
		{
			Jmin = J[i];
			min_pos = i;
		}
	}

	xmin = x[min_pos];
	Hmin = H[min_pos];

}


bool find_feas_xmin(MatD x,VecD J, VecD H, VecD &xmin, double &Jmin)	// Find set of parameters for current feasible incumbent solution
{
	// Assume H contains raw augmented constraint values- constraint satisfaction => H[i]<=0

	int num_pts = J.size();
	int min_pos = 0;
	Jmin = 1E10;
	for(int i=0;i<num_pts;i++)
	{
		if((Jmin > J[i])&&(H[i]<=0.0) )
		{
			Jmin = J[i];
			min_pos = i;
		}
	}

	if(Jmin==1E10)
	{
		std::cout<<"No feasible solution yet. \n";
		return false ;
	}
	else
		xmin = x[min_pos];

	std::cout<<"Curr feas minimum value :  "<<Jmin<<"\n\n";

	return true;

}


bool newpoint(int num_pts, int dim, MatD x0, VecD x_cand)
{
        double x_dist = 0.0;
        double x_dist_tol = 1e-7;
        for(int i =0;i<num_pts;i++)
        {
                x_dist = 0;
                for(int j=0;j<dim;j++)
                        x_dist += pow(x0[i][j] - x_cand[j],2);
                x_dist = pow(x_dist,0.5);
                if(x_dist < x_dist_tol)
                        return false;
        }
        return true;

}

void terminate_opt(int terminate)
{
        VecD term(1);
        term[0] = (double)terminate;

        print_vector_to_txt_file("stopfile",term,1);

}

bool succ_search(settings &set0,VecI xlog_index)
{
	int dim = set0.dim;
	filter filt;
	MatD x,xnew;
	VecD H,J,Hnew,Jnew;
	int num_pts,num_pts_y,num_pts_new;
	bool is_success;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	import_vector_from_txt_file("Hhist.dat",H,num_pts_y);


	import_matrix_from_txt_file("xnew.dat",xnew,num_pts_new,dim,xlog_index);
	import_vector_from_txt_file("Jnew.dat",Jnew,num_pts_new);
	import_vector_from_txt_file("Hnew.dat",Hnew,num_pts_new);


	int n_pts_f;
	import_filter_from_file("filter.txt",filt,n_pts_f,set0.dim,set0.xlog_index);


	// -------TODO--- Success effects for filter
	bool update_counter = filt_update(filt,xnew,Jnew,Hnew,set0);

	VecD xmin,xmin_new;
	double Jmin,Jmin_new,Hmin,Hmin_new;
	find_xmin(x,J,H,xmin,Jmin,Hmin);
	find_xmin(xnew,Jnew,Hnew,xmin_new,Jmin_new,Hmin_new);

	int indx_jmin_new;
	for(int i = 0; i < Jnew.size();i++)
	{
		if(Jnew[i]==Jmin_new)
			indx_jmin_new = i ;
	} 


	if(update_counter)
	{
		is_success = true;
		set0.succ_search += 1 ;

		int success_id = set0.npts_tot+indx_jmin_new-Jnew.size()+1 ;

		set0.succ_pts_id.push_back(success_id);
		set0.succ_delta.push_back(Jmin-Jmin_new);


	}

	else
		is_success = false;


	for (int i = 0; i < num_pts_new; i++)
	{
		x.push_back(xnew[i]);
		J.push_back(Jnew[i]);
		H.push_back(Hnew[i]);

	}


	num_pts += num_pts_new;
	num_pts_y += num_pts_new;

	print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	print_vector_to_txt_file("Jhist.dat",J,num_pts_y);
	print_vector_to_txt_file("Hhist.dat",H,num_pts_y);
	return is_success;


}

bool succ_poll(settings &set0,char polltype)
{

	VecI xlog_index = set0.xlog_index;
	int dim =set0.dim;
	filter filt;
    MatD x,xnew;
    VecD J,Jnew,H,Hnew;
    int num_pts,num_pts_y,num_pts_new;
	bool is_success;
	
	// Read history and candidate data

	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	import_vector_from_txt_file("Hhist.dat",H,num_pts_y);

	import_matrix_from_txt_file("xnew.dat",xnew,num_pts_new,dim,xlog_index);
	import_vector_from_txt_file("Jnew.dat",Jnew,num_pts_new);
	import_vector_from_txt_file("Hnew.dat",Hnew,num_pts_new);

	
	int n_pts_f;
	import_filter_from_file("filter.txt",filt,n_pts_f,dim,set0.xlog_index);


	// Check if new minimum is achieved

	// -------TODO Success condition for filter

	VecD xmin,xmin_new;
	double Jmin,Jmin_new,Hmin,Hmin_new;
	find_xmin(x,J,H,xmin,Jmin,Hmin);
	find_xmin(xnew,Jnew,Hnew,xmin_new,Jmin_new,Hmin_new);



	bool update_counter = filt_update(filt,xnew,Jnew,Hnew,set0);
	// For book keeping, trace-out location of current minimum in new points

	int indx_jmin_new;
	for(int i = 0; i < Jnew.size();i++)
	{
		if(Jnew[i]==Jmin_new)
			indx_jmin_new = i ;
	} 

	if(update_counter)
	{
		is_success = true;
		int success_id = set0.npts_tot+indx_jmin_new-Jnew.size()+1 ;

		if(polltype == 'I')
			set0.succ_pollinf += 1 ;
		else if(polltype == 'F')
			set0.succ_pollfeas += 1 ;

		set0.succ_pts_id.push_back(success_id);
		set0.succ_delta.push_back(Jmin-Jmin_new);

	}
	else
		is_success = false;

	for (int i = 0; i < num_pts_new; i++)
	{
		x.push_back(xnew[i]);
		J.push_back(Jnew[i]);
		H.push_back(Hnew[i]);
	}

	// Update counters

	num_pts += num_pts_new;
	num_pts_y += num_pts_new;

	print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	print_vector_to_txt_file("Jhist.dat",J,num_pts_y);
	print_vector_to_txt_file("Hhist.dat",H,num_pts_y);

	return is_success;
}



// --- FUNCTION FOR UPDATING HYPER-PARAMETERS AND MINIMIZING SURROGATE

bool minimize_surr(VecD &xkrig, settings &set0, int num_pts,MatD x,VecD J, VecD H, model &S)
{

	int dim = set0.dim;
	int dim_th = set0.dim_th;
	double likeli_max;
	VecD xmin(dim);

	double Jmin,Hmin;
	find_xmin(x,J,H,xmin,Jmin,Hmin);

	VecD theta = set0.theta;

	
// ============== GENERATE SURROGATE FOR THE POINTS ===========
	while(!param2model(S,num_pts,dim,dim_th,x,J,set0.th_index,theta))
	{
		std::cout<<" PARAMETER TO MODEL - FAILED: CHOLESKY. \n \n";
		for(int p = 0 ; p < dim_th ; p++)
	    {
		
			if(theta[p] > set0.thupp[p])
			{
				terminate_opt(-1);	
	       		return false;
			}
			else
				theta[p] = theta[p]*pow(2,(1.0/(dim-1)));
	    }
    }

	int nbrestarts = 0;
	double incpopsize = 4.0;
	double *theta_d;
	double tol_th = 1e-4 + 0.00001*num_pts;
	VecD theta_stddev(dim_th,0.1);
	std::cout<<" HYPER PARAMETER OPTIMIZATION : \n ";
	theta_d =cm_optimize(dim_th, S, loglikelihood, set0.thlow, set0.thupp, nbrestarts, incpopsize,theta,theta_stddev,tol_th);  
	std::cout<<" THETA OBTAINED : \n";

	for(int i=0;i<dim_th;i++)
	{
		theta[i] = theta_d[i];
		std::cout<<theta[i]<<"  ";
	}
	std::cout<<"\n \n";

    if(!param2model(S,num_pts,dim,dim_th,x,J,set0.th_index,theta))
    {
            std::cout<<" PARAMETER TO MODEL - FAILED: CHOLESKY. \n \n";
            return false;
    }


	double mse = msecheck(S);
	std::cout<<"\nMSE of new model=  "<<mse<<std::endl;

	set0.theta = theta ; 

//============== CALCULATE NEW CANDIDATE MINIMUM ==============

	// Multi-start search -- Do multiple searches by large scale perturbations to the current best
	
    nbrestarts = 0;
    incpopsize = 4.0;
    double *xkrig_d;
	double tolx = 1e-7;
	std::cout<<"SURROGATE MINIMIZATION SEARCH STEP \n";


	double currbest = 1E10 ;
	// Perturb aroudn current best point 
	VecD xmin_pert = xmin ; 
	double stdfact = 0.25 ;

	VecD stdvec(dim);

	for(int k=0; k< dim ;k++)
	{
		stdvec[k] = stdfact*(set0.xupp[k] - set0.xlow[k]);	
	}

	for(int k =0 ; k< set0.npts_surr_search ; k++)
	{

		VecD x_stddev = stdvec ;
		if(k>0)
		{
			
			for(int i=0;i<dim;i++)
			{
				xmin_pert[i] = xmin[i] +randd(-1*stdvec[i],1*stdvec[i]) ;  
			}
		}

		VecD tmpx(dim);
		VecD tmpout(1); 
   		xkrig_d =cm_optimize(dim, S, pred_fn, set0.xlow, set0.xupp, nbrestarts, incpopsize, xmin_pert, x_stddev, tolx);

 	  	for(int i=0;i<dim;i++)
 	   		tmpx[i] = xkrig_d[i];

 	   	MatD xmattmp ; 
 	   	xmattmp.push_back(tmpx) ;

 	   	predictor(S,1, xmattmp,tmpout);

 	   	if(tmpout[0]- currbest< -0.001*abs(currbest))
 	   	{
 	   		xkrig = tmpx;
 	   		currbest = tmpout[0]; 
 	   		std::cout<<"Multipoint CMAES : Improved candidate point found \n";
 	   	}
   	}

    push_to_grid(xkrig,set0.xlow,set0.xupp,set0.spc);

    return true;

}


// --- FUNCTION FOR CALCULATING POINTS WITH MAXIMUM ERRORS IN SURROGATE

void surr_err_minimize(settings &set0,model &S,MatD &xnew,int n_eval)
{

	int dim = set0.dim;
	int dim_th = set0.dim_th;
	double likeli_max;

//============== PARTITION DOMAIN ============================

//---- Number of divisions = 2^k, where k is the smallest int for which 2^k >= n_eval
// i.e. k = ceil(log(n_eval)/log2). A small factor is subtracted for roundoff in case of n_eval being a power of 2 

	int part_dim = ( ceil( log(n_eval)/log(2) - 1E-4) );
	int n_part = pow(2,part_dim);


	VecI perm,part_vec;
	for(int i=0;i<dim;i++)
	{
		perm.push_back(i);
	}

	std::random_shuffle (perm.begin(),perm.end()); 

// Use the first k dimensions to partition 

	for(int i=0;i<part_dim;i++)
		part_vec.push_back(perm[i]);


//============== CALCULATE CANDIDATE MINIMA ==============

	// Multi-start search -- Do multiple searches by large scale perturbations to the current best
	// Settings for CMAES 
        int nbrestarts = 0;
        double incpopsize = 4.0;
        double *xkrig_d;
	double tolx = 1e-7;
	std::cout<<"SURROGATE ERROR INFILL SEARCH STEP \n";

	double currbest = 1E10 ;
	VecD xstart(dim,0);
	VecD part_xupp;
	VecD part_xlow; 
	double stdfact = 0.25 ;

	VecD stdvec(dim);

	for(int k=0; k<dim ;k++)
	{
		stdvec[k] = stdfact*(set0.xupp[k] - set0.xlow[k]);	

	}

	// Create temporary point array and counter 

	MatD tmp_x0 = S.x0; 
	int tmp_npts = S.num_pts; 

	if(xnew.size()!=0)
	{
		for(int k=0 ; k<xnew.size();k++)
		{
			tmp_x0.push_back(xnew[k]);
			tmp_npts++;
		}

	}


	for(int k=0 ; k<n_eval ; k++)
	{

		VecD x_stddev = stdvec ;
		part_xupp=set0.xupp;
		part_xlow=set0.xlow; 
		int k2 = k ;

		for(int i=0; i<part_vec.size();i++)
		{
			if(  k2%2 == 0 ) 
			{
				part_xupp[part_vec[i]] = 0.5*(part_xlow[part_vec[i]] + part_xupp[part_vec[i]]);
//				std::cout<<part_xupp[part_vec[i]]<<"    upp \n" ; 
			}
			else
			{
				part_xlow[part_vec[i]] = 0.5*(part_xlow[part_vec[i]] + part_xupp[part_vec[i]]);
//				std::cout<<part_xlow[part_vec[i]]<<"    low \n" ;
			}
			k2 = k2/2 ;
		}


		for(int i=0;i<dim;i++)
		{
			xstart[i] = 0.5*(part_xlow[i] + part_xupp[i]) +randd(-1*stdvec[i],1*stdvec[i]) ;  
		}

		VecD tmpx(dim);
		VecD tmpout(1); 
//		std::cout<<"Surr before opt \n\n";
//		std::cout<<"neval = "<<n_eval<<"\n\n";

   		xkrig_d =cm_optimize(dim, S, err_fn, part_xlow, part_xupp, nbrestarts, incpopsize, xstart, x_stddev, tolx);

//		std::cout<<"Surr after opt \n\n";


 	  	for(int i=0;i<dim;i++)
 	   		tmpx[i] = xkrig_d[i];

	    push_to_grid(tmpx,set0.xlow,set0.xupp,set0.spc);

	    if(newpoint(tmp_npts,dim,tmp_x0,tmpx))
	    {
	    	xnew.push_back(tmpx);
		tmp_x0.push_back(tmpx);
		tmp_npts++;  

	    }
	    else
	    {
	    	MatD poll_points; 
	    	// Poll around until we find a non-repeated point nearby 
	    	while(poll_points.size()==0)
	    	{
			std::cout<<"Surrogate error yields repeating points. Poll now: \n";
			bool *mads_flag = new bool;

	    		int pollstg = 1+randn(2);
	    		poll_mads(tmpx,dim,set0.spc,part_xlow,part_xupp,pollstg,poll_points, mads_flag);
			delete mads_flag;

			if(poll_points.size()!=0)
			{
                                // Empty poll array if new point is not obtained
				if(!newpoint(tmp_npts,dim,tmp_x0,poll_points[0]))
					poll_points.clear() ; 

			}
	    	}
	    	xnew.push_back(poll_points[0]);
	    }

   	}
}


void search(settings &set0)
{

//			Search process breakdown : 
//				1 point from minimizing surrogate
//				(n) points from minimizing error from n different partitions	

	int dim = set0.dim; 
	int dim_th = set0.dim_th;
	model S;

	MatD x,xnew;
	VecD J,H,Jnew,Hnew,theta;

	VecD xkrig(dim);

	int  num_pts,num_pts_y,dim_th_chk;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,set0.xlog_index);
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
	import_vector_from_txt_file("Hhist.dat",H,num_pts_y);

	set0.nstage_search += 1 ; 

	theta = set0.theta;

	if(num_pts == num_pts_y)
		std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
	else
	{
		throw std::runtime_error("\nInput error : number of points in xhist dont match Jhist \n");
	}

//=============== CALCULATE CURRENT SURROGATE MINIMIZER ======

	int succ_minimize = 0;
	if(minimize_surr(xkrig,set0,num_pts,x,J,H,S))
	{
		if(newpoint(num_pts,dim,x,xkrig) == true)
		{
			xnew.push_back(xkrig);
			succ_minimize=1;
		}
		else
			std::cout<<"Minimizing the surrogate produced a pre-existing point \n\n";
	}


//=============== CALCULATE ERROR FROM SURROGATE ============
// ------- TODO : Add error based search part 

	int pts_err = dim+1-succ_minimize ;
	surr_err_minimize(set0,S,xnew,pts_err);


// xnew should now have n+1 candidates 
	if(xnew.size()>dim+1)
		throw std::runtime_error("Error function has some bugs. Please fix me.");	

	std::cout<<"Search points generated: \n\n";
	for(int j=0;j<xnew.size();j++)
	{
		std::cout<<"Point number"<<j+1<<"\n";
		for(int i=0; i<dim;i++)
		{
			std::cout<<xnew[j][i]<<"  "; 
		}
		std::cout<<"\n \n";

		set0.npts_search +=1 ; 
		set0.pt_hist += 'S' ;
	}

    print_theta_hist(S);

    print_matrix_to_txt_file("xnew.dat",xnew,xnew.size(),dim,set0.xlog_index);



}

void poll(settings &set0,char polltype, bool *poll_flag)
{
	int dim = set0.dim;

	VecI xlog_index = set0.xlog_index;
	filter filt;
	MatD x;
	VecD J,H;
	int  num_pts_y,num_pts;
	import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
	import_vector_from_txt_file("Hhist.dat",H,num_pts_y);	
	import_vector_from_txt_file("Jhist.dat",J,num_pts_y);

//	int n_pts_f;
//	import_filter_from_file("filter.txt",filt,n_pts_f,dim,set0.xlog_index);


	if(num_pts == num_pts_y)
		std::cout<<" \nNumber of function evaluations till now :   "<<num_pts<<"\n \n";
	else
	{
		std::cout<<" \nInput error : number of points in xhist dont match Jhist \n";
		return;
	}

    MatD poll_points,poll_to_eval;
    VecD xmin(dim);
    double Jmin,Hmin;
    int eval_size = 0;


// --- POLL BASED ON CHARACTER PROVIDED 
    std::cout<<"Before poll \n";

    int poll_stage = std::max(std::max(set0.nstage_pollfeas,set0.nstage_pollinf),1);

    if(polltype=='I')
	{
	    find_xmin(x,J,H,xmin,Jmin,Hmin);
		bool *mads_flag = new bool;
    	poll_mads(xmin,dim,set0.spc,set0.xlow,set0.xupp,poll_stage,poll_points, mads_flag);	
		if (*mads_flag == 1) {
			*poll_flag = 1;
			delete mads_flag;
			return;
		}
	    for(int i=0;i<dim+1;i++)
	    {
	    	if((newpoint(num_pts,dim,x,poll_points[i]))&&(bndchk(dim,poll_points[i],set0.xlow,set0.xupp)))
	    	{
	    		poll_to_eval.push_back(poll_points[i]);
	    		eval_size++;
	    	}
	    }
	}

	else
	{
		find_feas_xmin(x,J,H,xmin,Jmin);
		bool *mads_flag = new bool;
    	poll_mads(xmin,dim,set0.spc,set0.xlow,set0.xupp,poll_stage,poll_points, mads_flag);
		if (*mads_flag == 1) {
			*poll_flag = 1;
			delete mads_flag;
			return;
		}
	    for(int i=0;i<dim+1;i++)
	    {
	    	if((newpoint(num_pts,dim,x,poll_points[i]))&&(bndchk(dim,poll_points[i],set0.xlow,set0.xupp)))
	    	{
	    		poll_to_eval.push_back(poll_points[i]);
	    		eval_size++;
	    	}
	    }
	}



    std::cout<<"Poll eval size : "<<eval_size<<"\n";

    if(polltype=='I')
    {
		set0.npts_pollinf += eval_size ; 
		set0.nstage_pollinf += 1 ; 
		for(int i=0;i<eval_size;i++)
			set0.pt_hist += 'I';    	
    }
    else
    {
		set0.npts_pollfeas += eval_size ; 
		set0.nstage_pollfeas += 1 ; 
		for(int i=0;i<eval_size;i++)
			set0.pt_hist += 'F';    	

    }

    	print_matrix_to_txt_file("xnew.dat",poll_to_eval,eval_size,dim,xlog_index);

}

void refine(int dim,VecI &spc)		// FOR UNSUCCESSFUL SEARCH AND UNSUCCESSFUL POLL
{
	for(int i=0;i<dim;i++)
		spc[i] = (spc[i]-1)*4 + 1;
}



void coarsen(int dim,VecI &spc)		// FOR EITHER A SUCCESSFUL SEARCH OR A SUCCESSFUL POLL
{
	double coarsen_fact = 1;

	std::cout<<"SPC[0] BEFORE : "<<spc[0]<<std::endl;
	for(int i=0;i<dim;i++)
		spc[i] = (spc[i]-1)/coarsen_fact + 1;
	std::cout<<"SPC[0] AFTER COARSEN : "<<spc[0]<<std::endl;
}

bool filt_feas_chk(settings set0)
{
	int n_pts_f;
	filter filt;
	import_filter_from_file("filter.txt",filt,n_pts_f,set0.dim,set0.xlog_index);

	for(int i=0;i<filt.h.size();i++)
	{
		// Check for any feasible points
		if(filt.h[i]<=0)
			return true;
	}

	return false;

}


int main()
{


//	SEED RANDOM GENERATOR FIRST
	std::srand(std::time(NULL));

   	double Jmin=1e10,Jnew,Jkrig;
	int i,j,num_pts=0;
	MatD x;
	VecD J,H;


// 	READ SETTINGS FROM INPUT FILE
	settings set0 ;
	read_phase("smf.inp",set0);

//	SETUP FILTER VARIABLES
	filter filt;


	if(set0.dim ==0 )
	{
		throw std::runtime_error("The problem dimension was not initialized \n");
	}


//	DEFAULT FUNCTIONAL VALUES FOR NUM-LHS POINTS AND DIMENSION IN HYPER-PARAMETER SPACE

	if(set0.ndiv==0)
	{
		set0.ndiv = 3*set0.dim + 1;
	}


	if(set0.dim_th==0)
	{
		set0.dim_th = set0.dim;
	}


	int dim = set0.dim;
	int dim_th = set0.dim_th;

	set0.npts_tot =  set0.npts_lhs + set0.npts_search + set0.npts_pollfeas + set0.npts_pollinf ; 


//	THETA INDEX VECTOR
	VecI th_index(dim);

	set_th_index(th_index,dim);
	set0.th_index = th_index;

//	PARAMETER INDEX - TOGGLE FOR 0 - NORMAL; 1 -LOG SCALE

	VecI xlog_index(dim,0);
	set_xlog(xlog_index);	
	set0.xlog_index = xlog_index ;



//  SET BOUNDS FOR THETA

	if((set0.thlow.size() == 0)||(set0.thupp.size() == 0) )
	{	
		VecD thlow(dim_th);
		VecD thupp(dim_th);

		set_th_bounds(thlow,thupp,dim_th);

		set0.thupp = thupp ;
		set0.thlow = thlow ;

	}


	if(set0.theta.size() == 0)
	{

		//  INITIALIZE THETA
		VecD theta(dim_th);
		for(i=0;i<dim_th;i++)
			theta[i]= sqrt(set0.thupp[i]*set0.thlow[i]);

		set0.theta = theta ;


	}


	VecD xmin(dim),xkrig(dim);




	FILE *fp;
	fp = fopen("BestJHist.dat","a");



	int coarsen_flag=0;

//=============== DEFINE BOX CONSTRAINTS ======================


	std::cout<<"Past the read phase \n";

	bool dim_bds = ((set0.xupp.size()!=dim)||(set0.xlow.size()!=dim));

    if(dim_bds)
    {
        throw std::runtime_error("Dimensions of upper bounds/lower bounds are incorrect \n");
    }



//=============== SMF PARAMETERS ==============================

	double conv_tol = 0.00001;
	double x_del_sum = 0;
	VecD x_delta(dim);


//============== TERMINATION CONDITION ======================

    for(i=0;i<dim;i++)
    {
            x_delta[i] = (set0.xupp[i] - set0.xlow[i])/(set0.spc[i]-1);
            x_del_sum += x_delta[i];
    }


	if(x_del_sum <= dim*conv_tol)
	{
		std::cout<<"Insert post-processing stuff here \n";

		terminate_opt(0);
		return 0;
	}
	std::cout<<"X_DEL_SUM = "<<x_del_sum<<std::endl;


//=============== READ CURRENT STATE IN SMF LOOP =============
	bool *pflag = new bool;
	if (set0.phase==0)  //  FIRST EXECUTION OF COMMAND
	{

		std::cout<<"LHS initiated \n";
		num_pts = set0.ndiv;
		lhs(dim, num_pts, set0, 1, x, xlog_index)	;	
		
		print_matrix_to_txt_file("xnew.dat",x,num_pts,dim,xlog_index);
		
 
		set0.phase = 1;

		write_phase("smf.inp",set0);

		return 0;
	}

	else if (set0.phase==1)	// LHS COMPLETED  -- DO SEARCH
	{
		std::cout<<" \nLHS completed. Updating filter \n";

		int num_pts,num_pts_y;
		import_matrix_from_txt_file("xnew.dat",x,num_pts,dim,xlog_index);
		import_vector_from_txt_file("Jnew.dat",J,num_pts_y);
		import_vector_from_txt_file("Hnew.dat",H,num_pts_y);

		if(num_pts!=num_pts_y)
			std::cout<<"ERROR: Diff in xnew and Jnew = "<<num_pts-num_pts_y<<std::endl;

		print_matrix_to_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
		print_vector_to_txt_file("Jhist.dat",J,num_pts);
		print_vector_to_txt_file("Hhist.dat",H,num_pts);

		filt_update(filt,x,J,H,set0);

		search(set0);
		set0.phase = 2 ;

	}

	else if (set0.phase==2)	// SEARCH COMPLETED -- CHECK IF POLL IS NEEDED
	{

		int n_pts_f;
		bool search_flag = succ_search(set0,xlog_index);

		std::cout<<"Checking if search worked \n";

		if(search_flag == true)
		{
			std::cout<<"Search worked \n";
			coarsen_flag=1;
			search(set0);
			set0.phase = 2 ; 

		}
		else
		{
			std::cout<<"Search failed. Polling \n";
			bool feas_chk = filt_feas_chk(set0);

			if(feas_chk)		// Filter has feasible entries in it 
			{
				std::cout<<"Doing feasible poll \n";
				poll(set0,'F', pflag);
				set0.phase = 3;
			}
			else
			{
				std::cout<<"No current feasible points in filter. Doing infeasible poll \n";
				poll(set0,'I', pflag);
				set0.phase = 4;
			}
		}

	}
	else if (set0.phase==3)	// POLL FEASIBLE COMPLETED -- CHECK IF POLL WORKED
	{

		int n_pts_f;
		std::cout<<"Checking if feasible poll worked \n";
		bool poll_flag = succ_poll(set0,'F');
		if (*pflag == false) {
			poll_flag = false;
		}
		if(poll_flag == true)
		{
			std::cout<<"Poll worked. Do search \n";
			coarsen_flag=1;
			search(set0);
			set0.phase = 2; 

		}
		else
		{
			// Check if current best J is feasible 

			filter filt;
			int n_pts_f;
			import_filter_from_file("filter.txt",filt,n_pts_f,set0.dim,set0.xlog_index);
			VecD xmin;
			double Jmin,Hmin;
			find_xmin(filt.x,filt.y,filt.h,xmin,Jmin,Hmin);

			if(Hmin<=0)
			{

				std::cout<<"Feasible Poll failed. Refine and search\n";
				refine(dim,set0.spc);
	
				search(set0);
				set0.phase = 2; 

			}
			else
			{
				std::cout<<"Feasible Poll failed. Do infeasible poll\n";
				poll(set0,'I', pflag);
				set0.phase = 4;

			}
		}

	}

	else if (set0.phase==4)	// POLL INFEASIBLE COMPLETED -- CHECK IF POLL WORKED
	{

		int n_pts_f;

		std::cout<<"Checking if infeasible poll worked \n";
		bool poll_flag = succ_poll(set0,'I');
		if (*pflag == false) {
			poll_flag = false;
		}
		if(poll_flag == true)
		{
			std::cout<<"Infeasible Poll worked. Do search \n";
			coarsen_flag=1;
			search(set0);
			set0.phase = 2; 

		}
		else
		{
			std::cout<<"Infeasible Poll failed. Refine and search\n";
			refine(dim,set0.spc);
			search(set0);
			set0.phase = 2; 
		}
	
	}
	

	else if (set0.phase == 5) // RESTART CASE
	{

		std::cout<<"Restarted SMF! \n";

		search(set0);
		set0.phase = 2;

	}


	else
	{
		std::cout<<"Error. The value of phase is : "<<set0.phase<<"\n";
		terminate_opt(set0.phase);
	}
	delete pflag;

	if(coarsen_flag==1)
		coarsen(dim,set0.spc);


	// ------- WRITE PHASE FILE 

	write_phase("smf.inp",set0);



	if(set0.phase>=2)
	{
		VecD xmin_final;
		double Jmin_final,Hmin_final;
		int num_pts,num_pts_y;
		import_matrix_from_txt_file("xhist.dat",x,num_pts,dim,xlog_index);
		import_vector_from_txt_file("Jhist.dat",J,num_pts_y);
		import_vector_from_txt_file("Hhist.dat",H,num_pts_y);

		find_xmin(x,J,H,xmin_final,Jmin_final,Hmin_final);
		write_hist(fp,xmin_final,Jmin_final,xlog_index);

	    std::cout<<"Jmin till now : "<<Jmin_final<<" with constraint value: "<<Hmin_final<<"  at xmin : \n ";
       	for(int p =0 ; p < dim ; p++)
       		std::cout<<xmin_final[p]<<"  ";
       	std::cout<<std::endl;
	}


	fclose(fp);
}
