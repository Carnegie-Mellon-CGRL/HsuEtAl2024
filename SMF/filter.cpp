#include "krig.h"

// Update filter based on each simulation results
bool filt_update(filter &filt,MatD xsim,VecD ysim,VecD hsim,settings &set0)
{
// Assume filt.h & hsim contain sorted values of constraint violation aggregate
// These should not be converted to max(H,0) format-- H will be used for extracting infeasible poll candidates
// All other filter variables are stored corresponding to current filter 


// Store sizes of each 
	int npts_f = filt.y.size();
	int npts_n = ysim.size();


// Find non-dominated points in new set wrt old filter
	VecI succ_arr(npts_n,0);

	for(int i=0;i<npts_n;i++)
	{
		int curr_succ = 1;
		for(int j=0; j<npts_f;j++)
		{
			if((ysim[i]>filt.y[j])&&(std::max(hsim[i],0.0)>=std::max(filt.h[j],0.0)))
			{
				curr_succ = 0;
				break;
			}
		}
		succ_arr[i]=curr_succ ;
//		std::cout<<"Success at location "<<i+1<<" is : "<<succ_arr[i]<<" \n";
	}


// Check if there were any new non-dominated point and store their locations
// if not return

	VecI perm ; 

	int ndom_sum = 0;
	for(int i=0;i<npts_n;i++)
	{
		ndom_sum+=succ_arr[i];
		if(succ_arr[i]!=0)
			perm.push_back(i);

	}

	if(ndom_sum==0)
		return false;



// Add non-dominated points into filter


	for(int i=0;i<perm.size();i++)
	{

		VecD new_filt_h,new_filt_y;
		MatD new_filt_x ;

		int filtlen = filt.h.size();
		int loc = filtlen + 1; 

		for(int j=0;j<filt.h.size();j++)
		{


			if(std::max(filt.h[j],0.0) >= std::max(hsim[perm[i]],0.0))
			{

				loc = j;
				break;
			}

			new_filt_h.push_back(filt.h[j]);
			new_filt_y.push_back(filt.y[j]);
			new_filt_x.push_back(filt.x[j]);
		}

		new_filt_h.push_back(hsim[perm[i]]);
		new_filt_y.push_back(ysim[perm[i]]);
		new_filt_x.push_back(xsim[perm[i]]);


		for(int j=loc;j<filtlen;j++)
		{
			new_filt_h.push_back(filt.h[j]);
			new_filt_y.push_back(filt.y[j]);
			new_filt_x.push_back(filt.x[j]);

		}


		filt.h = new_filt_h ;
		filt.y = new_filt_y ;
		filt.x = new_filt_x ;

	}



// Remove dominated points from filter
// Step 1: Find non dominated points 

	VecI pts2rm; 

	for(int i=0;i<filt.h.size();i++)
	{
		int curr_ndom = 1; 

		for(int j=0;j<filt.h.size();j++)
		{
			if( (filt.y[i]>filt.y[j]) && (std::max(filt.h[i],0.0) >= std::max(filt.h[j],0.0)) && (i!=j))
			{
				curr_ndom = 0; 
				break ;
			}
		}

		if(curr_ndom == 0)
		{
			pts2rm.push_back(i);
		}

	}

	VecD new_filt_h,new_filt_y;
	MatD new_filt_x ;

// Step 2: remove non dominated points

	for(int i=0;i<filt.h.size();i++)
	{
		if(std::find(pts2rm.begin(),pts2rm.end(),i)== pts2rm.end())
		{
			// index i is not in pts2rm

			new_filt_x.push_back(filt.x[i]);
			new_filt_y.push_back(filt.y[i]);
			new_filt_h.push_back(filt.h[i]);
		}

	}

	filt.h = new_filt_h ;
	filt.y = new_filt_y ;
	filt.x = new_filt_x ;	

// The filter is non-dominated now, but needs to be sorted again

	// Extract permutation vector  

	VecI index(filt.h.size(),0);
	for(int i=0;i<index.size();i++)
	{
		index[i]=i;
	}


	sort(index.begin(),index.end(),[&]( const int& a,const int& b) {  return (filt.h[a] < filt.h[b]) ;  } ) ;


	// Sort x and y based on permutation vector from sorting h 
	for(int i=0;i<filt.h.size();i++)
	{
		new_filt_x[i] = filt.x[index[i]];
		new_filt_y[i] = filt.y[index[i]];		
	}

	print_filter_to_file("filter.txt",filt,set0.dim,filt.h.size(),set0.xlog_index);

	return true;

}
// // For debugging filter code 



// int main()
// {

// 	int dim = 6 ; 
// 	int n_init_pts = 12 ; 
// 	int n_new_pts = 5 ; 



// 	MatD filt_x;
// 	VecD filt_y,filt_h;

// 	MatD filt_x0;
// 	VecD filt_y0,filt_h0;

// 	settings set0;

// // Cook up some arbitrary initial values 


// 	std::cout<<"Initial points : \n";

// 	for(int i=0;i<n_init_pts;i++)
// 	{
// 		VecD xpts(dim,0);

// 		for(int j=0;j<dim;j++)
// 		{
// 			xpts[j] = 0.1*i + 0.01*j ;
// //			std::cout<<xpts[j]<<"  "; 
// 		}
// //		std::cout<<"\n\n";

// 		filt_x.push_back(xpts);

// 		std::cout<<"Point "<<i<<" : \n"; 

// 		for(int j= 0 ; j< dim; j++ )
// 		{
// 			std::cout<<filt_x[i][j]<<" "; 
// 		}

// 		std::cout<<"\n\n";


// 		double val = 4*sin(0.25*i*i + 0.05*i)  ;
// 		filt_y.push_back(val) ;

// 		double cons = cos(8.15/(i+0.5))*cos(8.15/(i+0.5)) +0.05*(i-8)*(i-9) -0.4;
// 		filt_h.push_back( cons ) ;

// 		std::cout<<"Value: "<<val<<"  Constraint: "<<cons<<" \n";


// 	}

// 	filter filt;

// 	bool is_succ =  filt_update(filt,filt_x,filt_y,filt_h,set0);


// 	std::cout<<"Value of filt_update: "<<is_succ<<"\n\n";
// 	std::cout<<"Filter size : "<<filt.h.size()<<" \n";

// 	for(int i=0;i<filt.h.size();i++)
// 	{
// 		std::cout<<"Value :   "<<filt.y[i]<<"  Constraint:  "<<filt.h[i]<<" \n";
// 	}


// 	// print_vector_to_txt_file("constraints.txt",filt.h, n_init_pts);

// 	// print_vector_to_txt_file("cost.txt",filt.y, n_init_pts);


// }
