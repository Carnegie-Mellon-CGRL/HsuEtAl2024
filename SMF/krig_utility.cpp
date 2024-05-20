#include <string>
#include <sstream>
#include "krig.h"

// ====================================================================================================
//........RANDOM NUMBER GENERATOR
// ====================================================================================================


int randn(int n) 
{ 
	// Time based seeding
//	std::srand(std::time(0));
	return std::rand()%n; 
}

double randd (double lower_bound, double upper_bound)
{
    // Time based seeding
//    std::srand(std::time(0));

    double f = (double)std::rand() / RAND_MAX;
    return lower_bound + f * (upper_bound - lower_bound);
}



// ====================================================================================================
//........FILE HANDLING FUNCTIONS
// ====================================================================================================

// void read_phase(const char* filename_p,int &phase)
// {
//     std::ifstream file_p;
//     std::string line;
//     file_p.open(filename_p);
//     if(file_p.is_open())
//     {
//         std::getline(file_p,line);
//         if (!line.compare("SEARCH"))
//             phase = 2;
//         else if (!line.compare("POLL"))
//             phase = 3;
//         else if (!line.compare("LHS"))
//             phase = 1;
//         else if (!line.compare("RESTART"))
//             phase = 4;
//         else
//         {
//             std::cout<<"Incorrect phase info \n";
//             phase=-1;
//         }
//     }
//     else
//     {        
//         std::cout<<"No phase file present \n ";
//         std::cout<<"Begin with LHS phase \n";
//         phase = 0;
// //        std::cout<<"Begin LHS phase then \n ";
//     }
//     file_p.close();
// }

// void write_phase(const char* filename_p, const char* pstr)
// {
//     std::ofstream file_p;

//     file_p.open(filename_p);

//     file_p << pstr ;
//     file_p <<"\n";
//     file_p.close();
// }


// -------- CASE INSENSITIVE COMPARISON 
// Pred = predicate function i.e. a function that :
// 1. returns a boolean output, based on if the inputs satisfy a certain criteria;
// 2. can take in iterators as input

bool compare_pred(unsigned char a, unsigned char b)
{
    return std::tolower(a) == std::tolower(b);
}

bool compare_in(std::string const& a, std::string const& b)
{

    // !!!!!!!!!!!!!!!!!!!!!ADD CODE : Remove any white spaces from start and end of line


    if (a.length()==b.length()) {
        return std::equal(b.begin(), b.end(),
                           a.begin(), compare_pred);
    }
    else {
        return false ;
    }
}

//=========== Read Character history line ===================================




//=========== Read numbers from a string ===================================

int ReadNumbers( const std::string & s, VecD & first ) 
{
    std::istringstream is( s );  // Input stream class to operate on strings
    double n;

    while( is >> n ) 
    {  // The >> is an extraction operator and parses input and places into the specified variable.
                        // It has a boolean return type
        first.push_back( n ); // Adds new element to the end of the vector.
    }
    return first.size();    // Returns the size of the final vector.
}


int ReadIntNumbers( const std::string & s, VecI & first ) 
{
    std::istringstream is( s );  // Input stream class to operate on strings
    int n;

    while( is >> n ) 
    {  // The >> is an extraction operator and parses input and places into the specified variable.
                        // It has a boolean return type
        first.push_back( n ); // Adds new element to the end of the vector.
    }
    return first.size();    // Returns the size of the final vector.
}


//=========================================== PHASE FILE FUNCTIONS ===========================================

void read_phase(const char* filename_p,settings &set0)
{
    std::ifstream file_p;
    std::string line;
    file_p.open(filename_p);

    if(file_p.is_open())
    {
        // File exists. Start read in loop

        while(std::getline(file_p,line))
        {   



            // read current state of the smf loop
            if(compare_in(line,"Phase")) 
            {
                std::cout<<line<<"\n";                
                std::getline(file_p,line);

                if (compare_in(line,"START"))
                    set0.phase = 0;
                else if (compare_in(line,"SEARCH"))
                    set0.phase = 2;
                else if (compare_in(line,"POLL-FEAS"))
                    set0.phase = 3;
                else if (compare_in(line,"POLL-INF"))
                    set0.phase = 4;
                else if (compare_in(line,"LHS"))
                    set0.phase = 1;
                else if (compare_in(line,"RESTART"))
                    set0.phase = 5;
                else
                {
                    std::cout<<"Incorrect phase option provided \nPlease select from : LHS, Search, Poll or Restart\n";
                    throw std::runtime_error("Unknown phase option selected \n");
                }

                std::cout<<line<<"\n";

            }

            // Read dimension and target number of LHS points
            if(compare_in(line,"Dimension"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.dim = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }


            if(compare_in(line,"Number of LHS divisions"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.ndiv = std::stoi(line);
                
                std::cout<<line<<"\n";
            }                        



            // read bounds
            if(compare_in(line,"Upper bound"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.xupp);
                
                std::cout<<line<<"\n";
            }

            if(compare_in(line,"Lower bound"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.xlow);
                
                std::cout<<line<<"\n";
            }

            // Read spacing 

            if(compare_in(line,"Spacing"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadIntNumbers(line,set0.spc);

                std::cout<<line<<"\n";
            }

            // Read number of points at each phase

            if(compare_in(line,"Number of lhs points"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.npts_lhs = std::stoi(line);
                 
                std::cout<<line<<"\n";
            }

            if(compare_in(line,"Number of search points"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.npts_search = std::stoi(line);
                
                std::cout<<line<<"\n";            }

            if(compare_in(line,"Number of feasible poll points"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.npts_pollfeas = std::stoi(line);
                
                std::cout<<line<<"\n";
            }            

            if(compare_in(line,"Number of infeasible poll points"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.npts_pollinf = std::stoi(line);
                
                std::cout<<line<<"\n";
            }

            // Read in number of search and poll stages executed till now
            if(compare_in(line,"Number of search stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.nstage_search = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }

            if(compare_in(line,"Number of feasible poll stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.nstage_pollfeas = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }


            if(compare_in(line,"Number of infeasible poll stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.nstage_pollinf = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }

            // Read in number of successful search and poll stages till now

            if(compare_in(line,"Number of successful search stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.succ_search = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }

            if(compare_in(line,"Number of successful feasible poll stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.succ_pollfeas = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }

            if(compare_in(line,"Number of successful infeasible poll stages"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line); 
                set0.succ_pollinf = std::stoi(line);
                
                std::cout<<line<<"\n";            
            }

            // Read in initial guess for hyper-parameter values

            if(compare_in(line,"Kriging theta"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.theta);
                
                std::cout<<line<<"\n";            
            }

            // Read in hyper-parameter bounds

            if(compare_in(line,"Kriging theta upper bound"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.thupp);
                
                std::cout<<line<<"\n";            
            }


            if(compare_in(line,"Kriging theta lower bound"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.thlow);
                
                std::cout<<line<<"\n";            
            }


            // Read in pointwise stage history

            if(compare_in(line,"Pointwise history"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                // Assume a maximum number of 5000 points stored in history
 //               std::size_t len_hist = line.copy(set0.hist_id,5000);
                set0.pt_hist = line ;               
                std::cout<<line<<"\n";            
            }

            if(compare_in(line,"IDs of points that were successful"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadIntNumbers(line,set0.succ_pts_id);

                std::cout<<line<<"\n";
            }

            if(compare_in(line,"Function values reduced"))
            {
                std::cout<<line<<"\n";
                std::getline(file_p,line);
                ReadNumbers(line,set0.succ_delta);
                
                std::cout<<line<<"\n";            
            }



        }
    }
    else
    {        
        std::cout<<"No phase file present \n ";

        // Return error
        throw std::runtime_error("Need a phase file for execution \n");

   }

    file_p.close();

}

void write_phase(const char* filename_p, settings set0)
{
    FILE *file_p;    // Input file stream class to operate on files
    file_p = fopen(filename_p,"w");

    // Write out all the options into the file

    // Phase
    fprintf(file_p,"Phase\n");
    if(set0.phase == 1)
        fprintf(file_p,"LHS\n");
    else if (set0.phase == 2)
        fprintf(file_p,"SEARCH\n");
    else if (set0.phase == 3)
        fprintf(file_p,"POLL-FEAS\n");
    else if (set0.phase == 4)
        fprintf(file_p,"POLL-INF\n");

    // Dimension 
    fprintf(file_p,"Dimension\n");
    fprintf(file_p,"%d\n",set0.dim);    

    // Target number of LHS points 
    fprintf(file_p,"Number of LHS divisions\n");
    fprintf(file_p,"%d\n",set0.ndiv);    


    // Effective hyper-parameter dimension
//    fprintf(file_p,"Effective hyper-parameter dimensions\n");
 //   fprintf(file_p,"%d\n",set0.dim_th);    



    // Upper and lower bounds
    fprintf(file_p,"Upper bound\n");
    // Use iterators here
    // std::vector<double>::iterator iter; 


    for(int i=0 ; i< set0.xupp.size(); i++ )
    {
        fprintf(file_p,"%.8f ",set0.xupp[i]);
    }
    fprintf(file_p, "\n");

    fprintf(file_p,"Lower bound\n");

    for(int i=0 ; i< set0.xlow.size(); i++ )
    {
        fprintf(file_p,"%.8f ",set0.xlow[i]);
    }
    fprintf(file_p, "\n");

    // Kriging hyper-parameters

    fprintf(file_p,"Kriging theta\n");

    for(int i=0 ; i< set0.theta.size(); i++ )
    {
        fprintf(file_p,"%.8f ",set0.theta[i]);
    }

    fprintf(file_p, "\n");


    fprintf(file_p,"Kriging theta upper bound\n");

    for(int i=0 ; i< set0.theta.size(); i++ )
    {
        fprintf(file_p,"%.8f ",set0.thupp[i]);
    }

    fprintf(file_p, "\n");



    fprintf(file_p,"Kriging theta lower bound\n");

    for(int i=0 ; i< set0.theta.size(); i++ )
    {
        fprintf(file_p,"%.8f ",set0.thlow[i]);
    }

    fprintf(file_p, "\n");


    fprintf(file_p,"Spacing\n");
    // // Use iterators here
    // std::vector<int>::iterator iterspc; 


    // std::cout<<"\n\nPrint check point \n ";

    // VecI spcprint = set0.spc; 

    // std::cout<<spcprint[*spcprint.begin()]<<"  \n";

    for(int i =0; i< set0.spc.size() ; i++)
    {
        fprintf(file_p,"%d ",set0.spc[i]);
    }
    fprintf(file_p, "\n");



    // SMF operational statistics


    fprintf(file_p,"\n\nSMF statistics : \n\n"); 


    fprintf(file_p,"Number of lhs points\n");
    fprintf(file_p,"%d\n",set0.npts_lhs);

    fprintf(file_p,"Number of search points\n");
    fprintf(file_p,"%d\n",set0.npts_search);

    fprintf(file_p,"Number of feasible poll points\n");
    fprintf(file_p,"%d\n",set0.npts_pollfeas);   

    fprintf(file_p,"Number of infeasible poll points\n");
    fprintf(file_p,"%d\n",set0.npts_pollinf);   

    fprintf(file_p,"Number of search stages\n");
    fprintf(file_p,"%d\n",set0.nstage_search);

    fprintf(file_p,"Number of infeasible poll stages\n");
    fprintf(file_p,"%d\n",set0.nstage_pollinf);   

    fprintf(file_p,"Number of feasible poll stages\n");
    fprintf(file_p,"%d\n",set0.nstage_pollfeas);   

    fprintf(file_p,"Number of successful search stages\n");
    fprintf(file_p,"%d\n",set0.succ_search);

    fprintf(file_p,"Number of successful feasible poll stages\n");
    fprintf(file_p,"%d\n",set0.succ_pollfeas);   

    fprintf(file_p,"Number of successful infeasible poll stages\n");
    fprintf(file_p,"%d\n",set0.succ_pollinf);   


    // Pointwise history

    if(set0.pt_hist.length()!=0)
    {
        fprintf(file_p,"Pointwise history\n");
        fprintf(file_p,"%s\n",set0.pt_hist.c_str());
    }

    // Successful step history

    if(set0.succ_pts_id.size()!=0)
    {
        fprintf(file_p,"IDs of points that were successful\n");
        for(int i =0; i< set0.succ_pts_id.size() ; i++)
        {
            fprintf(file_p,"%d ",set0.succ_pts_id[i]);
        }
        fprintf(file_p, "\n");

        fprintf(file_p,"Function values reduced\n");
        for(int i =0; i< set0.succ_delta.size() ; i++)
        {
            fprintf(file_p,"%.8f ",set0.succ_delta[i]);
        }
        fprintf(file_p, "\n");

    }


    fclose(file_p);
}



// ======== ROUTINES FOR READING AND WRITING FILTER FILES ===============

void import_filter_from_file(const char* filename_X, filter &filt, int &n_pts, int dim, VecI xindex)
{
    VecD first_line;  // Vector that is read in and added into first line 
    std::ifstream file_X;    // Input file stream class to operate on files
    std::string line;        // String class to operate on strings
    int MAX_INT = 32766;
    int cols_check;
    MatD v ; 
    file_X.open(filename_X);
    if (file_X.is_open())
    {
        int i=0;
        std::getline(file_X, line);  // Get line from the file


        cols_check =ReadNumbers( line, first_line ); // Read first line 
        if(dim!=cols_check-2)
        {
            std::cout<<"Input error - dimensions donot match \n";
            return;
        }
        v.push_back(first_line);

        for ( i=1;i<MAX_INT+1;i++)
        {         // Read lines continuously till we hit the end of the input or MAX_INT
            if ( !std::getline(file_X, line) ) break;
                VecD temp;      // Temporary vector to be read 
            ReadNumbers( line, temp );     // Read numbers from the line
            v.push_back(temp);
        }

        n_pts=i;
        if(n_pts >MAX_INT) std::cout<< "Number of rows must be smaller than MAX_INT";

        file_X.close();

        for(i=0;i<dim;i++)
        {
                if(xindex[i]==1)
                {
                        for(int r=0;r<n_pts;r++)
                        {
                                v[r][i] = log10(v[r][i]);
                        }
                }
        }

    }
    else
    {
        std::cout << "file open failed";
        std::cout<<"\n File name : "<<filename_X<<"\n";
    }

// Transfer contents from v to respective filter handles 

	for(int p = 0; p<n_pts; p++)
	{
		VecD tempx ; 
		for(int q= 0 ; q<dim; q++)	
			tempx.push_back(v[p][q]); 
		
        filt.x.push_back(tempx);
		filt.y.push_back(v[p][dim]); 
		filt.h.push_back(v[p][dim+1]);
	}

}





void print_filter_to_file(const char* filename_X, filter filt, int dim,int n_pts, VecI xindex)
{
    FILE *file_X;
    file_X = fopen(filename_X,"w");

    for(int j=0;j<n_pts;j++)
    {	
    	for(int k=0; k<dim;k++)
	{
            if(xindex[k] == 1)
			     fprintf(file_X,"%.6f	",pow(10,filt.x[j][k]));
            else
                 fprintf(file_X,"%.6f	",filt.x[j][k]);

	}
        fprintf(file_X,"%.6f	%.6f\n",filt.y[j],filt.h[j]);

    }
    fclose(file_X);
}




// ROUTINES TO READ AND WRITE VECTORS OR MATRICES FROM OTHER FILES  

void import_matrix_from_txt_file(const char* filename_X, MatD & v, int & rows, int cols, VecI xindex)
{
    VecD first_line;  // Vector that is read in and added into first line 
    std::ifstream file_X;    // Input file stream class to operate on files
    std::string line;        // String class to operate on strings
    int MAX_INT = 32766;
    int cols_check;
    file_X.open(filename_X);
    if (file_X.is_open())
    {
        int i=0;
        std::getline(file_X, line);  // Get line from the file
        
        
        cols_check =ReadNumbers( line, first_line ); // Read first line 
        if(cols!=cols_check)
        {
            std::cout<<"Input error - dimensions donot match \n";
            return;
        }
	    v.push_back(first_line);
        
        for ( i=1;i<MAX_INT+1;i++)
        {         // Read lines continuously till we hit the end of the input or MAX_INT
            if ( !std::getline(file_X, line) ) break;
	        VecD temp;  	// Temporary vector to be read 
            ReadNumbers( line, temp );     // Read numbers from the line
            v.push_back(temp);
        }
        
        rows=i; 
        if(rows >MAX_INT) std::cout<< "Number of rows must be smaller than MAX_INT";
        
        file_X.close();

	for(i=0;i<cols;i++)
	{
		if(xindex[i]==1)
		{
			for(int r=0;r<rows;r++)
			{
				v[r][i] = log10(v[r][i]);
			}
		}
	}
    }
    else
    {
        std::cout << "file open failed";
        std::cout<<"\n File name : "<<filename_X<<"\n";
    }

}

void import_vector_from_txt_file(const char* filename_Y, VecD & v, int& num_pts)
{
    std::ifstream file_Y;    // Input file stream class to operate on files
    std::string line;        // String class to operate on strings
    file_Y.open(filename_Y);
    if (file_Y.is_open())
    {
        num_pts = 0;
        while(std::getline(file_Y, line))
        {
                std::istringstream is( line );
                // is = Input stream class to operate on strings
                double n;
                while( is >> n )
                {
                        v.push_back( n );
                        num_pts++;
                }
        }
//	std::getline(file_Y, line);
//	num_pts = ReadNumbers(line,v);
        file_Y.close();
    }
    else
    {
        std::cout << "file open failed";
        std::cout<<"\n File name : "<<filename_Y<<"\n";
    }

}

void print_matrix_to_txt_file(const char* filename_X, MatD v, int num_pts, int dim,VecI xindex)
{
    FILE *file_X;    // Input file stream class to operate on files
    file_X = fopen(filename_X,"w");
    for(int i=0;i<num_pts;i++)
    {    
        for(int j=0;j<dim;j++)
	{
	    if(xindex[j] == 1)
            	fprintf(file_X,"%f  ",pow(10,v[i][j]));
	    else
		fprintf(file_X,"%f  ",v[i][j]);	
	}
        fprintf(file_X,"\n");
    }
        
    fclose(file_X);

}


void print_vector_to_txt_file(const char* filename_X, VecD v, int dim)
{
    FILE *file_X;    
    file_X = fopen(filename_X,"w");

    for(int j=0;j<dim;j++)
        fprintf(file_X,"%f\n",v[j]);
        
    fclose(file_X);
}

void write_hist(FILE *fp, VecD x, double J,VecI xindex)
{
    int n = x.size();
    fprintf(fp,"%f        ",J);
    for(int i=0;i<n;i++)
    {
	if(xindex[i]==1)
        	fprintf(fp,"%f  ",pow(10,x[i]));
	else
                fprintf(fp,"%f  ",x[i]);
	
    }
    fprintf(fp,"\n");
}

void print_theta_hist(model S)
{

	FILE *fp; 
	fp = fopen("theta_hist","a");

	for(int i=0;i<S.dim_th;i++)
		fprintf(fp,"%f  ",S.theta[i]);

	fprintf(fp,"%d \n",S.num_pts);


	fclose(fp);  

}

// ====================================================================================================
