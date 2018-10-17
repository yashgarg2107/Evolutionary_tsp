#include "tsp_2opt.h"
#include "tsp_3opt.h"

int n_iter, m;
double alpha = 1, beta = 1;		// default parameters
double Q = 1, rho = 0.5;		// rho - evaporation constant
double tau0, q0 = 0.5;

vect<vect<double> > prob_nums;
vect<vect<double> > pher;
vect<vect<double> > delp;
vect<vect<double> > visib;
vect<vect<int> > ant_paths;
vect<vect<bool> > tabus;
vect<double> ant_dist;

vector<int> ncount;		// for ants
vect<int> iter_path;
double iter_dist;

double tau_mx, tau_mi;	// for min max optimisation

void calc_visib(double beta)
{
	visib = vect<vect<double> > (num, vect<double> (num,0));

	for(int i=0;i<num;i++)
	{
		for(int j=i+1;j<num;j++)
		{
			if(distances[i][j]==0)		// to handle cases having same coordinates as in a280 dataset(170-171 pair)
			{
				visib[i][j]=1;
				visib[j][i]=1;
				continue;
			}
			visib[i][j] = 1/distances[i][j];
			visib[j][i] = pow(visib[i][j],beta);
			visib[i][j] = visib[j][i];
		}
	}
}

void calc_probab_numes(double alpha)
{
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			prob_nums[i][j] = pow(pher[i][j],alpha) * visib[i][j];
		}
	}
}