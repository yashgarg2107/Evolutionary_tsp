#include "tsp_incl.h"
#include "tsp_gen_incl.h"
#include "tsp_3opt.h"
#include "gen_selections.h"
#include "gen_crossovers.h"
#include "gen_mutations.h"

int main()
{
	char cities[] = "datasets/ex5.txt";

	// if distances directly provided without coordinates use follow_dist function
	// --------------------------
	// follow_dist(cities);
	// --------------------------
	// else run these commands
	// ---------------------------
	store_points(cities);
	calc_distances(0);
	// --------------------------
	
	clock_t start=clock();

	srand(0);
	construct_random();	
	min_dist = INT_MAX/2;

	prepare_neighbours();

	popu = 100; elit = 10;
	population = vect<gene> (popu);
	next_par = vect<pair<int,int> > (popu, pair<int,int>());
	
	construct_random_pop();
	calc_fitness();

	perform_elitism();
	roulette_select();

	for(int i=0;i<400;i++)
	{
		cout<<min_dist<<endl;

		update_positions_gen();

		uhx_crossover();
		mutate();

		if(i>200)			// can use other custom strategies for vgx, uhx heuristic crossovers	
			fill_elites();		// for pmx and ox crossover elitism is very important, so use it for every iteration

		calc_fitness();

		if(i>200)
			perform_elitism();

		if(i%2==0)
			roulette_select();
		else
			tournament_sel();
	}

	calc_fitness();
	clock_t end=clock();

	cout<<"Minimum distance = "<<min_dist<<endl;

	cout<<"Minimum Path Sequence : ";
	for(int i=0;i<num;i++)
		cout<<min_path[i]<<" ";
	cout<<endl;

	cout<<"Time Taken for solution = "<<((end-start)/double(CLOCKS_PER_SEC))<<endl;
	return 0;
}
