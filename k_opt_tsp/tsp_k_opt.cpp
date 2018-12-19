#include "tsp_incl.h"
#include "tsp_2opt.h"
#include "tsp_3opt.h"

int main()
{
	char cities[] = "datasets/ex4.txt";

	// if distances directly provided without coordinates use follow_dist function
	// --------------------------
	// follow_dist(cities);
	// --------------------------
	// else run these commands
	// ---------------------------
	store_points(cities);
	calc_distances(1);
	// --------------------------
	
	// print_distances();

	clock_t start=clock();

	prepare_neighbours();					// constructs nearest neighbour lists - O(n^2.log(n))
	construct_random();						// random path construction

	// construct_nn_opt(1);				
	// construct_nn(1);						// give either num(to get best nn) or 1(for any random nn) as range
	
	nn=20;
	improve_3opt_nn_dlb(min_path, positions);						// from observation at att48 2opt may not give most optimised 
											// on best nn construction, rather gives on random nn cons
	// improve_2opt_bf(min_path, positions);
	// improve_3opt_ff();
	// improve_3opt_bf();

	update_dist(min_path, min_dist);

	clock_t end=clock();


	cout<<"Minimum distance = "<<min_dist<<endl;

	cout<<"Minimum Path Sequence : ";
	for(int i=0;i<num;i++)
		cout<<min_path[i]<<" ";
	cout<<endl;

	cout<<"Time Taken for solution = "<<((end-start)/double(CLOCKS_PER_SEC))<<endl;
	return 0;
}