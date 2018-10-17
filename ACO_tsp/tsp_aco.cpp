#include "tsp_incl.h"
#include "tsp_aco_incl.h"
#include "ant_system.h"
#include "ant_colony.h"
#include "ant_minmax.h"

int main()
{
	char cities[] = "datasets/ex8.txt";

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

	prepare_neighbours();					// constructs nearest neighbour lists - O(n^2.log(n))
	construct_nn(1);

	// perform_ant_system();
	// perform_ant_colony();
	perform_ant_minmax();

	clock_t end=clock();

	cout<<"Minimum distance = "<<min_dist<<endl;

	cout<<"Minimum Path Sequence : ";
	
	for(int i=0;i<num;i++)
		cout<<min_path[i]<<" ";
	cout<<endl;

	cout<<"Time Taken for solution = "<<((end-start)/double(CLOCKS_PER_SEC))<<endl;
	return 0;
}