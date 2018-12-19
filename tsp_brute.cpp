#include "tsp_incl.h"

void find_path(int low, int high, double &disti)
{
	if(disti > min_dist)
		return;

	if(low==high)
	{
		disti += distances[idx[low-1]][idx[low]];
		disti += distances[idx[low]][idx[0]];
														// if instead of this distance at runtime, if we calculate distance at end, then we get
		if(min_dist>=disti)							// a bad timing as pruning can't be utilised and also will require O(n) at every permutation 
		{												// instead of just minimum ones											
			for(int i=0;i<num;i++)
			{
				cout<<idx[i]<<" ";
				min_path[i] = idx[i];
			}
			cout<<endl;
			min_dist = disti;
		}

		disti -= distances[idx[low-1]][idx[low]];
		disti -= distances[idx[low]][idx[0]];
	}

	for(int i=low;i<=high;i++)
	{
		swap(idx[low],idx[i]);
		disti += distances[idx[low-1]][idx[low]];

		find_path(low+1, high, disti);

		disti -= distances[idx[low-1]][idx[low]];
		swap(idx[low],idx[i]);
	}
}


int main()
{
	char cities[] = "datasets/ex2.txt";		// use only ex1,ex2,ex3 for this brute force algorithm - Bigger examples are not feasible

	// if distances directly provided without coordinates use follow_dist function
	// --------------------------
	// follow_dist(cities);
	// --------------------------
	// else run these commands
	// ---------------------------
	store_points(cities);
	calc_distances(0);
	// --------------------------

	// print_distances();
	
	double disti = 0;				// just to pass as reference

	min_dist = INT_MAX;
	min_path = vect<int> (num);

	clock_t start=clock();

	find_path(1,num-1,disti);		// complexity = O(n!) - as need to traverse array (while storing new minimums)
	
	clock_t end=clock();

	cout<<endl<<"Minimum distance = "<<min_dist<<endl;

	cout<<"Minimum Path Sequence : ";
	for(int i=0;i<num;i++)
		cout<<min_path[i]<<" ";
	cout<<endl;

	cout<<"Time taken for solution = "<<((end-start)/double(CLOCKS_PER_SEC))<<endl;
	return 0;
}