// header file for important functions of genetic algorithm for solving tsp

#include "tsp_2opt.h"
#include "tsp_3opt.h"

struct gene
{
	double dist;
	double fitness;
	vect<int> path;
	vect<int> posns;
};

bool compare(gene &a, gene &b)
{
	if(a.dist < b.dist)
		return true;
	return false;
}

struct mycompare
{
	bool operator()(gene a, gene b)
	{
		if(a.dist < b.dist)
			return true;
		return false;
	}
};

int popu, elit;
vect<gene> population;
vect<pair<int,int> > next_par;
priority_queue<gene, vect<gene>, mycompare> elites;

int myrandom (int i)
{ 
	return rand()%i;
}

void construct_random_pop()
{
	for(int i=0;i<popu;i++)
	{
		random_shuffle(min_path.begin()+1, min_path.end(), myrandom);	// 1 used to fix the initial vertex, to avoid redundant tours
		population[i].path = min_path;
	}
}

void update_positions_gen()
{
	for(int i=0;i<popu;i++)
	{
		population[i].posns = vect<int> (num);

		for(int j=0;j<num;j++)
		{
			population[i].posns[population[i].path[j]] = j;
		}
	}
}

void calc_fitness()
{
	vect<int> posns(num,0);
	for(int i=0;i<popu;i++)
	{
	 	double todo = (double) rand()/RAND_MAX;
		if(todo < 0.02)
	 		improve_3opt_nn_dlb(population[i].path, posns);
		
		update_dist(population[i].path, population[i].dist);

	 	population[i].fitness = 1.0/population[i].dist;

	 	// if(population[i].dist>2*min_dist)			// pruning strategy
	 	// 	population[i].fitness = 0;

	 	if(population[i].dist < min_dist)
	 	{
	 		min_path = population[i].path;
	 		min_dist = population[i].dist;
	 	}
	}
}