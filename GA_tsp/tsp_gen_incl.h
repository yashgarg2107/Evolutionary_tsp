// header file for important functions of genetic algorithm for solving tsp

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
	for(int i=0;i<popu;i++)
	{
	 	double tdist = 0;
	 	int s,e;

	 	for(int j=1;j<num;j++)
	 	{
	 		s = population[i].path[j-1];
	 		e = population[i].path[j];
	 		tdist += distances[s][e];
	 	}

	 	s = population[i].path[num-1];
	 	e = population[i].path[0];

	 	tdist += distances[s][e];
	 	population[i].dist = tdist;
	 	population[i].fitness = 1.0/tdist;

	 	// if(tdist>2*min_dist)			// pruning strategy
	 	// 	population[i].fitness = 0;

	 	if(population[i].dist < min_dist)
	 	{
	 		min_path = population[i].path;
	 		min_dist = population[i].dist;
	 	}
	}
}