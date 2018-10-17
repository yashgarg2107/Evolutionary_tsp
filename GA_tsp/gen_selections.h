// header file for selection(and elitism) operators of genetic crossover approach for tsp

void roulette_select()
{
	double deno = 0;
	vector<double> probabs(popu,0);

	for(int j=0;j<popu;j++)
		deno += population[j].fitness;

	for(int j=0;j<popu;j++)
	{
		probabs[j] = population[j].fitness/deno;
		probabs[j] += (j-1>=0 ? probabs[j-1] : 0);			// accumulation for roulette wheel selection (cdf basically)
	}

	int k = 0;

	while(k < next_par.size())
	{
		int next = -1;
		double val = (double) rand()/RAND_MAX;

		for(int j=0;j<popu;j++)
		{					
			if(val < probabs[j]){	// linear search working faster than lower bound method of stl
				next = j; 
				break;
			}	
		}

		next_par[k].first = next;

		while(next == next_par[k].first)		// to get a different second parent
		{
			val = (double) rand()/RAND_MAX;
			for(int j=0;j<popu;j++)
			{					
				if(val < probabs[j]){
					next = j; 
					break;
				}	
			}
		}	

		next_par[k].second = next;
		k++;
	}
}


void tournament_sel()		// study reservior sampling also, although not neede here
{
	sort(population.begin(), population.end(),compare);		// take the best parent out of k samples, population sorted for this purpose

	int k = 5; // k-way tournament selection
	vect<int> vals(k);

	for(int i=0;i<next_par.size();i++)
	{
		for(int j=0;j<k;j++)
			vals[j] = rand()%popu;

		sort(vals.begin(),vals.end());
		next_par[i].first = vals[0];

		for(int j=0;j<k;j++)
			vals[j] = rand()%popu;

		sort(vals.begin(),vals.end());
		next_par[i].second = vals[0];		
	}
}

void perform_elitism()
{
	for(int i=0;i<popu;i++)
	{
		if(elites.size()<elit)
		{
			elites.push(population[i]);
		}
		else
		{
			gene temp = elites.top();
			if(compare(population[i],temp))
			{
				elites.pop();
				elites.push(population[i]);
			}
		}
	}

}

void fill_elites()
{
	sort(population.begin(), population.end(), compare);

	vect<gene> els;

	while(!elites.empty())
	{
		els.add(elites.top());
		elites.pop();
	}

	int k=0;
	for(int i=popu-elit;i<popu;i++)
	{
		elites.push(els[k]);
		population[i] = (els[k++]);
	}

}