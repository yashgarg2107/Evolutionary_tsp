// header file for functions of ant system heuristic algorithm with 3opt optimisations for tour improvements

void perform_ant_system()
{
	alpha = 1; Q = 1; 		// in later versions by same authors, no alpha and Q parameters, hence default 1
	beta = 4; rho = 0.5;

	m = 25; n_iter = 500;

	srand(0);

	ant_paths = vect<vect<int> > (m, vect<int> (num,0));
	pher = vect<vect<double> > (num,vect<double> (num,0.001));
	prob_nums = vect<vect<double> > (num, vect<double> (num,0));

	calc_visib(beta);
	min_dist = INT_MAX;   // so that update rule dont get biased using earlier nn opt solution
	
	
	for(int t=1;t<=n_iter;t++)
	{
		int s=0;

		delp = vect<vect<double> > (num,vect<double> (num,0));
		tabus = vect<vect<bool> > (m, vect<bool> (num,0)); 
		ant_dist = vect<double> (m,0);

		for(int i=0;i<m;i++)
		{
			int posn = rand()%num;
			ant_paths[i][0] = posn;		// placed m ants on n nodes randomly
			tabus[i][ant_paths[i][s]]=1;	// 1 marks visited node
		}			

		s++;
		calc_probab_numes(alpha);	// as updates happen only at end of all ants traversal, so we only need to calculate 
									// probability numerators at start of each iteration, this is different from case of ACS

		while(s<num)
		{
			for(int i=0;i<m;i++)
			{
				double deno = 0;
				vector<double> probabs(num,0);

				for(int j=0;j<num;j++)
					if(tabus[i][j]==0)
						deno += prob_nums[ant_paths[i][s-1]][j];

				for(int j=0;j<num;j++)
				{
					if(tabus[i][j]==0)
					{
						probabs[j] = prob_nums[ant_paths[i][s-1]][j]/deno;
					}

					probabs[j] += (j-1>=0 ? probabs[j-1] : 0);
				}

				int next = -1;
				double maxi = 0;
				double val = (double) rand()/RAND_MAX;

				for(int j=0;j<num;j++)
				{					
					if(val < probabs[j])	// linear search working faster than lower bound method of stl
					{
						next = j; 
						break;
					}	
				}

				ant_paths[i][s] = next;
				tabus[i][next] = 1;
			}
			s++;
		}

		for(int j=0;j<m;j++)
		{
			for(int i=0;i<num;i++)
			{
				int s = ant_paths[j][i], e = ant_paths[j][(i+1)%num];
				ant_dist[j] += distances[s][e];
			}

			for(int i=0;i<num;i++)		// updating delta pheromone because of each ant's traversal
			{
				int s = ant_paths[j][i], e = ant_paths[j][(i+1)%num];
				delp[s][e] += Q/ant_dist[j];
				delp[e][s] = delp[s][e];
			}

			if(ant_dist[j] < min_dist)
			{
				min_dist = ant_dist[j];
				min_path = ant_paths[j];
			}
		}
		// cout<<min_dist<<endl;

		for(int i=0;i<num;i++)		 // updating final pheromone trails
		{
			for(int j=0;j<num;j++)
			{
				pher[i][j] = pher[i][j]*rho + delp[i][j];
			}
		}
	}
}