// header file for functions of ant colony system heuristic algorithm with 3opt optimisations for tour improvements

void perform_ant_colony()
{
	beta = 4; rho = 0.9;
	m = 25; n_iter = 100; nn=20;

	tau0 = 1.0/(min_dist);	// nn approx here // works better than (n*dist)^(-1)
	q0 = 0.7;  // exploitation vs exploration tradeoff

	srand(0);

	ant_paths = vect<vect<int> > (m, vect<int> (num,0));
	pher = vect<vect<double> > (num,vect<double> (num,tau0));
	vect<int> posns(num,0);

	calc_visib(beta);
	min_dist = INT_MAX;   // so that global update rule dont get biased using earlier nn opt solution

	while(n_iter--)
	{
		int s=0;

		delp = vect<vect<double> > (num,vect<double> (num,0));
		tabus = vect<vect<bool> > (m, vect<bool> (num,0)); 
		ant_dist = vect<double> (m,0);

		for(int i=0;i<m;i++)			
		{
			int posn = rand()%num;
			ant_paths[i][0] = posn;		// placed m ants on n nodes randomly
			tabus[i][ant_paths[i][s]]=1;
		}			

		s++;

		while(s<num)
		{
			int i = 0;

			for(i=0;i<m;i++)
			{
				double deno = 0;
				int prev = ant_paths[i][s-1];
				vector<double> probabs(num,0);

				for(int j=0;j<num;j++)
				{
					if(tabus[i][j]==0)				// for local updates // prob calculation
					{
						probabs[j] = pher[prev][j] * visib[prev][j];
						deno += probabs[j];
					}
				}

				for(int j=0;j<num;j++)
				{
					if(tabus[i][j]==0)
						probabs[j] = probabs[j]/deno;
				}

				int next = 0;
				double maxi = probabs[0];

				for(int j=1;j<num;j++)
				{
					if(maxi < probabs[j])			// initialising next, based on max probab first
					{
						maxi = probabs[j];
						next = j;					// exploitation choice
					}	

					probabs[j] += probabs[j-1];		// genearting cdf by accumulation
				}	

				double q = (double) rand()/RAND_MAX;

				if(q>q0)		// if such probab then explore based on probabilities
				{
					double val = (double) rand()/RAND_MAX;
					for(int j=0;j<num;j++)
					{
						if(val < probabs[j])
						{
							next = j;
							break;
						}	
					}
				}		

				ant_paths[i][s] = next;
				tabus[i][next] = 1;

				pher[prev][next] = pher[prev][next]*rho + (1-rho)*tau0;
				pher[next][prev] = pher[prev][next];		// local updates
			
				if(s==num-1)
				{
					prev = ant_paths[i][s];		// last step of local update
					next = ant_paths[i][0];
					pher[prev][next] = pher[prev][next]*rho + (1-rho)*tau0;
					pher[next][prev] = pher[prev][next];
				}
			}
			s++;
		}

		for(int j=0;j<m;j++)
		{
			improve_3opt_nn_dlb(ant_paths[j], posns);
			update_dist(ant_paths[j], ant_dist[j]);

			if(ant_dist[j] < min_dist)
			{
				min_dist = ant_dist[j];
				min_path = ant_paths[j];
			}
		}

		cout<<min_dist<<endl;

		for(int i=0;i<num;i++)		// calc del from global maxima
		{
			int s = min_path[i], e = min_path[(i+1)%num];

			delp[s][e] += 1/min_dist;
			delp[e][s] = delp[s][e];
		}

		for(int i=0;i<num;i++)		// global update rule
		{
			for(int j=0;j<num;j++)
			{
				pher[i][j] = pher[i][j]*rho + (1-rho)*delp[i][j];
			}
		}
	}
}