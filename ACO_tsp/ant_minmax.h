// header file for functions of min-max ant system heuristic algorithm with 3opt optimisations for tour improvements

void improve_tour(double &ant_d, vect<int> &ant_p)
{
	vect<int> t_p = min_path;
	double t_d = min_dist;

	min_dist = ant_d;
	min_path = ant_p;

	update_positions();
	// vector<int> temp = positions;

	// sort(temp.begin(),temp.end());
	// for(int i=0;i<num;i++)
	// 	cout<<temp[i]<<" ";
	// cout<<endl<<endl;

	improve_3opt_nn_dlb();

	ant_d = min_dist;
	ant_p = min_path;

	min_path = t_p;
	min_dist = t_d;
}

void perform_ant_minmax()
{
	alpha = 1; beta = 4; rho = 0.8;
	m = 25; n_iter = 250; nn = 20;

	tau_mx = (1.0/(1-rho))*(1/min_dist);
	tau_mi = tau_mx/(2*num);				// corresponds to p_best = 0.005 = probability of constructing best solution after MMAS defined convergence
	
	// general formula for tau_mi = tau_mx * ((1 - c)/((avg - 1)*c)) , where avg = num/2 and c = (p_best)^(1/num)

	srand(0);

	ant_paths = vect<vect<int> > (m, vect<int> (num,0));
	pher = vect<vect<double> > (num,vect<double> (num,tau_mx*10));
	prob_nums = vect<vect<double> > (num, vect<double> (num,0));
	
	calc_visib(beta);
	min_dist = INT_MAX;   // so that update rule dont get biased using earlier nn opt solution


	for(int t=1;t<=n_iter;t++)
	{
		iter_dist = INT_MAX;
		int s = 0, flag = 0;

		delp = vect<vect<double> > (num,vect<double> (num,0));
		tabus = vect<vect<bool> > (m, vect<bool> (num,0)); 
		ant_dist = vect<double> (m,0);
		ncount = vect<int> (m,0);

		for(int i=0;i<m;i++)
		{
			int posn = rand()%num;
			ant_paths[i][0] = posn;		// placed m ants on n nodes randomly
			tabus[i][ant_paths[i][s]]=1;
		}

		s++;
		calc_probab_numes(alpha);

		while(s<num)
		{
			for(int i=0;i<m;i++)
			{
				// calculate next town

				int prev = ant_paths[i][s-1];
				int next = -1; ncount[i] = 0;				// use from nearest 20 of prev if available

				double deno = 0;

				for(int k=0;k<nn;k++)
				{
					int j = nearest_neighbours[prev][k];

					if(tabus[i][j]==0)
					{
						deno += prob_nums[prev][j];
						ncount[i]++;
					}
				}

				if(ncount[i]>0)
				{
					vector<double> probabs(nn,0);

					for(int k=0;k<nn;k++)
					{
						int j = nearest_neighbours[prev][k];
						
						if(tabus[i][j]==0)
						{
							probabs[k] = prob_nums[prev][j]/deno;
						}
						probabs[k] += (k-1>=0 ? probabs[k-1] : 0);
					}

					double val = (double) rand()/RAND_MAX;

					for(int k=0;k<nn;k++)
					{
						int j = nearest_neighbours[prev][k];
						if(val < probabs[k])
						{
							next = j;
							break;
						}	
					}
				}

				else
				{
					double maxi = 0;

					for(int j=0;j<num;j++)
					{
						if(tabus[i][j] == 0 && maxi < prob_nums[prev][j])			// initialising next, based on max probab first
						{
							maxi = prob_nums[prev][j];
							next = j;
						}	
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

			improve_tour(ant_dist[j], ant_paths[j]);	// local search - 3opt
			cout<<ant_dist[j]<<endl;
			
			if(ant_dist[j] < iter_dist)
			{
				iter_dist = ant_dist[j];
				iter_path = ant_paths[j];
			}
		}

		if(iter_dist < min_dist)
		{
			min_dist = iter_dist;
			min_path = iter_path;
		}
		cout<< min_dist<<endl;
		

		// tau limit updates afer every iteration
		
		tau_mx = (1.0/(1-rho))*(1/min_dist);
		tau_mi = tau_mx/(2*num);


		// updation using iteration best method

		for(int i=0;i<num;i++)
		{
			int s = iter_path[i], e = iter_path[(i+1)%num];
			delp[s][e] = 1/iter_dist;
			delp[e][s] = delp[s][e]; 	// choice here
		}

		for(int i=0;i<num;i++)
		{
			for(int k=0;k<nn;k++)
			{
				int j = nearest_neighbours[i][k];
				pher[i][j] = pher[i][j]*rho + delp[i][j];

				if(pher[i][j] < tau_mi) pher[i][j] = tau_mi;		// checking limits
				else if(pher[i][j] > tau_mx) pher[i][j] = tau_mx;

				pher[j][i] = pher[i][j];
			}
		}



		// increasing frequency of global update rule as iter increases

		if(t>25)
		{
			if(t<=75)
			{if(t%5==0) flag=1;}

			else if(t<=125)
			{if(t%3==0) flag=1;}

			else if(t<=250)
			{if(t%2==0) flag=1;}

			else if(t>250) flag = 1;
		}
		

		if(flag==1) // pheromone trail update using global best
		{
			delp = vect<vect<double> > (num,vect<double> (num,0));

			for(int i=0;i<num;i++)
			{
				int s = min_path[i], e = min_path[(i+1)%num];
				delp[s][e] += 1/min_dist;
				delp[e][s] = delp[s][e]; 	// choice here
			}

			for(int i=0;i<num;i++)
			{
				for(int k=0;k<nn;k++)
				{
					int j = nearest_neighbours[i][k];
					pher[i][j] = pher[i][j]*rho + delp[i][j];

					if(pher[i][j] < tau_mi) pher[i][j] = tau_mi;		// checking limits
					else if(pher[i][j] > tau_mx) pher[i][j] = tau_mx;

					pher[j][i] = pher[i][j];
				}
			}
		}

	}
}