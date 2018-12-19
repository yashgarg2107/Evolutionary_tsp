// header file for crossover operators of genetic crossover approach for tsp

void pmx_crossover()	// partially mapped crossover
{
	vect<gene> next_pop;
	int i = 0;

	while(next_pop.size()<popu)
	{
		gene par1 = population[next_par[i].first];
		gene par2 = population[next_par[i].second];

		gene child1, child2;
		child1.path = vect<int> (num,-1);
		child2.path = vect<int> (num,-1);

		vect<bool> vis1(num, 0);
		vect<bool> vis2(num, 0);

		int p1 = rand()%num; 
		int p2 = p1;
		
		while(p2==p1)
			p2 = rand()%num;

		if(p1>p2)
			swap(p1,p2);

		for(int j=p1;j<=p2;j++)
		{
			child1.path[j] = par1.path[j];
			vis1[child1.path[j]] = 1;

			child2.path[j] = par2.path[j];
			vis2[child2.path[j]] = 1;
		}

		for(int j=p1;j<=p2;j++)
		{
			if(vis1[par2.path[j]]==0)
			{
				int val = par2.path[j];
				int p = par2.posns[par1.path[j]];

				while(child1.path[p]!=-1)
				{
					int temp = child1.path[p];
					p = par2.posns[par1.path[p]];
				}

				child1.path[p] = val;
				vis1[val] = 1;
			}

			if(vis2[par1.path[j]]==0)
			{
				int val = par1.path[j];
				int p = par1.posns[par2.path[j]];

				while(child2.path[p]!=-1)
				{
					int temp = child2.path[p];
					p = par1.posns[par2.path[p]];
				}

				child2.path[p] = val;
				vis2[val] = 1;
			}
		}

		for(int j=0;j<num;j++)
		{
			if(child1.path[j]==-1)				// copy the remaining of par2's path
				child1.path[j] = par2.path[j];

			if(child2.path[j]==-1)
				child2.path[j] = par1.path[j];
		}

		i++;
		next_pop.add(child1);		// if you want only one child from one pair of parent chromosomes comment below 3 lines

		if(next_pop.size()==popu)
			break;
		next_pop.add(child2);
	}
	population = next_pop;
}


void ox_crossover()	// ordered crossover
{
	vect<gene> next_pop;
	int i = 0;

	while(next_pop.size()<popu)
	{
		gene par1 = population[next_par[i].first];
		gene par2 = population[next_par[i].second];

		gene child1, child2;
		child1.path = vect<int> (num,-1);
		child2.path = vect<int> (num,-1);

		vect<bool> vis1(num, 0);
		vect<bool> vis2(num, 0);

		int p1 = rand()%num; 
		int p2 = p1;
		
		while(p2==p1)
			p2 = rand()%num;

		if(p1>p2)
			swap(p1,p2);

		for(int j=p1;j<=p2;j++)
		{
			child1.path[j] = par1.path[j];
			vis1[child1.path[j]] = 1;

			child2.path[j] = par2.path[j];
			vis2[child2.path[j]] = 1;
		}

		p1 = 0; p2 = 0;	// redefining role of p1 and p2
		
		while(p1<num && child1.path[p1]!=-1)
			p1++;
		while(p2<num && child2.path[p2]!=-1)
			p2++;

		for(int j=0;j<num;j++)
		{
			if(vis1[par2.path[j]]==0)				// copy the remaining of par2's path
			{
				child1.path[p1] = par2.path[j];
				while(p1<num && child1.path[p1]!=-1)
					p1++;
			}

			if(vis2[par1.path[j]]==0)
			{
				child2.path[p2] = par1.path[j];
				while(p2<num && child2.path[p2]!=-1)
					p2++;
			}
		}

		i++;
		next_pop.add(child1);		// if you want only one child from one pair of parent chromosomes comment below 3 lines
		
		if(next_pop.size()==popu)
			break;
		next_pop.add(child2);
	}
	population = next_pop;
}


void vgx_crossover()	// very greedy crossover
{
	vect<gene> next_pop;

	for(int i=0;i<next_par.size();i++)
	{
		gene par1 = population[next_par[i].first];
		gene par2 = population[next_par[i].second];

		gene child;
		vect<bool> vis(num, 0);

		int p1, p2, ch;
		int n1l, n1r, n2l, n2r, nx;

		p1 = rand()%num;
		ch = par1.path[p1];
		child.path.add(ch);
		vis[ch] = 1;

		while(child.path.size()!=num)
		{
			nx = -1;

			ch = child.path[child.path.size()-1];

			p1 = par1.posns[ch];
			n1l = par1.path[(p1+num-1)%num];
			n1r = par1.path[(p1+1)%num];
			
			p2 = par2.posns[ch];
			n2l = par2.path[(p2+num-1)%num];
			n2r = par2.path[(p2+1)%num];

			vect<pair<double,int> > near_pts;
			near_pts.add(make_pair(distances[ch][n1l],n1l));
			near_pts.add(make_pair(distances[ch][n1r],n1r));
			near_pts.add(make_pair(distances[ch][n2l],n2l));
			near_pts.add(make_pair(distances[ch][n2r],n2r));

			sort(near_pts.begin(),near_pts.end());
			
			for(int j=0;j<4;j++)
			{
				if(!vis[near_pts[j].second])
				{
					nx = near_pts[j].second;
					break;
				}
			}

			if(nx==-1)
			{
				for(int j=0;j<num-1;j++)
				{
					if((!vis[nearest_neighbours[ch][j]]) && (ch!=nearest_neighbours[ch][j]))
					{
						nx = nearest_neighbours[ch][j];
						break;
					}
				}
			}

			child.path.add(nx);
			vis[nx] = 1;
		}

		next_pop.add(child);
	}
	population = next_pop;
}


void uhx_crossover()	// unnamed heuristic crossover
{						// was making a big mistake in this function, was incrementing the city number instead of city pointers
	vect<gene> next_pop;

	for(int i=0;i<next_par.size();i++)
	{
		gene par1 = population[next_par[i].first];
		gene par2 = population[next_par[i].second];

		gene child;
		vect<bool> vis(num, 0);

		int p1r,p1l,p2l, p2r;		// positional indexes
		int n1l, n1r, n2l, n2r;		// city numbers
		int ch, nx;

		p1r = rand()%num; 
		ch = par1.path[p1r];
		p2r = par2.posns[ch];

		child.path.add(ch);
		vis[ch] = 1;

		p1l = p1r; p2l = p2r;

		p1l = (p1l+num-1)%num;
		p1r = (p1r+1)%num;
		p2l = (p2l+num-1)%num;
		p2r = (p2r+1)%num;

		n1l = par1.path[p1l];
		n1r = par1.path[p1r];
		n2l = par2.path[p2l];
		n2r = par2.path[p2r];

		while(child.path.size()!=num)
		{
			nx = -1;
			ch = child.path[child.path.size()-1];

			vect<pair<double,int> > near_pts;
			near_pts.add(make_pair(distances[ch][n1l],n1l));
			near_pts.add(make_pair(distances[ch][n1r],n1r));
			near_pts.add(make_pair(distances[ch][n2l],n2l));
			near_pts.add(make_pair(distances[ch][n2r],n2r));

			sort(near_pts.begin(),near_pts.end());
			
			for(int j=0;j<4;j++)
			{
				if(!vis[near_pts[j].second])
				{
					int v = near_pts[j].second;

					if(v==n1l)
						nx = n1l; 
					else if(v==n1r)
						nx = n1r; 
					else if(v==n2l)
						nx = n2l; 
					else if(v==n2r)
						nx = n2r; 

					break;
				}
			}

			if(nx!=-1)
			{
				child.path.add(nx);
				vis[nx] = 1;
			}

			if(child.path.size()==num)
				break;

			while(vis[n1l])
			{
				p1l = (p1l+num-1)%num;
				n1l =  par1.path[p1l];
			}
			while(vis[n1r])
			{
				p1r = (p1r+1)%num;
				n1r = par1.path[p1r];
			}
			while(vis[n2l])
			{
				p2l = (p2l+num-1)%num;
				n2l =  par2.path[p2l];
			}
			while(vis[n2r])
			{
				p2r = (p2r+1)%num;
				n2r = par2.path[p2r];
			}
		}

		next_pop.add(child);	
	}
	population = next_pop;
}

// Reference : https://arxiv.org/pdf/1504.02590.pdf