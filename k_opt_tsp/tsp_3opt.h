// header file for functions of 3opt local search algorithm and its optimized variations

void make_best_move_pos(vect<int> &curr_path, vect<int> &posns, int best_case, int i, int j, int k)
{
	switch(best_case)
	{
		case 0: reverse_tour_pos(curr_path,posns,(i+1)%num,j);
				break;

		case 1: reverse_tour_pos(curr_path,posns,(i+1)%num,j);
				reverse_tour_pos(curr_path,posns,(j+1)%num,k);
				break;

		case 2: reverse_tour_pos(curr_path,posns,(i+1)%num,j);
				reverse_tour_pos(curr_path,posns,(j+1)%num,k);
				reverse_tour_pos(curr_path,posns,(k+1)%num,i);
				break;
	}
}

void make_best_move(vect<int> &curr_path, int best_case, int i, int j, int k)
{
	switch(best_case)
	{
		case 0: reverse_tour(curr_path,(i+1)%num,j);
				break;

		case 1: reverse_tour(curr_path,(i+1)%num,j);
				reverse_tour(curr_path,(j+1)%num,k);
				break;

		case 2: reverse_tour(curr_path,(i+1)%num,j);
				reverse_tour(curr_path,(j+1)%num,k);
				reverse_tour(curr_path,(k+1)%num,i);
				break;
	}
}

pair<double,int> find_best_case(int x1, int x2, int y1, int y2, int z1, int z2)
{
	double current, new_dis[3], gain[3], best_gain;
	
	current = distances[x1][x2] + distances[y1][y2] + distances[z1][z2]; 

	new_dis[0] = distances[x1][y1] + distances[x2][y2] + distances[z1][z2];
	gain[0] = current - new_dis[0];		// represents cyclic single 2-opt cases 1,2,3

	new_dis[1] = distances[x1][y1] + distances[x2][z1] + distances[y2][z2];
	gain[1] = current - new_dis[1];		// represents cyclic doube 2-opt cases 4,5,6

	new_dis[2] = distances[x1][y2] + distances[y1][z2] + distances[z1][x2];
	gain[2] = current - new_dis[2];		// represents cyclic triple 2-opt case 7

	best_gain = max(gain[0],max(gain[1],gain[2]));

	if(best_gain<=0.0000001)			// to handle precision problems
		return make_pair(0,-1);

	else if(best_gain==gain[0])
		return make_pair(best_gain,0);

	else if (best_gain==gain[1])
		return make_pair(best_gain,1);
	else 
		return make_pair(best_gain,2);
}

void improve_3opt_nn_dlb(vect<int> &curr_path, vect<int> &posns)
{
	int dlb[num] = {false}; // set to off initially
	double current, new_dis;
	int mi, mj, mk, break_flag=0;
	int x1,x2,y1,y2,z1,z2,val;
	pair<double,int> gp;
	bool bet = false;

	update_posns(curr_path, posns);

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			if(dlb[curr_path[i]]==true)			// if dlb set true, dont look at this one
				continue;

			val = curr_path[i];

			for(int dir=0;dir<2;dir++)			// direction - forward(0) and backward(1)
			{
				if(dir==0)
					mi = i;
				else
					mi = (num+i-1)%num;				// considering links between both the neighbours
				
				x1 = curr_path[mi];			
				x2 = curr_path[(mi+1)%num];

				if(dlb[x1]==true || dlb[x2]==true)  // additional speedup with less effect to accuracy
					continue;

				for(int j=0;j<nn;j++)				
				{									
					y1 = nearest_neighbours[x1][j];
					mj = posns[y1];
					y2 = curr_path[(mj+1)%num];
					
					if ((x1 == y1) || (x2 == y1) || (y2 == x1))
          				continue;

					for(int k=0;k<nn;k++)
					{
						z1 = nearest_neighbours[y1][k];
						mk = posns[z1];
						z2 = curr_path[(mk+1)%num];

						if((x1==z1) || (y1==z1))
							continue;

						if(mk>mi)						// is cyclic permutation or not
							bet = (mj>mi) && (mj<mk);
						else if(mk<mi)
							bet = (mj>mi) || (mj<mk);
						else
							bet = false;
						if(!bet)
							continue;

						gp = find_best_case(x1,x2,y1,y2,z1,z2);

						if(gp.first > 0)					// if 3opt move gain found - make move
						{
							dlb[x1]=false; dlb[x2]=false;
							dlb[y1]=false; dlb[y2]=false;
							dlb[z1]=false; dlb[z2]=false;

							make_best_move_pos(curr_path,posns,gp.second,mi,mj,mk);
							break_flag = 1;
							break;		
						}
					}	
					if(break_flag) break;
				}
				if(break_flag)break;
			}

			if(break_flag) 
				break;
			else
				dlb[val] = true;
		}
		if(gp.second==-1)
			break;
	}
	// References : https://tsp-basics.blogspot.com/2017/03/3-opt-basic-with-dlb.html
}

void improve_3opt_dlb(vect<int> &curr_path)		// using don't look bits
{
	int dlb[num] = {false}; // set to off initially
	double current, new_dis;
	int mi, mj, mk, break_flag=0;
	int x1,x2,y1,y2,z1,z2,val;
	pair<double,int> gp;
	bool bet = false;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			if(dlb[curr_path[i]]==true)			// if dlb set true, dont look at this one
				continue;

			val = curr_path[i];
			for(int dir=0;dir<2;dir++)			// direction - forward(0) and backward(1)
			{
				if(dir==0)
					mi = i;
				else
					mi = (num+i-1)%num;				// considering links between both the neighbours
				
				x1 = curr_path[mi];			
				x2 = curr_path[(mi+1)%num];

				if(dlb[x1]==true || dlb[x2]==true)  // additional speedup with less effect to accuracy
					continue;

				for(int j=0;j<num;j++)					
				{									
					y1 = curr_path[j];
					y2 = curr_path[(j+1)%num];
					mj = j;

					if ((x1 == y1) || (x2 == y1) || (y2 == x1))
          				continue;

					for(int k=0;k<num;k++)
					{
						z1 = curr_path[k];
						z2 = curr_path[(k+1)%num];
						mk = k;

						if((x1==z1) || (y1==z1))
							continue;

						if(mk>mi)						// is cyclic permutation or not
							bet = (mj>mi) && (mj<mk);
						else if(mk<mi)
							bet = (mj>mi) || (mj<mk);
						else
							bet = false;
						if(!bet)
							continue;

						gp = find_best_case(x1,x2,y1,y2,z1,z2);

						if(gp.first > 0)					// if 3opt move gain found - make move
						{
							dlb[x1]=false; dlb[x2]=false;
							dlb[y1]=false; dlb[y2]=false;
							dlb[z1]=false; dlb[z2]=false;

							make_best_move(curr_path,gp.second,mi,mj,mk);
							break_flag = 1;
							break;		
						}
					}
					if(break_flag) break;
				}
				if(break_flag) break;
			}
			
			if(break_flag)
				break;
			else
				dlb[val] = true;
		}
		if(gp.second==-1)
			break;
	}
	
	// References : https://tsp-basics.blogspot.com/2017/03/3-opt-basic-with-dlb.html
}

void improve_3opt_nn(vect<int> &curr_path, vect<int> &posns)
{
	double current, new_dis;
	int mi, mj, mk, break_flag=0;
	int x1,x2,y1,y2,z1,z2;
	pair<double,int> gp;
	bool bet = false;

	update_posns(curr_path, posns);

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			mi=i;
			x1 = curr_path[mi];			
			x2 = curr_path[(mi+1)%num];

			for(int j=0;j<nn;j++)				
			{				
				y1 = nearest_neighbours[x1][j];
				mj = posns[y1];
				y2 = curr_path[(mj+1)%num];

				if ((x2 == y1) || (y2 == x1))
      				continue;
				
				for(int k=0;k<nn;k++)
				{
					z1 = nearest_neighbours[y1][k];
					mk = posns[z1];
					z2 = curr_path[(mk+1)%num];
					
					if((x1==z1) || (y1==z1))
						continue;

					if(mk>mi)						// is cyclic permutation or not
						bet = (mj>mi) && (mj<mk);
					else if(mk<mi)
						bet = (mj>mi) || (mj<mk);
					else
						bet = false;
					if(!bet)
						continue;

					gp = find_best_case(x1,x2,y1,y2,z1,z2);

					if(gp.first > 0)					// if 3opt move gain found - make move
					{
						make_best_move_pos(curr_path,posns,gp.second,mi,mj,mk);
						break_flag = 1;
						break;		
					}
				}
				if(break_flag) break;
			}
			if(break_flag) break;
		}
		if(gp.second==-1)
			break;
	}
	
	// References : https://tsp-basics.blogspot.com/2017/03/3-opt-basic-with-dlb.html
}

void improve_3opt_ff(vect<int> &curr_path)
{
	int break_flag = 0;
	int mi, mj, mk;
	int x1,x2,y1,y2,z1,z2;
	pair<double,int> gp;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)					// differing from 2opt here we can have +1 indices - as that will result in 2opt moves
		{
			x1 = curr_path[i];
			x2 = curr_path[(i+1)%num];		

			for(int j=1;j<num-2;j++)			// take careful look at the loop ranges
			{
				y1 = curr_path[(j+i)%num];
				y2 = curr_path[(j+i+1)%num];

				for(int k=j+1;k<num;k++)
				{
					z1 = curr_path[(k+i)%num];
					z2 = curr_path[(k+i+1)%num];

					gp = find_best_case(x1,x2,y1,y2,z1,z2);

					if(gp.first > 0)
					{
						mi = i;
						mj = (j+i)%num;		// take care of these indexes that are passed // dont pass i,j,k instead
						mk = (k+i)%num;
						make_best_move(curr_path,gp.second,mi,mj,mk);

						break_flag = 1;
						break;
					}
				}
				if(break_flag) break;
			}
			if(break_flag) break;
		}
		if(gp.second==-1)
			break;
	}
	// References : https://tsp-basics.blogspot.com/2017/03/3-opt-iterative-basic-algorithm.html
}

void improve_3opt_bf(vect<int> &curr_path)
{
	double best=0,best_case=-1;
	int mi, mj, mk;
	int x1,x2,y1,y2,z1,z2;
	pair<double,int> gp;

	while(true)
	{
		best = 0;
		best_case = -1;

		for(int i=0;i<num;i++)					// differing from 2opt here we can have +1 indices - as that will result in 2opt moves
		{
			x1 = curr_path[i];
			x2 = curr_path[(i+1)%num];		

			for(int j=1;j<num-2;j++)			// take careful look at the loop ranges
			{
				y1 = curr_path[(j+i)%num];
				y2 = curr_path[(j+i+1)%num];

				for(int k=j+1;k<num;k++)
				{
					z1 = curr_path[(k+i)%num];
					z2 = curr_path[(k+i+1)%num];

					gp = find_best_case(x1,x2,y1,y2,z1,z2);

					if(gp.first > best)
					{
						best = gp.first;
						best_case = gp.second;
						mi = i;
						mj = (j+i)%num;		// take care of these indexes that are passed // dont pass i,j,k instead
						mk = (k+i)%num;
					}
				}
			}
		}

		if(best_case==-1)
			break;

		make_best_move(curr_path,best_case,mi,mj,mk);
	}

	/*
	References : 1) https://tsp-basics.blogspot.com/2017/03/3-opt-move.html
			   : 2) https://tsp-basics.blogspot.com/2017/03/3-opt-iterative-general-idea.html
			   : 3) https://tsp-basics.blogspot.com/2017/03/3-opt-iterative-basic-algorithm.html
	*/
}

