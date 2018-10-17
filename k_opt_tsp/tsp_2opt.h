// header file for functions of 2opt local search algorithm and its optimized variations

void improve_2opt_nn_dlb() 	// fixed radius + nearest neighbour
{
	int dlb[num] = {false}; // set to off initially
	double current, new_dis, gain=0;
	int mi, mj, break_flag=0;
	int x1,x2,y1,y2,val;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			if(dlb[min_path[i]]==true)			// if dlb set true, dont look at this one
				continue;

			val = min_path[i];
			for(int dir=0;dir<2;dir++)			// direction - forward(0) and backward(1)
			{
				if(dir==0)
					mi = i;
				else
					mi = (num+i-1)%num;				// considering links between both the neighbours
				
				x1 = min_path[mi];			
				x2 = min_path[(mi+1)%num];
				
				for(int j=0;j<nn;j++)				// no need to specify perfect ranges - handeled in later if conditions	
				{									// needed here for proper functioning of implementation which uses symmetric forms of 2opt move
					y1 = nearest_neighbours[x1][j];
					y2 = min_path[(positions[y1]+1)%num];
					mj = positions[y1];

					if ((x2 == y1) || (y2 == x1))
          				continue;
					
					current = distances[x1][x2] + distances[y1][y2];
					new_dis = distances[x1][y1] + distances[x2][y2];

					gain = current - new_dis;	
					
					if(gain>0)					// first improvement found and performed
					{
						dlb[x1]=false; dlb[x2]=false;
						dlb[y1]=false; dlb[y2]=false;
						reverse_tour_pos(min_path,(mi+1)%num,mj);
						min_dist -= gain;
						break_flag = 1;
						break;		
					}
				}
				if(break_flag) break;
			}

			if(break_flag)
				break;
			else
				dlb[val] = true;
		}
		if(gain<=0)
			break;
	}
	update_positions();
	
	// References : https://tsp-basics.blogspot.com/2017/03/2-opt-with-dlb.html
}


void improve_2opt_dlb()		// using don't look bits
{
	int dlb[num] = {false}; // set to off initially
	double current, new_dis, gain=0;
	int mi, mj, break_flag=0;
	int x1,x2,y1,y2,val;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			if(dlb[min_path[i]]==true)			// if dlb set true, dont look at this one
				continue;

			val = min_path[i];
			for(int dir=0;dir<2;dir++)			// direction - forward(0) and backward(1)
			{
				if(dir==0)
					mi = i;
				else
					mi = (num+i-1)%num;				// considering links between both the neighbours
				
				x1 = min_path[mi];			
				x2 = min_path[(mi+1)%num];
				
				for(int j=0;j<num;j++)				// no need to specify perfect ranges - handeled in later if conditions	
				{									// needed here for proper functioning of implementation which uses symmetric forms of 2opt move
					y1 = min_path[j];
					y2 = min_path[(j+1)%num];
					mj = j;

					if ((x1 == y1) || (x2 == y1) || (y2 == x1))
          				continue;
					
					current = distances[x1][x2] + distances[y1][y2];
					new_dis = distances[x1][y1] + distances[x2][y2];

					gain = current - new_dis;	
					
					if(gain>0)					// first improvement found and performed
					{
						dlb[x1]=false; dlb[x2]=false;
						dlb[y1]=false; dlb[y2]=false;

						reverse_tour(min_path,(mi+1)%num,mj);
						min_dist -= gain;
						break_flag = 1;
						break;		
					}
				}
				if(break_flag) break;
			}

			if(break_flag)
				break;
			else
				dlb[val] = true;
		}
		if(gain<=0)
			break;
	}
	update_positions();
	
	// References : https://tsp-basics.blogspot.com/2017/03/2-opt-with-dlb.html
}


void improve_2opt_nn()		// using nearest neighbours lists
{
	nn = min(nn, num-1);
	double current, new_dis,gain=0;
	int mi, mj, break_flag=0,range;
	int x1,x2,y1,y2;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			x1 = min_path[i];
			x2 = min_path[(i+1)%num];			

			for(int j=0;j<nn;j++)	
			{
				y1 = nearest_neighbours[x1][j];
				y2 = min_path[(positions[y1]+1)%num];

				if ((x2 == y1) || (y2 == x1))	// edges selected should not be adjacent
          			continue;

				current = distances[x1][x2] + distances[y1][y2];
				new_dis = distances[x1][y1] + distances[x2][y2];

				gain = current - new_dis;	// we have to maximise the distance
				
				if(gain>0)					// first improvement found and performed
				{
					mi = i; mj = positions[y1];
					reverse_tour_pos(min_path,(mi+1)%num,mj);	// modified to adjust positions vector also
					min_dist -= gain;
					break_flag = 1;
					break;		
				}
			}
			if(break_flag) break;
		}
		if(gain<=0)
			break;
	}
	
	// References : https://tsp-basics.blogspot.com/2017/03/2-opt-with-neighbor-lists.html
}

void improve_2opt_fr()		// using fixed radius search idea 
{
	double current, new_dis, gain=0, radius;
	int mi, mj, break_flag=0;
	int x1,x2,y1,y2;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num;i++)
		{
			for(int dir=0;dir<2;dir++)			// direction - forward(0) and backward(1)
			{
				if(dir==0)
				{
					mi = i;
					x1 = min_path[mi];
					x2 = min_path[(mi+1)%num];
				}
				else if(dir==1)
				{
					mi = (num+i-1)%num;
					x1 = min_path[(mi+1)%num];
					x2 = min_path[mi];
				}
				radius = distances[x1][x2];

				for(int j=0;j<num;j++)				// no need to specify perfect ranges - handeled in later if conditions	
				{									// needed here for proper functioning of implementation which uses symmetric forms of 2opt move
					y1 = min_path[j];
					y2 = min_path[(j+1)%num];
					mj = j;

					if(dir==1)
						swap(y1,y2);

					if ((x1 == y1) || (x2 == y1) || (y2 == x1))
          				continue;
					
					if(distances[x2][y2]>radius)	// fixed radius search check - checks forward and backward(later - using symmetric nature of loop)
						continue;					// optimizes by avoiding below steps of checking

					current = distances[x1][x2] + distances[y1][y2];
					new_dis = distances[x1][y1] + distances[x2][y2];

					gain = current - new_dis;	
					
					if(gain>0)					// first improvement found and performed
					{
						reverse_tour(min_path,(mi+1)%num,mj);
						min_dist -= gain;
						break_flag = 1;
						break;		
					}
				}
				if(break_flag) break;
			}
			if(break_flag) break;
		}
		if(gain<=0)
			break;
	}
	update_positions();
	
	// References : https://tsp-basics.blogspot.com/2017/03/2-opt-basic-with-fixed-radius-search.html
}

void improve_2opt_ff()		// first found improvement executed
{
	double current, new_dis,gain=0;
	int mi, mj, break_flag=0,range;
	int x1,x2,y1,y2;

	while(true)
	{
		break_flag = 0;

		for(int i=0;i<num-2;i++)
		{
			x1 = min_path[i];
			x2 = min_path[i+1];			// a small optimisation // accessed one if put out of j loop

			if(i==0)
				range = num-1;			// just a optimisation step, avoided itself as current = new_dis in case of overlap as x1=y2
			else 
				range = num;

			for(int j=i+2;j<range;j++)	// edges selected should not be adjacent
			{
				y1 = min_path[j];
				y2 = min_path[(j+1)%num];

				current = distances[x1][x2] + distances[y1][y2];
				new_dis = distances[x1][y1] + distances[x2][y2];

				gain = current - new_dis;	// we have to maximise the distance
				
				if(gain>0)					// first improvement found and performed
				{
					mi = i; mj = j;
					reverse_tour(min_path,mi+1,mj);
					min_dist -= gain;
					break_flag = 1;
					break;		
				}
			}
			if(break_flag) break;
		}
		if(gain<=0)
			break;
	}
	update_positions();
	
	// References : https://tsp-basics.blogspot.com/2017/03/2-opt-basic-iterative-algorithm.html
}

void improve_2opt_bf()		// best found improvement executed
{
	double current, new_dis,gain=0,best=0;
	int mi, mj, range;
	int x1,x2,y1,y2;

	while(true)
	{
		best = 0;

		for(int i=0;i<num-2;i++)
		{
			x1 = min_path[i];
			x2 = min_path[i+1];			// a small optimisation // accessed one if put out of j loop

			if(i==0)
				range = num-1;			// just a optimisation step, avoided itself as current = new_dis in case of overlap as x1=y2
			else 
				range = num;

			for(int j=i+2;j<range;j++)	// edges selected should not be adjacent
			{
				y1 = min_path[j];
				y2 = min_path[(j+1)%num];

				current = distances[x1][x2] + distances[y1][y2];
				new_dis = distances[x1][y1] + distances[x2][y2];

				gain = current - new_dis;	// we have to maximise the distance
				if(gain>best)
				{
					best = gain;
					mi = i; mj = j;			// best improvement stored and used at last 
				}							// there is another variant also - first improving 2-opt move (takes more time generally)
			}
		}
		if(best<=0)
			break;

		reverse_tour(min_path,mi+1,mj);
		min_dist -= best;

	}
	update_positions();

	/*
	References : 1) http://on-demand.gputechconf.com/gtc/2014/presentations/S4534-high-speed-2-opt-tsp-solver.pdf
			   : 2) https://tsp-basics.blogspot.com/2017/03/2-opt-basic-iterative-algorithm.html
			   : 3) https://tsp-basics.blogspot.com/2017/03/2-opt-move.html
	*/
}
