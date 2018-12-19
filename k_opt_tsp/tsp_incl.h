#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <queue>
#include <omp.h>
#include <unordered_set>

#define vect vector
#define add push_back

using namespace std;

int num,nn;

vect<pair<double,double> > points;
vect<vect<int> > nearest_neighbours;			// nearest neighbours lists
vect<vect<double> > distances;

double min_dist;
vect<int> idx;									// for brute force soln
vect<int> min_path;
vect<int> positions;

void update_posns(vect<int> &curr_path, vect<int> &posns)
{
	for(int i=0;i<num;i++)
		posns[curr_path[i]]=i;
}

void update_dist(vect<int> &curr_path, double &dist)
{
	dist = 0;

	for(int i=1;i<num;i++)										
			dist += distances[curr_path[i-1]][curr_path[i]];

	dist += distances[curr_path[num-1]][curr_path[0]];
}

void reverse_tour_pos(vect<int> &curr_path, vect<int> &posns, int i, int j)		// general implementation - takes care even is i>j
{
	int invs = ((num+j-i+1) % num)/2;
	while(invs--)				// reversing subarray tour
	{
		swap(curr_path[i],curr_path[j]);
		posns[curr_path[i]]=i;		// also changes positions as needed at runtime
		posns[curr_path[j]]=j;

		i = (i+1)%num;		// cause tours are cyclic
		j = (num+j-1)%num;
	}
}

void reverse_tour(vect<int> &curr_path, int i, int j)		// general implementation - takes care even is i>j
{
	int invs = ((num+j-i+1) % num)/2;

	while(invs--)				// reversing subarray tour
	{
		swap(curr_path[i],curr_path[j]);

		i = (i+1)%num;		// cause tours are cyclic
		j = (num+j-1)%num;
	}
}


void construct_nn(int range)	// initialise path using nearest unvisited neighbour approach (all/one vertex chosen at starting)
{
	vect<bool> visited(num);
	vect<int> temp_path(num,0);
	double temp_dist;

	min_dist = INT_MAX;

	int near,ix, curr=0;
	double near_dis = 0;

	for(int j=0;j<range;j++)	// time complexity of this best nn construction in O(n^3) if range = num, or O(n^2) is range = 1
	{
		for(int i=0;i<num;i++)
			visited[i] = false;

		ix=0;
		curr = j;
		temp_path[ix++] = curr;
		visited[curr] = true;

		while(true)			
		{
			near=-1;
			near_dis = INT_MAX;

			for(int i=0;i<num;i++)
			{
				if(i!=curr && (!visited[i]) && (distances[curr][i] < near_dis) )
				{
					near_dis = distances[curr][i];
					near = i;
				}
			}

			if(near==-1)
				break;
			curr = near;
			visited[curr] = true;
			temp_path[ix++] = curr;
		}

		temp_dist = 0;	

		for(int i=1;i<num;i++)										// can also calculate this during above loop
			temp_dist += distances[temp_path[i-1]][temp_path[i]];

		temp_dist += distances[temp_path[num-1]][temp_path[0]];

		if(temp_dist<min_dist)
		{
			min_dist = temp_dist;
			min_path = temp_path;
		}
	}	
	update_posns(min_path, positions);
}

void construct_nn_opt(int range)	// initialise path using nearest unvisited neighbour approach (all/one vertex chosen at starting)
{									// optimised nearest neighbour construction method - uses preprocessing of sorting earlier // fast
	vect<bool> visited(num);
	vect<int> temp_path(num,0);
	double temp_dist;

	min_dist = INT_MAX;

	int near,ix,num_vis=0, curr=0;
	double near_dis = 0;

	for(int j=0;j<range;j++)	// time complexity of this best nn construction in O(n^3) if range = num, or O(n^2) is range = 1
	{
		for(int i=0;i<num;i++)
			visited[i] = false;

		num_vis = 0; ix=0;

		curr = j;
		temp_path[ix++] = curr;
		visited[curr] = true;

		while(num_vis!=num)			
		{
			for(int i=0;i<num-1;i++)
			{
				if(!visited[nearest_neighbours[curr][i]])
				{
					near_dis = distances[curr][nearest_neighbours[curr][i]];
					near = nearest_neighbours[curr][i];
					break;
				}
			}

			curr = near;
			visited[curr] = true;
			num_vis++;
			temp_path[ix++] = curr;
		}

		temp_dist = 0;	

		for(int i=1;i<num;i++)										// can also calculate this during above loop
			temp_dist += distances[temp_path[i-1]][temp_path[i]];

		temp_dist += distances[temp_path[num-1]][temp_path[0]];

		if(temp_dist<min_dist)
		{
			min_dist = temp_dist;
			min_path = temp_path;
		}
	}	
	update_posns(min_path, positions);
}

void construct_random()		// generates random walk initial path
{							// in the order of given nodes
	min_dist = 0;
	
	for(int i=1;i<num;i++)
	{
		min_dist += distances[i-1][i];
		min_path.add(i-1);
	}
	min_dist += distances[num-1][0];
	min_path.add(num-1);

	update_posns(min_path, positions);
}


bool compare_nn(pair<int,double> &a, pair<int,double> &b)
{
	if(a.second<b.second)
		return true;
	return false;
}

void prepare_neighbours()
{
	nearest_neighbours = vect<vect<int> > (num,vect<int>(num-1,0));	// num-1 to avoid storing the same point's index

	vect<pair<int,double> > temp(num);

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			if(i==j)
				temp[j] = make_pair(j, INT_MAX);			// had to do this to handle some cases having same points
			else
				temp[j] = make_pair(j, distances[i][j]);
		}

		sort(temp.begin(), temp.end(), compare_nn);

		for(int j=0;j<num-1;j++)
			nearest_neighbours[i][j] = temp[j].first;
	}
}


void print_distances()
{
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
			cout<<distances[i][j]<<" ";
		cout<<endl;
	}
}

double eucl_dist(pair<double,double> &a, pair<double,double> &b, int type)
{
	double x_dist = (a.first-b.first)*(a.first-b.first);
	double y_dist = (a.second-b.second)*(a.second-b.second);
	
	if(type==0)
		return sqrt(x_dist + y_dist);

	else if(type==1)								// for att type format
		return ceil(sqrt((x_dist + y_dist)/10.0));
}

void calc_distances(int type)
{
	distances = vect<vect<double> > (num, vect<double> (num,0));	// 2d array size allocation

	for(int i=0;i<num;i++)
	{
		for(int j=i+1;j<num;j++)
		{
			distances[i][j] = eucl_dist(points[i],points[j],type);
			distances[j][i] = distances[i][j];						// symmetry property for optimization
		}				
	}

	positions = vect<int> (num);
}


void store_points(char dataset_file[])
{
	freopen(dataset_file, "r", stdin);
	cin>>num;

	pair<double,double> temp;		

	for(int i=0;i<num;i++)
	{
		cin >> temp.first >> temp.second;
		points.add(temp);
		idx.add(i);
	}
}

void follow_dist(char dataset_file[])
{
	freopen(dataset_file, "r", stdin);
	cin>>num;

	distances = vect<vect<double> > (num, vect<double> (num,0));

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
			cin>>distances[i][j];

		idx.add(i);
	}
	positions = vect<int> (num);
}