// header file for mutation operators of genetic crossover approach for tsp
// swap_2 and swap_3 will be selected with prob = 15% each, whereas seg_rev and rand_shuff will be selected with prob = 40 and 30% each (as they gve better results)

void point_swap_2(vect<int> &path)
{
	int a = rand()%num;
	int b = a;

	while(b==a)
		b = rand()%num;

	swap(path[a],path[b]);
}

void point_swap_3(vect<int> &path)
{
	int a = rand()%num;
	int b = a;

	while(b==a)
		b = rand()%num;

	int c = b;
	while(c==b || c==a)
		c = rand()%num;

	swap(path[a],path[b]);		//a->b, b->a
	swap(path[a],path[c]);		//b->c, c->a
}


void seg_rev(vect<int> &path)	// reverse a random length segment	// reverse sequence mutation (RSM)
{
	int a = rand()%num;
	int b = a;

	while(b==a)
		b = rand()%num;

	if(a>b)
		swap(a,b);

	reverse(path.begin()+a, path.begin()+b);
}

void part_shuff(vect<int> &path)	// shuffle a random length segment	// Partial Shuffle mutation (PSM)
{
	int a = rand()%num;
	int b = a;

	while(b==a)
		b = rand()%num;

	if(a>b)
		swap(a,b);

	random_shuffle(path.begin()+a, path.begin()+b, myrandom);
}

void uhx_mutation(vect<int> &path)		// works very good if every offspring goes through only this, no probabilistic selection
{										// can add probabilistic for uhx crossover
	int a = rand()%num;
	int ne = nearest_neighbours[a][0];
	int b;

	for(int i=0;i<num;i++)
	{
		if(path[i]==ne)
		{
			b=i;
			break;
		}
	}

	if(b<a)
		swap(a,b);

	reverse(path.begin()+a, path.begin()+b);
}

void mutate()
{
	for(int i=0;i<popu;i++)
	{
		double todo = (double) rand()/RAND_MAX;
		
		// if(todo < 0.2)
		// {
		// 	double val = (double) rand()/RAND_MAX;

		// 	if(val < 0.15)
		// 		point_swap_2(population[i].path);
		// 	else if(val < 0.30)
		// 		point_swap_3(population[i].path);
		// 	if(val < 0.70)
		// 		seg_rev(population[i].path);
		// 	else
		// 		part_shuff(population[i].path);
		// }
		
		uhx_mutation(population[i].path);		// use with vgx, uhx only // comment above block // but with pmx,ox use above block and comment this line
	}
}
