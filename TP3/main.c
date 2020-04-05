#include <stdio.h>
#include <stdlib.h>

#define MAX_NODES 100000 //number of nodes
#define MAX_CLUSTERS 100000 //number of clusters


typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

typedef struct {
	unsigned long id;
	unsigned long label;
} node;

typedef struct {
    int exists;
	node *nodes;
	unsigned long degree;
	unsigned long size;
} community;

//edge list structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	edge *edges;//list of edges
	node * nodes;
	unsigned long *d;
	unsigned long *cd;
	unsigned long *adj;
	community *communities;
} graph;

// Generate a graph and store it into file
graph* generate_graph(unsigned long n_nodes, unsigned long n_clusters, double p, double q, char* input){


    unsigned long i,j;
    double temp;
    char* filename = input;

    FILE * fp;
    graph *g=malloc(sizeof(graph));
    g->e = 0;
    g->n = n_nodes;
    g->nodes = malloc(n_nodes * sizeof(node));
    g->edges = malloc(n_nodes * n_nodes * sizeof(edge));



    // creating nodes
    for (i=0; i<n_nodes; i++)
        g->nodes[i].id = i;

    // creating edges
    fp = fopen (filename,"w");
    for (i=0; i<n_nodes; i++){
        for (j=0; j<n_nodes; j++){
            // nodes from the same cluster
            if ((g->nodes[i].id / (n_nodes / n_clusters)) == (g->nodes[j].id / (n_nodes / n_clusters))){
                temp = (double) rand() / (RAND_MAX);
                if (temp<p){
                    g->edges[g->e].s = i;
                    g->edges[g->e].t = j;
                    g->e++;
                    fprintf (fp, "%lu%s%lu \n", i," ", j);
                }
            }
            else{
                temp = (double) rand() / (RAND_MAX);
                if (temp<q){
                    g->edges[g->e].s = i;
                    g->edges[g->e].t = j;
                    g->e++;
                    fprintf (fp, "%lu%s%lu \n", i," ", j);
                }
            }
        }
    }
    fclose(fp);
    return g;
}

void mkadjlist(graph* g){
	unsigned long i,u,v;
	g->d = calloc(g->n,sizeof(unsigned long));

	for (i=0;i<g->e;i++) {
		g->d[g->edges[i].s]++;
        if (g->edges[i].s != g->edges[i].t){
		g->d[g->edges[i].t]++;}
	}
	g->cd=malloc((g->n+1)*sizeof(unsigned long));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+g->d[i-1];
        g->d[i-1]=0;
	}
	g->adj=malloc(2*g->e*sizeof(unsigned long));

	for (i=0;i<g->e;i++) {
		u=g->edges[i].s;
		v=g->edges[i].t;
		g->adj[g->cd[u] + g->d[u]++]=g->edges[i].t;
		if (u!=v){
		g->adj[ g->cd[v] + g->d[v]++ ]=g->edges[i].s;}
	}
	free(g->d);
	free(g->edges);
}

void FisherYates(unsigned long *X, unsigned long n) {
     unsigned long i, j, tmp;

     for (i = n - 1; i > 0; i--) {
         j = rand() % (i + 1);
         tmp = X[j];
         X[j] = X[i];
         X[i] = tmp;
     }
}

void label_propagation(graph *g){
    unsigned long i, no, j, k;
    unsigned long max_id;
    int r, STOP;
    unsigned long iter;
    unsigned long *X;
    unsigned long *freq;
    X = malloc(g->n * sizeof(unsigned long));
    // Step 1 : give unique label to each node randomly from neighbors
    for (i=0; i<g->n; i++){
        r = rand() % (g->cd[i+1]-g->cd[i]);
        g->nodes[i].label = g->adj[g->cd[i]+r];
        X[i] = i; // Initializing the order
    }

    STOP = 0;
    iter = 0;
    while ((iter < 1000) & (STOP == 0))
    {
        iter++;
        //Step 2 : Shuffle the order of the nodes
        FisherYates(X, g->n);
        //Step 3
        STOP = 1;
        for (i=0; i<g->n; i++)
        {
            no = X[i];
            freq = calloc(g->n, sizeof(unsigned long));
            for (j=g->cd[no]; j<g->cd[no+1]; j++)
            {
                // Calculating frequencies of the adj nodes labels
                freq[g->nodes[g->adj[j]].label]++;
            }
            max_id = 0;
            for (k=1; k<g->cd[no+1]-g->cd[no]; k++)
            {
                if (freq[k]> freq[max_id])
                    max_id = k;
            }
            if (g->nodes[no].label != max_id)
            {
                g->nodes[no].label = max_id;
                STOP = 0;
            }
            free(freq);
        }
    }
    printf("%s%lu\n", "Label propagation number of iteration = ",iter);

}

void algo_community(graph *g){
    unsigned long i, no1, no2 , j, k, l, cpt;
    unsigned long max_id;
    int STOP;
    unsigned long iter;
    unsigned long *X;
    X = malloc(g->n * sizeof(unsigned long));
    // Step 1 : Initializing labels for each node
    for (i=0; i<g->n; i++)
    {
        g->nodes[i].label = i;
        X[i] = i; // Initializing the order
    }

    STOP = 0;
    iter = 0;
    while ((iter < 100) & (STOP == 0))
    {
        iter++;
        //Step 2 : Shuffle the order of the nodes
        FisherYates(X, g->n);
        //Step 3
        STOP = 1;
        cpt = 0;
        for (i=0; i<g->n-1; i++)
        {
            no1 = X[i];
            for (j=i+1; j<g->n; j++)
            {
                no2 = X[j];
                // Counter
                for (k=g->cd[no1]; k< g->cd[no1+1]; k++)
                {
                    for (l=g->cd[no2]; l< g->cd[no2+1]; l++)
                    {
                        if (g->adj[l] == g->adj[k])
                            cpt++;
                    }
                }
            if ((cpt >= 15) & (g->nodes[no1].label != g->nodes[no2].label))
            {
                g->nodes[no1].label = g->nodes[no2].label;
                STOP = 0;
            }
            cpt=0;
            }
        }
    }
    printf("%s%lu\n", "proposed algorithm number of iteration = ",iter);

}


void free_graph(graph *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g);
}


int main(){

    time_t t1,t2;
    unsigned long j;
    printf("\n---- Generating graph -----\n");

    t1=time(NULL);
    graph* g;
    g = generate_graph(400, 4, 0.2, 0.001, "label.txt");
    t2=time(NULL);
    printf("\nNumber of nodes: %lu, number of edges: %lu\n", g->n, g->e);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    printf("\n---- Creating adjacency list -----\n");

    t1=time(NULL);
    mkadjlist(g);
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    printf("\n---- Label propagation -----\n");

    t1=time(NULL);
    //label_propagation(g);
    algo_community(g);
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));


	return 0;
}
