#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>//to estimate the runing time

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed
#define NNODES 100000000 //maximum number of edges for memory allocation, will increase if needed
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAXL 1000
#include <windows.h> // notice this! you need it! (windows)

// ********* graph structures and loading functions ***********
//edge structure:
typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

//edge list structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	edge *edges;//list of edges
} edgelist;

//adjmatrix structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	bool *mat;//adjacency matrix
	unsigned long *nodes;
} adjmatrix;

//adjarray structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	unsigned long *cd;//cumulative degree cd[0]=0 length=n+1
	unsigned long *adj;//concatenated lists of neighbors of all nodes
	unsigned long *nodes;

} adjarray;

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
edgelist* readedgelist(char* input){
	unsigned long e1=NLINKS;
	unsigned long e2=NNODES;
	FILE *file=fopen(input,"r");
    char line[MAXL];
	edgelist *g=malloc(sizeof(edgelist));
	unsigned long * nodes;
	nodes = malloc(e2 * sizeof(unsigned long));
	g->n=0;
	g->e=0;
	g->edges=malloc(e1*sizeof(edge));//allocate some RAM to store edges

	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2){
            if (MAX(g->edges[g->e].s, g->edges[g->e].t) > e2){
                    e2+=NLINKS;
                    nodes=(unsigned long *)realloc(nodes, e2 * sizeof(unsigned long));
            }
            nodes[g->edges[g->e].s] = 1;
            nodes[g->edges[g->e].t] = 1;

            if (g->e++==e1) {//increase allocated RAM if needed
                e1+=NLINKS;
                g->edges=realloc(g->edges,e1*sizeof(edge));

            }
        }
	}
	fclose(file);
	int i;
    for(i=0;i<e2;i++)
	{
		//g->n+=1;
		g->n+=*(nodes+i);
	}

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

// free edgelist from memory
void free_edgelist(edgelist *g){
	free(g->edges);
	free(g);
}

//reading the edgelist from file and store as adjmatrix
adjmatrix* readgraphtomatrix(char* input){
	unsigned long e1=NLINKS;
	unsigned long e2=NNODES;
	unsigned long k,i,j;
	char line[MAXL];

	adjmatrix *g=malloc(sizeof(adjmatrix));

	g->n=0;
	g->e=0;
	g->nodes=malloc(e2 * sizeof(unsigned long));

    FILE *file=fopen(input,"r");
	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(i), &(j))==2){
            if (MAX(i, j) > e2){
                    e2+=NLINKS;
                    g->nodes=(unsigned long*)realloc(g->nodes, e2 * sizeof(unsigned long));
            }
            g->nodes[i] = 1;
            g->nodes[j] = 1;

            if (g->e++==e1) {//increase allocated RAM if needed
                e1+=NLINKS;

            }
        }
	}
	fclose(file);

    for(i=0;i<e2;i++)
	{
		g->n+=*(g->nodes+i);
	}

	i=0;
	k=0;
    while(i<e2 & k<g->n)
    {
        if(g->nodes[i] == 1)//if the node exists
        {
            g->nodes[i] = k;
            k++;
        }
        else{
            g->nodes[i] = -1;
        }
        i++;
    }

	g->mat=calloc(g->n*g->n,sizeof(bool));
	file=fopen(input,"r");
	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(i), &(j))==2){
            g->mat[g->nodes[j]+g->n*g->nodes[i]] = 1;
            g->mat[g->nodes[i]+g->n*g->nodes[j]] = 1;
        }
	}
	fclose(file);

	return g;
}

// free adjmatrix from memory
void free_adjmatrix(adjmatrix *g){
	free(g->nodes);
	free(g->mat);
	free(g);
}

//reading the edgelist from file and store as adjarray
adjarray* readgraphtoarray(char* input){
	unsigned long e1=NLINKS;
	unsigned long e2=NNODES;
	unsigned long k,i,j,u,v;;
	char line[MAXL];

	adjarray *g=malloc(sizeof(adjarray));

	g->n=0;
	g->e=0;
	g->nodes=malloc(e2 * sizeof(unsigned long));

    FILE *file=fopen(input,"r");
	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(i), &(j))==2){
            if (MAX(i, j) > e2){
                    e2+=NLINKS;
                    g->nodes=(unsigned long *)realloc(g->nodes, e2 * sizeof(unsigned long));
            }
            g->nodes[i] = 1;
            g->nodes[j] = 1;

            if (g->e++==e1) {//increase allocated RAM if needed
                e1+=NLINKS;

            }
        }
	}
	fclose(file);

    for(i=0;i<e2;i++)
	{
		g->n+=*(g->nodes+i);
	}

	i=0;
	k=0;
    while(i<e2 & k<g->n)
    {
        if(g->nodes[i] == 1)//if the node exists
        {
            g->nodes[i] = k;
            k++;
        }
        else{
            g->nodes[i] = -1;
        }
        i++;
    }

	unsigned long *d=calloc(g->n,sizeof(unsigned long));
    g->cd=malloc((g->n+1)*sizeof(unsigned long));
	g->cd[0]=0;
	g->adj=malloc(2*g->e*sizeof(unsigned long));

    file=fopen(input,"r");
	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(i), &(j))==2){
            d[g->nodes[i]]++;
            d[g->nodes[j]]++;
        }
	}
	fclose(file);

	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		d[i-1]=0;
	}

    file=fopen(input,"r");
	while (fgets(line, sizeof line, file)){
        if (*line == '#')
            continue;
        if (sscanf(line,"%lu %lu", &(i), &(j))==2){
            g->adj[g->cd[g->nodes[i]] + d[g->nodes[i]]++] = g->nodes[j];
            g->adj[g->cd[g->nodes[j]] + d[g->nodes[j]]++] = g->nodes[i];
        }
	}
	fclose(file);
	free(d);

	return g;
}

//free adjarray from memory
void free_adjarray(adjarray *g){
	free(g->nodes);
	free(g->cd);
	free(g->adj);
	free(g);
}

// ************ Queue structures and functions************

// queue element structure
typedef struct queue_element queue_element;
typedef struct queue_element{
    unsigned long node;
    queue_element *next;
};

// queue structure
typedef struct {
    queue_element *first;
} queue;

// initialize queue
queue *initialization (){
    queue *q=malloc(sizeof(queue));
    queue_element *element=malloc(sizeof(queue_element));
    if (q==NULL || element==NULL)
        exit(0);
    element=NULL;
    q->first=element;
    return(q);
}

// add element to queue
void add(queue *q, unsigned long v){
    queue_element *element=malloc(sizeof(queue_element));
    if (q==NULL || element==NULL)
        exit(0);
    element->node=v;
    element->next=NULL;
    queue_element* last=q->first;
    if(last==NULL)
        q->first=element;
    else {
        while (last->next!=NULL)
            last=last->next;
        last->next=element;
    }
}

// pop element from queue
void pop(queue *q){
    if (q==NULL)
        exit(0);
    q->first=q->first->next;
}

// check if queue if empty
int isempty(queue *q){
    if (q==NULL || q->first==NULL)
        return(1);
    else
        return(0);
}


// ***** Connected  components ********
bool* BFS(adjarray* g, unsigned long s, bool* marked){
    unsigned long n = g->n;
    unsigned long u=0;
    unsigned long size=0;
    unsigned long i=0;

    queue  *q = NULL;
    q = initialization();
    add(q,s);
    marked[s]=1;
    while (isempty(q)==0){
        size++;
        u=q->first->node;
        pop(q);
        for (i=g->cd[u];i<g->cd[u+1];i++){
            if (marked[g->adj[i]]==0){
                add(q, g->adj[i]);
                marked[g->adj[i]]=1;
            }
        }
    }
    return size;
}

void graphconnectedcomponents(adjarray* g){
    bool* marked=malloc(g->n*sizeof(bool));
    unsigned long i=0;
    unsigned long size, count;
    count = 0;
    for (i=0;i<g->n;i++)
        marked[i]=0;

    for(i = 0;i < g->n;++i)
    {
        if(marked[i] == 0){
            size = BFS(g,i, marked);
            count++;
            float percent = (float) size/g->n;
            printf("\nConnected component %lu%s%lu%s%f\n",count,", size = ", size, ", fraction = ", percent);
        }
    }
    free(marked);
}


// ****** diameter ******

int* diameter_bfs(adjarray *g, unsigned long s){
	unsigned long compt=2;
	unsigned long i,j,v;

	int *diameter;
	unsigned long *list;
	diameter=malloc((g->n+1)*sizeof(int));
	list=malloc((g->n+1)*sizeof(unsigned long));

	for (i=0;i<=g->n;i++){
		diameter[i]=-1;
    }
	list[1]=s;
	diameter[s]=0;

	for (i=1;i<compt;i++) {
		v=list[i];
		for (j=g->cd[v];j<g->cd[v+1];j++) {
			if (diameter[g->adj[j]]==-1) {
				list[compt++]=g->adj[j];
				diameter[g->adj[j]]=diameter[v]+1;
            }
		}
	}
    free(list);
	return diameter;
}

unsigned long max_of_vector(unsigned long n, int* vect){
	unsigned long i,imax;
	unsigned vmax=0;
	for (i=1;i<=n;i++){
		if ((vect[i]!=-1) && (vect[i]>=vmax)){
			imax=i;
			vmax=vect[i];
		}
	}
	return imax;
}

int lowerbound(adjarray* g){
    int *diameter;
    unsigned long max_diam=0;
    unsigned long pos=1;
    unsigned long i;
    int *nodes;
    time_t t1;

    nodes=calloc((g->n+1),sizeof(int));
    t1=time(NULL);

    for (i=1;i<g->n;i++) {
        if (nodes[pos]==0)
            i--;
        else{
            pos=i;
            if (nodes[pos]==1)
                continue;
        }
        nodes[pos]=1;
        diameter=diameter_bfs(g,pos);
        pos=max_of_vector(g->n,diameter);
        max_diam=MAX(max_diam,diameter[pos]);
        free(diameter);
        if ((time(NULL)-t1) > 6*60){
            break;
        }
    }

    return max_diam;
}

// ****** Triangles ******

int countTriangle(adjarray* g, char* input){
    int count_Triangle = 0;

    int k=0;
    int l=0;

    FILE* file = fopen (input, "r");
    unsigned long i = 0;
    unsigned long j =0;

    while (fscanf(file,"%lu %lu", &(i), &(j))==2)
    {
        for (k=g->cd[g->nodes[i]]; k<g->cd[g->nodes[i]+1]; k++){
            for (l=g->cd[g->nodes[j]]; l<g->cd[g->nodes[j]+1];l++){
                if ( (g->adj[k]==g->adj[l]) && (g->nodes[i]<g->nodes[j]) && (g->nodes[j]<g->adj[k])) {count_Triangle++;}
            }
        }

    }
    fclose (file);
    return count_Triangle;

}


int main(int argc,char** argv){

    const char *filenames[4];
    filenames[0] = "D:/Projects/RO/TP1/email-Eu-core.txt";
    filenames[1] = "D:/Projects/RO/TP1/com-amazon.ungraph.txt";
    filenames[2] = "D:/Projects/RO/TP1/com-lj.ungraph.txt";
    filenames[3] = "D:/Projects/RO/TP1/com-orkut.ungraph.txt";
    filenames[4] = "D:/Projects/RO/TP1/com-friendster.ungraph.txt";
    filenames[5] = "D:/Projects/RO/TP1/test.txt";

    time_t t1,t2;
    unsigned long i,j;
    int idx;
    idx=2;

    printf("File: %s\n", filenames[idx]);

    printf("\n----- Reading graph from file -----\n");

    edgelist* g;
    t1=time(NULL);
    g = readedgelist(filenames[idx]);
    printf("Number of nodes in graph: %lu\n", g->n);
    printf("Number of edges in graph: %lu\n", g->e);
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
    free_edgelist(g);


    printf("\n----- Reading graph and storing it in matrix-----\n");

    adjmatrix* gmat;
    t1=time(NULL);
    gmat = readgraphtomatrix(filenames[idx]);
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
    free_adjmatrix(gmat);

    printf("\n----- Reading graph and storing it in array-----\n");

    adjarray* garray;
    t1=time(NULL);
    garray = readgraphtoarray(filenames[idx]);
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    printf("\n---- Connected componenets -----\n");

    t1=time(NULL);
    graphconnectedcomponents(garray);
    t2=time(NULL);

    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    printf("\n---- Triangles -----\n");

    t1=time(NULL);
    printf("Number of triangles: %d\n", countTriangle(garray, filenames[idx]));
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    printf("\n---- Lower bound -----\n");

    t1=time(NULL);
    printf("Lower bound: %d\n", lowerbound(garray));
    t2=time(NULL);
    printf("Elapsed time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

    return 0;
}
