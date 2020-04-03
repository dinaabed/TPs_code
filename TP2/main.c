#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>//to estimate the runing time
#include <math.h>
#define NLINKS 100000000

// Structure d'un arret
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

//compute the maximum of three unsigned long
unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

unsigned long graph_size(char* name_file)
{
    FILE* file = fopen (name_file, "r");
    unsigned long i = 0;
    unsigned long j =0;
    unsigned long l = 0;

    while (fscanf(file,"%lu %lu", &i, &j)==2)
    {
        l=max3(l,i,j);
        //e++;
    }
    fclose (file);
    l++;

    return l;


}
int* degree_out (char* name_file)
{
    unsigned long l = graph_size(name_file);
    int* nodes = malloc(l*sizeof(int));
    int i=0;
    // intialisation à zéro
    for (i=0; i<l ; i++)
    {
        nodes[i]=0;
    }
    FILE* file = fopen (name_file, "r");

    unsigned long j = 0;
    unsigned long k =0 ;
    while (fscanf(file,"%lu %lu", &j, &k)==2)
    {
        nodes[j] = nodes[j] + 1;
    }
    fclose (file);
    return nodes ;
}

int* degree_in(char* name_file)
{
    unsigned long l = graph_size(name_file);
    int* nodes = malloc(l*sizeof(int));
     // intialisation à zéro
    int i=0;
    for (i=0; i<l ; i++)
    {
        nodes[i]=0;
    }
    FILE* file = fopen (name_file, "r");
    unsigned long j = 0;
    unsigned long k =0 ;
    while (fscanf(file,"%lu %lu", &j, &k)==2)
    {
        nodes[k] = nodes[k] + 1;
    }
    fclose (file);
    return nodes ;
}
double sum(double* v, int n )
{
    double s =0;
    int i =0;
    for (i=0; i<n ; i++)
        s += fabs(v[i]) ;
    return s;
}
double* PageRank(char* input, double alpha, unsigned int t)
{
    int* degrees_out = degree_out(input);
    unsigned long n = graph_size(input);
    double* P0 = malloc(n*sizeof(double));
    double* P = malloc(n*sizeof(double));
    int k=0;
    int l=0;
    int b=0;
    unsigned long i = 0;
    unsigned long j =0;
    // INITIALIZATION
    for (k=0 ; k<n ; k++)
    {
        P[k] = 0;
        P0[k] = 1/((double)n);
    }

    k=0;
    printf("start_bis alpha= %lf\n",alpha);
    int stop=0;
    float error = 0.0;
    while (!stop && (k<t))
    {
        for (b=0 ; b<n ; b++)
        {
            P[b] = 0;
        }
        FILE* file = fopen (input, "r");

        while (fscanf(file,"%lu %lu", &i, &j)==2)
        {
            if (degrees_out[i]!=0)
            {
                P[j] = P[j] + P0[i]/((double)degrees_out[i]);
            }
            else
            {
                P[j] = P[j] +  P0[i]/((double)n);
            }
        }

        fclose (file);

        for (l=0 ; l<n ; l++)
        {
            P[l] = (1-alpha)*P[l] + alpha/((double)n);
        }
        double S=sum(P,n);
        l=0;
        for (l=0 ; l<n ; l++)
        {
            P[l] = P[l] + (1-S)/((double)n);

        }
        error = 0.0;
        i=0;
        for(i=0; i<n; i++)
        {
            error =  error + fabs(P[i] - P0[i]);
        }
        if (error < 0.000001)
        {
            stop = 1;
            printf("number of iterations: %d\n ", k);
        }
        l=0;
        for (l=0 ; l<n ; l++)
        {
            P0[l] = P[l] ;
        }
        k++;
    }

    return P;
}

int* count_degree (char* input)
{
    unsigned long n = graph_size(input);
    int* nodes = malloc(n*sizeof(int));
    int k=0;
    for (k=0; k<n ; k++)
    {
        nodes[k]=0;
    }
    FILE* file = fopen (input, "r");
    unsigned long i = 0;
    unsigned long j =0 ;
    while (fscanf(file,"%lu %lu", &i, &j)==2)
    {
        nodes[i] = nodes[i] + 1;
        nodes[j] = nodes[j] +1;
    }
    fclose (file);
    return nodes ;
}

int argmin( int* table,int n )
{
    int argmin_table = 0;
    int j=0;
    for( j=1; j<n; j++)
    {
        if (table[j]< table[argmin_table])
        {
            argmin_table = j ;
        }
    }

    return  argmin_table;
}


int* k_core_decomposition(char* input )
{
    int* degrees = count_degree(input);
    unsigned long n = graph_size(input);
    unsigned long i =n ;
    int c =0;

    int* C = malloc( n* sizeof(int));
    int argmin_degrees =0;
    unsigned long k =0;
    unsigned long l=0;
    time_t t1,t2;
    t1=time(NULL);
    while(i!=0)
    {
        argmin_degrees = argmin( degrees, n );
        if (c < degrees[argmin_degrees])
        {
            c = degrees[argmin_degrees];
        }
        degrees[argmin_degrees] = 2*n ; // on essaye d'éliminer l'élément ) argmin_degrees


        C[argmin_degrees] = c;
        FILE* file = fopen (input, "r");
        while (fscanf(file,"%lu %lu", &l, &k)==2)
        {
            if( l == argmin_degrees )
            {
                degrees[k]-=1 ;
            }
            else if ( k == argmin_degrees)
            {
                degrees[l]-=1 ;
            }

        }
        fclose (file);
        i-=1;
    }
    t2=time(NULL);
    printf(" K_core decomposition running time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));



    return C ;
}



int main(int argc,char** argv)
{


    printf("Reading edgelist from file %s\n",argv[1]);

    unsigned long n = graph_size(argv[1]);

    printf("size of graph: %lu \n",n);

    bool PageRANK;
    PageRANK = false;

    if (PageRANK)
    {


        double* rankes15 = PageRank(argv[1], 0.15, 1000);
        double* rankes01 = PageRank(argv[1], 0.1, 1000);
        double* rankes02 = PageRank(argv[1], 0.2, 1000);
        double* rankes05 = PageRank(argv[1], 0.5, 1000);
        double* rankes09 = PageRank(argv[1], 0.9, 1000);
        printf("PageRank runned");
        int* degrees_in = degree_in(argv[1]);
        int* degrees_out = degree_out(argv[1]);

        FILE * fp;

        fp = fopen ("C:\\Users\\abedd\\Desktop\\code\\output_PageRank.txt","w");



        int k=0;
        printf("writing PageRank \n");
        for (k=0; k<n ; k++)
        {
            fprintf (fp, "%f,%f,%f,%f,%f,%d,%d \n ",rankes15[k],rankes09[k],rankes05[k],rankes02[k],rankes01[k],degrees_in[k],degrees_out[k]);
        }
        fclose (fp);
        }
    else
    {

        FILE * fk;


        int* k_core = k_core_decomposition(argv[1]);
        int* degrees = count_degree(argv[1]);
        printf("writing k_Core \n");
        fk = fopen ("C:\\Users\\abedd\\Desktop\\code\\output_Kcore.txt","w");
        int k=0;
        for (k=0; k<n ; k++)
        {
            fprintf (fk, "%d,%d \n ",degrees[k],k_core[k]);
        }
        fclose (fk);
        }








    return 0;
}
