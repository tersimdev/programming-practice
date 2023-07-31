#include<stdio.h>
#include <stdlib.h>

// Author: Terence
// This code uses HLD to do multiple queries on a tree, 
// where each query is a list of nodes (with <= value) from node A to node B.
// Has many comments documenting each step, breaking down the LCA, Segment Tree & HLD.





int* city_population (int N, int* population, int** road, int Q, int** cities) ;

int main() {
    int N;
    scanf("%d", &N);
    int i_population;
    int *population = (int *)malloc(sizeof(int)*(N));
    for(i_population = 0; i_population < N; i_population++)
    	scanf("%d", &population[i_population]);
    int i_road, j_road;
    int **road = (int **)malloc((N-1)*sizeof(int *));
    for(i_road = 0; i_road < N-1; i_road++)
    {
    	road[i_road] = (int *)malloc((2)*sizeof(int));
    }
    for(i_road = 0; i_road < N-1; i_road++)
    {
    	for(j_road = 0; j_road < 2; j_road++)
    	{
    		scanf("%d", &road[i_road][j_road]);
    	}
    }
    int Q;
    scanf("%d", &Q);
    int i_cities, j_cities;
    int **cities = (int **)malloc((Q)*sizeof(int *));
    for(i_cities = 0; i_cities < Q; i_cities++)
    {
    	cities[i_cities] = (int *)malloc((3)*sizeof(int));
    }
    for(i_cities = 0; i_cities < Q; i_cities++)
    {
    	for(j_cities = 0; j_cities < 3; j_cities++)
    	{
    		scanf("%d", &cities[i_cities][j_cities]);
    	}
    }

    int* out_ = city_population(N, population, road, Q, cities);
    printf("%d", out_[0]);
    int i_out_;
    for(i_out_ = 1; i_out_ < Q; i_out_++)
    	printf("\n%d", out_[i_out_]);
}

#include <math.h> // for logarithm func

typedef struct _array
{
    int size;
    int* get;
} Array;

void create_adj(int N, int** road, Array** adj)
{
    //create temp array to store number of edges for vertex i
    //calloc to init to 0
    int* counts = (int*) calloc((N+1), sizeof(int));
    counts[0] = 0;
    //find the counts first for N-1 edges
    for (int i = 0; i < N-1; ++i)
    {
        ++(counts[road[i][0]]);
        ++(counts[road[i][1]]);
    }

    //allocate memory for each neighbour array
    for (int i = 0; i < N+1; ++i)
    {
        if (counts[i] <= 0)
            continue;
        (*adj)[i].size = counts[i];
        (*adj)[i].get = (int*) malloc((counts[i]) * sizeof(int));
    }

    //repurpose counts to use as indices to the arrays
    for (int i = 0; i < N+1; ++i)
    {
        --(counts[i]); //make it point to last index
    }

    //then finally fill the array
    for (int i = 0; i < N-1; ++i)
    {
        int u = road[i][0];
        int v = road[i][1];

        if (counts[u] >= 0)
        {
            (*adj)[u].get[(counts[u])--] = v;
        }
        if (counts[v] >= 0)
        {
            (*adj)[v].get[(counts[v])--] = u;
        }
    }
    //return memory
    free(counts);
}

//merge two arrays, in sorted order
void merge_sorted(Array* dest, Array arr1, Array arr2)
{
    int size1 = arr1.size;
    int size2 = arr2.size;
    dest->size = size1 + size2;
    if (dest->size == 0)
        return; //nothing to merge
    dest->get = (int*) malloc(dest->size * sizeof(int));
    int i = 0, j = 0;
    int d = 0;
    while(d < dest->size)
    {
        if (i >= size1)
        {
            //add the rest arr2
            while (j < size2)
                dest->get[d++] = arr2.get[j++];
            return;
        }
        if (j >= size2)
        {
            //add the rest of arr1
            while (i < size1)
                dest->get[d++]= arr1.get[i++];
            return;
        }

        if (arr1.get[i] < arr2.get[j])
            dest->get[d++] = arr1.get[i++];
        else
            dest->get[d++] = arr2.get[j++];
    }
}

//DFS for tree size for HLD, and precompute parent for LCA
void DFS(Array* adj, int** parents, int* level, int* size, int lca_MAX, int v, int par, int l)
{
    parents[v][0] = par;
    //precomputation
    //explanation: parent[v][i] represents teh vertex 2^i jumps from v
    //i.e. p[5][0], what is vertex 5's parent
    //i.e. p[5][1], what is vertex 5's parent's parent
    for (int i=1; i <= lca_MAX; ++i)
    {
        if (parents[v][i-1]!=0) //=0 would mean its jumping past the root's parent
            parents[v][i] = parents[parents[v][i-1]][i-1];
    }
    //set subtree size
    size[v] += 1;
    //set node height/depth
    level[v] = l;
    Array nbrs = adj[v];
    for (int i = 0; i < nbrs.size; ++i)
    {
        if (nbrs.get[i] != par)
        {
            //recursive step on children of v
            DFS(adj,parents,level,size,lca_MAX,nbrs.get[i],v,l+1);
            size[v]+=size[nbrs.get[i]];
        }
    }
}

//Lowest common ancestor
//return lca as return value, and the direct child in lca_children
//0th index will be u's direct child from lca, and 1st will be v
//ASSUMES V is further from root
//int LCA(int** parent, int* level, int lca_MAX, int u, int v, int* lca_children)
//{
//    int level_diff = level[v] - level[u];
//    //move v to where u is
//    for (int i = lca_MAX; i>=0; --i) //start from biggest possible jump
//    {
//        //a bit operation, explanation:
//        //
//        if((level_diff & (1<<i)) != 0)
//            v = parent[v][i];
//    }
//    // this means that u is along path from v to root, LCA is u
//    if (u==v)
//    {
//        lca_children[0] = u;
//        lca_children[1] = v;
//        return u;
//    }
//
//    //next jump both to their parents until parents same
//    for (int i = lca_MAX; i>=0; --i)
//    {
//        if (parent[u][i] != parent[v][i])
//        {
//            u = parent[u][i];
//            v = parent[v][i];
//        }
//    }
//    lca_children[0] = u;
//    lca_children[1] = v;
//    return parent[u][0]; //one last jump to get parent of u and v
//}

//heavy light decomposition
//very long param list but its mostly pass by reference arrays
//similar to DFS in concept
void HLD(int* chain, int* chain_head, int* position, int* subtree_size, int* pop_arr, Array* adj, int* population, int* chain_id, int* pos, int v, int par)
{
    int heavy_child = -1, heavy_size = 0;
    chain[v] = *chain_id;
    position[v] = (*pos)++;

    //first find heaviest child of this node
    Array nbrs = adj[v];
    for (int i = 0; i < nbrs.size; ++i)
    {
        if (nbrs.get[i] != par)
        {
            if (subtree_size[nbrs.get[i]] > heavy_size)
            {
                heavy_size = subtree_size[nbrs.get[i]];
                heavy_child = nbrs.get[i];
            }
        }
    }
    //if there is a heavy child, its still in same chain
    if (heavy_child!=-1)
    {
        pop_arr[*pos] = population[heavy_child-1];
        HLD(chain, chain_head, position, subtree_size, pop_arr, adj, population, chain_id, pos, heavy_child, v);
    }
    //else for the other children, increase chain id and HLD
    nbrs = adj[v]; //just in case
    for (int i = 0; i < nbrs.size; ++i)
    {
        if (nbrs.get[i] != par && nbrs.get[i] != heavy_child)
        {
            (*chain_id)++;
            chain_head[*chain_id] = nbrs.get[i];
            pop_arr[*pos] = population[nbrs.get[i]-1];
            HLD(chain, chain_head, position, subtree_size, pop_arr, adj, population, chain_id, pos, nbrs.get[i], v);
        }
    }
}

//build segment tree
//similar to mergesort
void build_tree(Array** tree, int* pop_arr, int idx, int low, int high)
{
    if (low == high)
    {
        (*tree)[idx].get = (int*) malloc(sizeof(int));
        (*tree)[idx].get[0] = pop_arr[low];
        (*tree)[idx].size = 1;
    }
    else
    {
        int mid = (low+high)>>1; //divide by two using efficient bitshift
        //recurse on half the arr
        build_tree(tree, pop_arr, 2*idx, low, mid);
        build_tree(tree, pop_arr, 2*idx+1, mid+1, high);
        //idx linked list is formed by subtree's linked list
        merge_sorted(&((*tree)[idx]), (*tree)[2*idx], (*tree)[2*idx+1]);
    }
}

//return idx of first element > w
//acts like lower bound in cpp
int binary_search(Array arr, int w)
{
    int start = 0, end = arr.size;
    int mid;
    while(start < end) //terminates when start==end, or when start jumps past end
    {
        //jump to middle
        mid = start + ((end - start) >> 1);
        //if less jump to mid of right side
        //else jump to left
        if (w >= arr.get[mid])
            start = mid + 1;
        else
            end = mid;
    }
    //fix the start jumping past
    if (start < arr.size && arr.get[start] <= w)
        ++start;

    return start;
}

//query segment tree
int query_tree(Array* tree, int idx, int low, int high, int l, int r, int w)
{
    // out of bound
    if (r < low || l > high)
        return 0;

    //complete overlap
    //low and high is within query range
    if (low >= l && high <= r)
    {
        //printf("idx %d, %d to %d (%d) : %d\n", idx, l, r, w, binary_search(tree[idx], w));
        return binary_search(tree[idx], w);
    }

    //else partial overlap
    int mid = (low+high)>>1;
    return query_tree(tree, 2*idx, low, mid, l, r, w)
        + query_tree(tree, 2*idx+1, mid+1, high, l, r, w);
}

//query the HLD tree
int query_HLD(int* chain, int* chain_head, int* position, int** parent, Array* tree, int* depth, int N, int u, int v, int w)
{
    //printf("query %d to %d\n", u, v);
    // v should be smaller, as we expect ito be the lca which is higher up the chain.
    int curr = u;
    int ans = 0;
    //need to traverse to same chain as destination v
    while (chain[curr] != chain[v])
    {
        //make curr's chain be the further one from root
        if (depth[chain_head[chain[curr]]] <= depth[chain_head[chain[v]]])
        {
            int tmp = curr;
            curr = v;
            v = tmp;
        }
        //do a segment tree query from current position to chain head
        ans += query_tree(tree, 1,0,N-1, position[chain_head[chain[curr]]], position[curr], w);
        curr = parent[chain_head[chain[curr]]][0]; //move curr to parent of current chain head
    }
    //now curr should be on same chain as v
    //left with just segment tree query
    //make curr be the further one from root
    if (depth[curr] <= depth[v])
    {
        int tmp = curr;
        curr = v;
        v = tmp;
    }
    ans += query_tree(tree, 1,0,N-1, position[v], position[curr], w);
    return ans;
}

int* city_population (int N, int* population, int** road, int Q, int** cities)
{
    //declare ret array
    int* ret = (int*) malloc(Q * sizeof(int));

    //create adjacency list of graph
    //list of arrays
    Array* adj = (Array*) calloc((N+1) , sizeof(Array));
    create_adj(N, road, &adj);

    //declare lca stuff
    int lca_MAX = (int)(log(N)/log(2)); //max for binary lifting
    int** parent = (int**) malloc((N+1) * sizeof(int*)); //precomputed parents
    for (int i = 0; i < N+1; ++i)
        parent[i] = (int*) calloc((lca_MAX+1), sizeof(int));
    int* graph_depth = (int*) calloc((N+1) , sizeof(int)); //to store depth of each vertex
    //int* lca_children = (int*) calloc((2) , sizeof(int));

    //declare hld stuff
    int* subtree_size = (int*) calloc((N+1) , sizeof(int)); //to store subtree size at vertex
    int* chain = (int*) calloc((N+1) , sizeof(int));
    int* chain_head = (int*) calloc((N+1) , sizeof(int));
    int* position = (int*) calloc((N+1) , sizeof(int));
    int chain_id = 0;
    int pos = 0;

    //declare segment tree stuff
    int* pop_arr = (int*) calloc((N) , sizeof(int)); //store population values in HLD order
    int segment_SIZE = (int)pow(2, ceil(log(N)/log(2)) + 1);
    Array* seg_tree = (Array*) malloc((segment_SIZE) * sizeof(Array));
    for (int i = 0; i < segment_SIZE; ++i)
    {
        seg_tree[i].size = 0;
        seg_tree[i].get = NULL;
    }

    //dfs with 1st node as root
    DFS(adj,parent,graph_depth,subtree_size,lca_MAX, 1,0,0);
    chain_head[0] = 1;
    pop_arr[0] = population[0];
    HLD(chain, chain_head, position, subtree_size, pop_arr, adj, population, &chain_id, &pos, 1,0);
    build_tree(&seg_tree, pop_arr, 1, 0, N-1);
    for (int i = 0; i < Q; ++i)
    {
        int u = cities[i][0];
        int v = cities[i][1];
        int w = cities[i][2];
        //ues HLD to traverse chains and query segment tree, no need lca
        ret[i] = query_HLD(chain, chain_head, position, parent, seg_tree, graph_depth, N, u, v, w);
    }

    return ret;
}
