#include <stdio.h>
#include <stdlib.h>

// Author: Terence
// This code does Ford-Fulkerson Algorithm with BFS,
// to do matching of Student to Project to Teacher.
// Date structures and helper functions boilerplate were provided.




typedef struct _listnode {
  int vertex;
  struct _listnode *next;
} ListNode;
typedef ListNode StackNode;
typedef ListNode QueueNode;

typedef struct _graph {
  int V;
  int E;
  //ListNode **list;
  int** adj;
} Graph;

typedef struct _queue {
  int size;
  QueueNode *head;
  QueueNode *tail;
} Queue;

typedef struct _stack {
  int size;
  StackNode *head;
} Stack;


void insertAdjVertex(ListNode **AdjList, int vertex);
void removeAdjVertex(ListNode **AdjList, int vertex);
int matching(Graph g);

void enqueue(Queue *qPtr, int item);
int dequeue(Queue *qPtr);
int getFront(Queue q);
int isEmptyQueue(Queue q);
void removeAllItemsFromQueue(Queue *qPtr);
void printQ(QueueNode *cur);
//////STACK///////////////////////////////////////////
void push(Stack *sPtr, int vertex);
int pop(Stack *sPtr);
int peek(Stack s);
int isEmptyStack(Stack s);
void removeAllItemsFromStack(Stack *sPtr);
//////////////////////////////////

// helpers
void print_graph(Graph g) {
  printf("Printing Graph\nNum Edges = %d\n", g.E);
  for (int i = 0; i < g.V; ++i) {
    printf("%d: ", i);
    int* arr = g.adj[i];
    for (int j = 0; j < g.V; ++j)
    {
      if (arr[j] > 0)
        printf("%d, ", j);
    }
    printf("\n");
  }
}
void copy_graph(Graph *dest, Graph source) {
  dest->E = source.E;
  dest->V = source.V;
  dest->adj = (int**) malloc(source.V * sizeof(int*));
  for (int i = 0; i < source.V; ++i)
  {
    dest->adj[i] = (int*) malloc(source.V * sizeof(int));
    for (int j = 0; j < source.V; ++j)
    {
      dest->adj[i][j] = source.adj[i][j];
    }
  }
}
//check if edge u,v is in graph
// int edge_in_graph(Graph g, int u, int v) {
//    ListNode* tmp = g.list[u];
//    while(tmp != NULL)
//    {
//       if (tmp->vertex == v)
//          return 1;
//       tmp = tmp->next;
//    }
//    return 0;
// }

int main() {
  int Prj, Std, Mtr; // Project, Student and Mentor;
  int maxMatch;
  scanf("%d %d %d", &Std, &Prj, &Mtr);
  int np, nm; // number of projects and number of mentors

  // build graph
  Graph g;
  // Write your code
  g.E = Std + Prj * Mtr; // source + sink edges
  // 0 and last idx will be source and sink
  // Prj->Std->Std->Mtr
  g.V = 1 + Prj + Std + Std + Mtr + 1;
  int prjIdx = 1;
  int std1Idx = Prj + 1;
  int std2Idx = Prj + Std + 1;
  int mtrIdx = Prj + Std + Std + 1;
  int sinkIdx = g.V - 1;
  //g.list = (ListNode**)malloc(g.V * sizeof(ListNode*));
  //g.list[g.V - 1] = NULL; //sink will have no outgoing edges
  g.adj = (int**) malloc(g.V * sizeof(int*));
  for (int i = 0; i < g.V; ++i)
    g.adj[i] = (int*) calloc(g.V, sizeof(int));
  // create adjanceny matrix
  for (int i = 0; i < Std; ++i) {
    scanf("%d %d", &np, &nm);
    //from project to student1
    for (int j = 0; j < np; ++j) {
      int p = 0;
      scanf("%d", &p);
      g.adj[prjIdx+p-1][std1Idx+i] = 1;
    }
    //from student1 to student2
    g.adj[std1Idx+i][std2Idx+i] = 1;
    //from student2 to mentor
    for (int j = 0; j < nm; ++j) {
      int m = 0;
      scanf("%d", &m);
      g.adj[std2Idx+i][mtrIdx+m-1] = 1;
    }
  }
  // make mentors point to sink
  for (int i = 0; i < Mtr; ++i) {
    g.adj[mtrIdx+i][sinkIdx] = 1;
  }
  //add source to projects
  for (int i = 0; i < Prj; ++i) {
    g.adj[0][prjIdx+i] = 1;
  }
  //print_graph(g);
  // apply Ford Fulkerson algorithm
  // use DFS or BFS to find a path
  maxMatch = matching(g);
  printf("%d\n", maxMatch);

  return 0;
}

// if reached destination i.e. sink, return 1
int find_aug_path(Graph res_g, int* parent) {
  Queue q;
  q.head = q.tail = NULL;
  q.size = 0;
  enqueue(&q, 0);
  int target = res_g.V - 1;
  int *visited = (int *)calloc(res_g.V, sizeof(int));
  visited[0] = 1;
  while (!isEmptyQueue(q)) {
    int front = getFront(q);
    dequeue(&q);
    // get the neighbours and add to queue
    int* nbrs = res_g.adj[front];
    for (int i = 0; i < res_g.V; ++i) {
      if (nbrs[i] <= 0) //if there is path
        continue;
      int nbr = i;
      if (nbr == target) 
      {
         parent[nbr] = front;
         return 1;
      }
      if (!visited[nbr])
      {
         parent[nbr] = front;
         visited[nbr] = 1;
         enqueue(&q, nbr);
      }
    }
  }
  free(visited);
  return 0;
}

int matching(Graph g) {
  if (g.adj == NULL)
    return 0;
  // create residual graph
  int count = 0;
  Graph res_g;
  copy_graph(&res_g, g);
  //print_graph(res_g);
  int* parent = (int*) malloc(g.V * sizeof(int));
  parent[0] = -1;
  // do BFS find augmenting path
  while (find_aug_path(res_g, parent)) {
    int s = 0, t = g.V-1; //source, sink 
    //printf("Path: ");
    //reverse path in res_g and find min flow
    int minflow = 999999;
    while(s != t)
    {
      //printf("<-%d",t);
      int flowf = res_g.adj[parent[t]][t];
      if (flowf < minflow)
        minflow = flowf;
      t = parent[t]; //traverse path
    }
    //printf("\n");

    s = 0, t = g.V-1; //source, sink
    while(s != t)
    {
      res_g.adj[parent[t]][t] -= minflow;
      res_g.adj[t][parent[t]] += minflow;
      t = parent[t]; //traverse path
    }
    //print_graph(res_g);
    count += minflow;
  }
  return count;
}

void removeAdjVertex(ListNode **AdjList, int vertex) {
  ListNode *temp, *preTemp;
  if (*AdjList != NULL) {
    if ((*AdjList)->vertex == vertex) // first node
    {
      temp = *AdjList;
      *AdjList = (*AdjList)->next;
      free(temp);
      return;
    }
    preTemp = *AdjList;
    temp = (*AdjList)->next;
    while (temp != NULL && temp->vertex != vertex) {
      preTemp = temp;
      temp = temp->next;
    }
    preTemp->next = temp->next;
    free(temp);
  }
}

void insertAdjVertex(ListNode **AdjList, int vertex) {
  if (*AdjList == NULL) {
    *AdjList = (ListNode *)malloc(sizeof(ListNode));
    (*AdjList)->vertex = vertex;
    (*AdjList)->next = NULL;
  } else {
    ListNode *temp;
    temp = (ListNode *)malloc(sizeof(ListNode));
    temp->vertex = vertex;
    temp->next = *AdjList;
    *AdjList = temp;
  }
}

/*BELOW HERE IS BOILERPLATE, provided functions*////////////////////////

void enqueue(Queue *qPtr, int vertex) {
  QueueNode *newNode;
  newNode = (QueueNode *)malloc(sizeof(QueueNode));
  if (newNode == NULL)
    exit(0);

  newNode->vertex = vertex;
  newNode->next = NULL;

  if (isEmptyQueue(*qPtr))
    qPtr->head = newNode;
  else
    qPtr->tail->next = newNode;

  qPtr->tail = newNode;
  qPtr->size++;
}

int dequeue(Queue *qPtr) {
  if (qPtr == NULL || qPtr->head == NULL) // Queue is empty or NULL pointer
  {
    return 0;
  } else {
    QueueNode *temp = qPtr->head;
    qPtr->head = qPtr->head->next;
    if (qPtr->head == NULL) // Queue is emptied
      qPtr->tail = NULL;

    free(temp);
    qPtr->size--;
    return 1;
  }
}

int getFront(Queue q) { return q.head->vertex; }

int isEmptyQueue(Queue q) {
  if (q.size == 0)
    return 1;
  else
    return 0;
}

void removeAllItemsFromQueue(Queue *qPtr) {
  while (dequeue(qPtr))
    ;
}

void printQ(QueueNode *cur) {
  if (cur == NULL)
    printf("Empty");

  while (cur != NULL) {
    printf("%d ", cur->vertex);
    cur = cur->next;
  }
  printf("\n");
}

void push(Stack *sPtr, int vertex) {
  StackNode *newNode;
  newNode = (StackNode *)malloc(sizeof(StackNode));
  newNode->vertex = vertex;
  newNode->next = sPtr->head;
  sPtr->head = newNode;
  sPtr->size++;
}

int pop(Stack *sPtr) {
  if (sPtr == NULL || sPtr->head == NULL) {
    return 0;
  } else {
    StackNode *temp = sPtr->head;
    sPtr->head = sPtr->head->next;
    free(temp);
    sPtr->size--;
    return 1;
  }
}

int isEmptyStack(Stack s) {
  if (s.size == 0)
    return 1;
  else
    return 0;
}

int peek(Stack s) { return s.head->vertex; }

void removeAllItemsFromStack(Stack *sPtr) {
  while (pop(sPtr))
    ;
}
