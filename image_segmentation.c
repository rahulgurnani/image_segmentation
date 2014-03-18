/*
	Image Segmentation using Prim's Algorithm
	Author:-Rahul Gurnani
*/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
//represents each pixel
typedef struct {
    int x, y;
} pixel;
//represents edgethat joins A and B
typedef struct {
    pixel A, B;
    int weight;
} edge;
//adjacency list representation of graph
typedef struct {
    int n;
    edge adj[8];
} graph;
//heap data structure
typedef struct {
    edge* a;
    int n;
} heap;
void give_colors(int* finalR, int* finalG, int* finalB, int splits);		//to give finalcolours after segmentation
int val(int a, int b);														
int inbound(int i, int j, int width, int height);							//to check whether (i,j) is in bound of image
void add_edges(graph **, pixel **, int**, int**, int**, int , int , int , int ) ;	//to add edges
void dfs(graph** , int** , pixel* , int , int , int *, int*) ;					//depth first search for blob colouring
void prim(graph** Gr, heap* h, edge* tree, int** marked, int max, int* size) ;	//prims algorithm to draw MST
int val(int a, int b) {
    return (a >= b) ? a - b : b - a;
}

int inbound(int i, int j, int width, int height) {
    if (i < height && i >= 0 && j < width && j >= 0)return 1;
    else return 0;
}

void add_edges(graph **Gr, pixel **image, int**R, int**G, int**B, int i, int j, int width, int height) {
    int n = 0, rt, gt, bt, r, g, b;
    r = R[i][j];
    g = G[i][j];
    b = B[i][j];

    if (inbound(i - 1, j - 1, width, height))//1
    {
        rt = R[i - 1][j - 1];
        gt = G[i - 1][j - 1];
        bt = B[i - 1][j - 1];

        Gr[i][j].adj[n].A = image[i - 1][j - 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i - 1, j, width, height))//2
    {
        rt = R[i - 1][j ];
        gt = G[i - 1][j ];
        bt = B[i - 1][j ];

        Gr[i][j].adj[n].A = image[i - 1][j];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i - 1, j + 1, width, height))//3
    {
        rt = R[i - 1][j + 1];
        gt = G[i - 1][j + 1];
        bt = B[i - 1][j + 1];
        Gr[i][j].adj[n].A = image[i - 1][j + 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i, j + 1, width, height))//4
    {
        rt = R[i][j + 1];
        gt = G[i][j + 1];
        bt = B[i][j + 1];

        Gr[i][j].adj[n].A = image[i][j + 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i + 1, j + 1, width, height))//5
    {
        rt = R[i + 1][j + 1];
        gt = G[i + 1][j + 1];
        bt = B[i + 1][j + 1];

        Gr[i][j].adj[n].A = image[i + 1][j + 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i + 1, j, width, height))//6
    {
        rt = R[i + 1][j ];
        gt = G[i + 1][j ];
        bt = B[i + 1][j ];

        Gr[i][j].adj[n].A = image[i + 1][j];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i + 1, j - 1, width, height))//7
    {
        rt = R[i + 1][j - 1];
        gt = G[i + 1][j - 1];
        bt = B[i + 1][j - 1];

        Gr[i][j].adj[n].A = image[i + 1][j - 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    if (inbound(i, j - 1, width, height))//8
    {
        rt = R[i][j - 1];
        gt = G[i][j - 1];
        bt = B[i][j - 1];

        Gr[i][j].adj[n].A = image[i][j - 1];
        Gr[i][j].adj[n].B = image[i][j];
        Gr[i][j].adj[n].weight = val(r, rt) + val(g, gt) + val(b, bt);
        n++;
    }
    Gr[i][j].n = n;
}

void dfs(graph** Gr, int** marked, pixel* s, int x, int y, int *size, int*color) {
    marked[x][y] = *color;
    pixel p;
    int k, i;
    s[*size].x = x;
    s[*size].y = y;
    (*size)++;
    while (*size > 0) {
        (*size)--;
        p = s[*size];
        k = Gr[p.x][p.y].n;
        for (i = 0; i < k; i++) {
            if (marked[Gr[p.x][p.y].adj[i].A.x][Gr[p.x][p.y].adj[i].A.y] != 0)continue;
            else {
                marked[Gr[p.x][p.y].adj[i].A.x][Gr[p.x][p.y].adj[i].A.y] = *color;
                s[*size].x = Gr[p.x][p.y].adj[i].A.x;
                s[*size].y = Gr[p.x][p.y].adj[i].A.y;
                (*size)++;
            }
        }
        marked[p.x][p.y] = *color;
    }
}

void prim(graph** Gr, heap* h, edge* tree, int** marked, int max, int* size) {
    visit(Gr, h, marked, 0, 0);
    edge e;
    while (h->n > 0 && *size < max - 1) {
        e = del_min(h);
        if ((marked[e.A.x][e.A.y] == 1)&&(marked[e.B.x][e.B.y] == 1))continue;
        tree[*size] = e;
        (*size)++;
        if (marked[e.A.x][e.A.y] == 0)visit(Gr, h, marked, e.A.x, e.A.y);
        if (marked[e.B.x][e.B.y] == 0)visit(Gr, h, marked, e.B.x, e.B.y);

    }
}

void visit(graph** Gr, heap* h, int** marked, int x, int y) {
    marked[x][y] = 1;
    int k = Gr[x][y].n;
    int i;
    for (i = 0; i < k; i++) {
        if (marked[Gr[x][y].adj[i].A.x][Gr[x][y].adj[i].A.y] == 0)insert_minheap(h, Gr[x][y].adj[i]);
    }
}

void insert_minheap(heap*h, edge e) {
    int i = h->n;
    h->a[i] = e;
    h->n++;
    while (i > 0 && (h->a[(i - 1) / 2].weight > h->a[i].weight)) {
        swap(&(h->a[(i - 1) / 2]), &(h->a[i]));
        i = (i - 1) / 2;
    }
}

void swap(edge*a, edge * b) {
    edge temp = *a;
    *a = *b;
    *b = temp;
    return;
}

edge del_min(heap * h) {
    edge e = h->a[0];
    h->n--;
    swap(&(h->a[0]), &(h->a[h->n]));
    min_heapify(h, 0);
    return e;
}

void min_heapify(heap* h, int i) {
    int l, r, smallest;
    while (1) {
        l = i + i + 1;
        r = i + i + 2;
        if (l < h->n && h->a[l].weight < h->a[i].weight)
            smallest = l;
        else smallest = i;
        if (r < h->n && h->a[r].weight < h->a[smallest].weight)
            smallest = r;
        if (smallest == i)
            return;
        swap(&(h->a[i]), &(h->a[smallest]));
        i = smallest;
    }
}

void make_maxheap(heap * p) {
    int i;
    for (i = (p->n) / 2 - 1; i >= 0; i--) {
        max_heapify(p, i);
    }
}

void max_heapify(heap* h, int i) {
    int l, r, largest;
    while (1) {
        l = i + i + 1;
        r = i + i + 2;
        if (l < h->n && h->a[l].weight > h->a[i].weight)
            largest = l;
        else largest = i;
        if (r < h->n && h->a[r].weight > h->a[largest].weight)
            largest = r;
        if (largest == i)
            return;
        swap(&(h->a[i]), &(h->a[largest]));
        i = largest;
    }
}

edge del_max(heap * h) {
    edge e = h->a[0];
    h->n--;
    swap(&(h->a[0]), &(h->a[h->n]));
    max_heapify(h, 0);
    return e;
}

void build_graph(heap*p, graph**L, int width, int height) {
    int i, j;
    pixel temp;
    for (i = 0; i < height; i++)
        for (j = 0; j < width; j++) {
            L[i][j].n = 0;
        }
    for (i = 0; i < p->n; i++) {

        L[p->a[i].B.x][p->a[i].B.y].adj[L[p->a[i].B.x][p->a[i].B.y].n] = p->a[i];
        L[p->a[i].B.x][p->a[i].B.y].n++;
        temp = p->a[i].A;
        p->a[i].A = p->a[i].B;
        p->a[i].B = temp;
        L[p->a[i].B.x][p->a[i].B.y].adj[L[p->a[i].B.x][p->a[i].B.y].n] = p->a[i];
        L[p->a[i].B.x][p->a[i].B.y].n++;
    }
}

void initialize(int** a, int width, int height) {
    int i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            a[i][j] = 0;
        }
    }
}

void give_colors(int* finalR, int* finalG, int* finalB, int splits) {
    int i;
    for (i = 0; i < splits; i++) {
        finalR[i] = rand() % 256;
        finalG[i] = rand() % 256;
        finalB[i] = rand() % 256;
    }
}


int main() {
    srand((unsigned int) time(NULL));
    FILE *fp;
    fp = fopen("original.ppm", "r");
    char magic[3];
    fscanf(fp, "%s", magic);
    int width, height, i, j, max_value;
    fscanf(fp, "%d %d", &width, &height);
    fscanf(fp, "%d", &max_value);

    pixel **image;
    image = (pixel**) malloc(sizeof (pixel*) * height);
    for (i = 0; i < height; i++) {
        image[i] = (pixel*) malloc(sizeof (pixel) * width);
    }

    int** marked = (int**) malloc(sizeof (int*) * height);
    for (i = 0; i < height; i++) {
        marked[i] = (int*) malloc(sizeof (int) * width);
    }

    int** R = (int**) malloc(sizeof (int*) * height);
    int** G = (int**) malloc(sizeof (int*) * height);
    int** B = (int**) malloc(sizeof (int*) * height);
    for (i = 0; i < height; i++) {
        R[i] = (int*) malloc(sizeof (int) * width);
        G[i] = (int*) malloc(sizeof (int) * width);
        B[i] = (int*) malloc(sizeof (int) * width);
    }

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            fscanf(fp, " %d %d %d", &R[i][j], &G[i][j], &B[i][j]);
            image[i][j].x = i;
            image[i][j].y = j;
        }
    }

    fclose(fp);

    int max = width*height;
    graph** Gr = (graph**) malloc(sizeof (graph*) * height);
    for (i = 0; i < height; i++) {
        Gr[i] = (graph*) malloc(sizeof (graph) * width);
    }
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            add_edges(Gr, image, R, G, B, i, j, width, height);
        }
    }


    initialize(marked, width, height);

    edge* tree = (edge*) malloc(sizeof (edge)*(max));
    heap* h = (heap*) malloc(sizeof (heap));
    h->a = (edge*) malloc(sizeof (edge)*8 * max);
    h->n = 0;
    int size = 0;
    prim(Gr, h, tree, marked, max, &size);

    heap* p = (heap*) malloc(sizeof (heap));
    p->a = tree;
    p->n = size;

    make_maxheap(p);

    int splits;
    printf("No. of segments = ");			//number of segments splits
    scanf("%d",&splits);
    for (i = 0; i < splits - 1; i++) {
        del_max(p);
        size--;
    }

    graph** L = (graph**) malloc(sizeof (graph*) * height);
    for (i = 0; i < height; i++) {
        L[i] = (graph*) malloc(sizeof (graph) * width);
    }
    build_graph(p, L, width, height);

    print_graph(L,width,height);

    initialize(marked, width, height);
    pixel* d = (pixel*) malloc(sizeof (pixel) * max);
    int depth = 0, color = 1;
    for (i = 0; i < height; i++)
        for (j = 0; j < width; j++)
            if (marked[i][j] == 0) {
                dfs(L, marked, d, i, j, &depth, &color);
                color++;
            }
    int* finalR = (int*) malloc(sizeof (int) * splits);
    int* finalG = (int*) malloc(sizeof (int) * splits);
    int* finalB = (int*) malloc(sizeof (int) * splits);
    give_colors(finalR, finalG, finalB, splits);
    FILE* fo;
    fo = fopen("output5.ppm", "w");
    fprintf(fo, "%s\n%d %d\n%d", magic, width, height, max_value);
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            fprintf(fo, "\n%d\n%d\n%d", finalR[marked[i][j]-1], finalG[marked[i][j]-1], finalB[marked[i][j]-1]);
        }
    }
    fclose(fo);
    
    for(i=0;i<height;i++)
    {
        free(image[i]);free(marked[i]);free(R[i]);free(G[i]);free(B[i]);free(Gr[i]);free(L[i]);
    }
    free(image);free(marked);free(R);free(G);free(B);free(Gr);free(L);free(tree);free(h->a);free(h);free(p);free(d);free(finalR);
    free(finalG);free(finalB);
    
    return 0;
}












