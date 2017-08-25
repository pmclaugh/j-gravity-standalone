#include "standalone.h"
#include <math.h>
#include <limits.h>

int sbod_comp(const void *a, const void *b)
{
    uint64_t ma = ((t_sortbod *)a)->morton;
    uint64_t mb = ((t_sortbod *)b)->morton;
    if (ma > mb)
        return 1;
    if (ma == mb)
        return 0;
    else
        return -1;
}

#define SORT_NAME mort
#define SORT_TYPE t_sortbod
#define SORT_CMP(x, y) (sbod_comp(&x, &y))
#include "sort.h"

#define THETA 1
#define LEAF_THRESHOLD pow(2, 10)
#define THREADCOUNT 8


/*
    utility and debug functions
*/
void print_bounds(t_bounds b)
{
    printf("bounds: %f %f %f %f %f %f\n", b.xmin, b.xmax, b.ymin, b.ymax, b.zmin, b.zmax);
}

int body_in_bounds(t_body body, t_bounds bounds)
{
    if (body.position.x > bounds.xmax || body.position.x < bounds.xmin)
        return (0);
    if (body.position.y > bounds.ymax || body.position.y < bounds.ymin)
        return (0);
    if (body.position.z > bounds.zmax || body.position.z < bounds.zmin)
        return (0);
    return (1);
}

t_tree *new_tnode(cl_float4 *positions, cl_float4 *velocities, int count, t_tree *parent)
{
    t_tree *node = (t_tree *)calloc(1, sizeof(t_tree));
    //node->bodies = bodies;
    node->positions = positions;
    node->velocities = velocities;
    node->count = count;
    node->parent = parent;
    node->children = NULL;
    return (node);
}

int node_depth(t_tree *node)
{
    int depth = 0;
    while (node->parent)
    {
        node = node->parent;
        depth++;
    }
    return depth;
}

/*
    generation of morton codes
*/

// method to seperate bits from a given integer 3 positions apart
//from http://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
uint64_t splitBy3(const unsigned int a)
{
    uint64_t x = a & 0x1fffff; // we only look at the first 21 bits
    x = (x | x << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 16) & 0x1f0000ff0000ff;  // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}
 
uint64_t mortonEncode_magicbits(const unsigned int x, const unsigned int y, const unsigned int z)
{
	//interweave the rightmost 21 bits of 3 unsigned ints to generate a 63-bit morton code
    uint64_t answer = 0;
    answer |= splitBy3(z) | splitBy3(y) << 1 | splitBy3(x) << 2;
    return answer;
}

uint64_t morton64(float x, float y, float z)
{
	//x, y, z are in [0, 1]. multiply by 2^21 to get 21 bits, confine to [0, 2^21 - 1]
    x = fmin(fmax(x * 2097152.0f, 0.0f), 2097151.0f);
    y = fmin(fmax(y * 2097152.0f, 0.0f), 2097151.0f);
    z = fmin(fmax(z * 2097152.0f, 0.0f), 2097151.0f);
    return (mortonEncode_magicbits((unsigned int)x, (unsigned int)y, (unsigned int)z));
}

uint64_t mortonize(cl_float4 pos, t_bounds bounds)
{
    float distance =  1.0 / (bounds.xmax - bounds.xmin);
    return morton64((pos.x - bounds.xmin) * distance, (pos.y - bounds.ymin) * distance, (pos.z - bounds.zmin) * distance);
}

///Center of gravity calculations///

static cl_float4 vadd(cl_float4 a, cl_float4 b)
{
    //add two vectors.
    return ((cl_float4){a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w});
}

static cl_float4 center_add(cl_float4 total, cl_float4 add)
{
	//method for tallying running total for center of gravity
    add.w = fabs(add.w);
    add.x *= add.w;
    add.y *= add.w;
    add.z *= add.w;
    return (cl_float4){total.x + add.x, total.y + add.y, total.z + add.z, total.w + add.w};
}

static t_body COG_from_children(t_tree **children)
{
    // we calculate centers of gravity starting at leaf cells.
    // therefore if a cell has children, all their COGs have been calculated.
    // it is much faster to just use 8 values than potentially millions.
    // this reduces COG calculate from O(nlog(n)) to very nearly O(n).
	cl_float4 center = (cl_float4){0,0,0,0};
    float real_total = 0;
    for (int i = 0; i < 8; i++)
    {
        center.x += children[i]->as_single->positions[0].x * children[i]->as_single->velocities[0].w;
        center.y += children[i]->as_single->positions[0].y * children[i]->as_single->velocities[0].w;
        center.z += children[i]->as_single->positions[0].z * children[i]->as_single->velocities[0].w;
        center.w += children[i]->as_single->velocities[0].w;
        real_total += children[i]->as_single->positions[0].w;
    }
    t_body b;
    b.velocity = (cl_float4){0, 0, 0, center.w};
    center.x /= center.w;
    center.y /= center.w;
    center.z /= center.w;
    center.w = real_total;
    b.position = center;
    return b;
}

static t_body COG_from_bodies(cl_float4 *positions, cl_float4 *velocities, int count)
{
    cl_float4 center = (cl_float4){0,0,0,0};
    if (count == 0)
        return (t_body){center, center};
    float real_total = 0;
    for (int i = 0; i < count; i++)
    {
        center = center_add(center, positions[i]);
        real_total += positions[i].w;
    }
    t_body b;
    b.velocity = (cl_float4){0, 0, 0, center.w};
    center.x /= center.w;
    center.y /= center.w;
    center.z /= center.w;
    center.w = real_total;
    b.position = center;
    return b;
}

static t_tree *make_as_single(t_tree *c)
{
    t_body b;
    if (!c->children)
	   b = COG_from_bodies(c->positions, c->velocities, c->count);
    else
        b = COG_from_children(c->children);
    t_tree *s = calloc(1, sizeof(t_tree));
    s->parent = NULL;
    s->children = NULL;
    s->count = 1;
    s->as_single = NULL;
    s->positions = calloc(1, sizeof(cl_float4));
    s->velocities = calloc(1, sizeof(cl_float4));
    s->positions[0] = b.position;
    s->velocities[0] = b.velocity;
    return (s);
}

static int count_tree_array(t_tree **arr)
{
    int count;
     if (!arr)
        return 0;
    for (count = 0; arr[count]; count++)
        ;
    return (count);
 }

typedef struct s_list{
    struct s_list *next;
    void *data;
    int id;
}               t_list;

static void lst_push(t_list **lst, t_tree *cell)
{
    t_list *l = calloc(1, sizeof(t_list));
    l->data = cell;
    if (*lst)
    {
        l->next = *lst;
        l->id = (*lst)->id + 1;
    }
    else
        l->next = NULL;
    *lst = l;
}

static void enumerate_worker(t_tree *root, t_list **q)
{
    if (!root)
        return ;
    if (!root->children)
    {
        if (root->count)
            lst_push(q, root);
    }
    else
        for (int i = 0; i < 8; i++)
            enumerate_worker(root->children[i], q);
}

static t_tree **enumerate_better(t_tree *root)
{
    t_list *q = NULL;
    enumerate_worker(root, &q);
    int count = q->id + 1;
    printf("enum_better found %d leaves\n", count);
    t_tree **ret = calloc(count + 1, sizeof(t_tree *));
    for (int i = 0; i < count; i++)
    {
        ret[i] = q->data;
        t_list *tmp = q;
        q = q->next;
        free(tmp);
    }
    return ret;
}

static t_tree **enumerate_leaves(t_tree *root)
{
    //return a linear t_tree** that's all the leaf nodes (ie childless nodes) in the tree
    //this is an excellent opportunity to very quickly do centers of gravity/as_single
    t_tree **ret;

    if (!root->children)
    {
        //root->as_single = make_as_single(root);
        ret = (t_tree **)calloc(2, sizeof(t_tree *));
        ret[0] = root->count ? root : NULL; //we do not bother enumerating empty leaves (empty cells cant have children so this covers all)
        ret[1] = NULL;
        return (ret);
    }
    t_tree ***returned = (t_tree ***)calloc(8, sizeof(t_tree **));
    int total = 0;
    for (int i = 0; i < 8; i++)
    {
        returned[i] = enumerate_leaves(root->children[i]);
        total += count_tree_array(returned[i]);
    }
    //root->as_single = make_as_single(root);
    ret = (t_tree **)calloc(total + 1, sizeof(t_tree *));
    for (int i = 0; i < total;)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; returned[j][k]; k++, i++)
            {
                ret[i] = returned[j][k];
            }
            free(returned[j]);
        }
        free(returned);
    }
    ret[total] = NULL;
    return (ret);
}

/*
    assemble neighborhood and related functions. This is where workunits are built.
*/

static cl_float4 midpoint_from_bounds(t_bounds b)
{
    return (cl_float4){(b.xmax - b.xmin) / 2, (b.ymax - b.ymin) / 2, (b.zmax - b.zmin) / 2};
}

static float multipole_acceptance_criterion(t_tree *us, t_tree *them)
{
    //assess whether a cell is "near" or "far" for the sake of barnes-hut
    //if the value returned is less than THETA, that cell is far
    float s;
    float d;
    cl_float4 r;
    cl_float4 us_midpoint;

    if (us == them)
        return (THETA);
    us_midpoint = midpoint_from_bounds(us->bounds);
    s = them->bounds.xmax - them->bounds.xmin;
    r.x = them->as_single->positions[0].x - us_midpoint.x;
    r.y = them->as_single->positions[0].y - us_midpoint.y;
    r.z = them->as_single->positions[0].z - us_midpoint.z;
    d = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);

    //in normal Barnes-Hut, the MAC is evaluated for every body vs every cell. 
    //we're evaluating it cell to cell for speed and to build good workunits.
    //this could result in some comparisons that should be near being done as far.
    //subtracting half of our cell's width from the measured distance compensates for this adequately.
    //basically, we are measuring as if all the bodies in our cell are right up against the nearest side of our cell.
    d -= (us->bounds.xmax - us->bounds.xmin) / 2;
    return (s/d);
}

static void assemble_worker(t_tree *cell, t_tree *root, t_list **q)
{
    if (!root)
        return ;
    if (root->parent && multipole_acceptance_criterion(cell, root) < THETA)
    {
        lst_push(q, root->as_single);
    }
    else
    {
        if (root->children)
            for (int i = 0; i < 8; i++)
                assemble_worker(cell, root->children[i], q);
        else if (root->count > 0)
            lst_push(q, root);
    }
}

static t_tree **assemble_better(t_tree *cell, t_tree *root)
{
    t_list *q = NULL;
    assemble_worker(cell, root, &q);
    int count = q->id + 1; //could even keep track of particle counts too
    t_tree **ret = calloc(count + 1, sizeof(t_tree *));
    for (int i = 0; i < count; i++)
    {
        ret[i] = q->data;
        t_list *tmp = q;
        q = q->next;
        free(tmp);
    }
    return ret;
}

size_t nearest_256(size_t n)
{
    return (((n / 256) + 1) * 256);
}

t_workunit *cell_to_wu(t_tree *cell, int idx)
{
    t_workunit *wu = calloc(1, sizeof(t_workunit));
    wu->idx = idx;
    wu->N = calloc(nearest_256(cell->count), sizeof(cl_float4));
    wu->V = calloc(nearest_256(cell->count), sizeof(cl_float4));
    memcpy(wu->N, cell->positions, cell->count * sizeof(cl_float4));
    memcpy(wu->V, cell->velocities, cell->count * sizeof(cl_float4));
    wu->n_count = cell->count;
    int mcount = 0;
    for (int i = 0; cell->neighbors[i]; i++)
        mcount += cell->neighbors[i]->count;
    wu->M = calloc(nearest_256(mcount), sizeof(cl_float4));
    mcount = 0;
    for (int i = 0; cell->neighbors[i]; i++)
    {
        memcpy(&(wu->M[mcount]), cell->neighbors[i]->positions, cell->neighbors[i]->count * sizeof(cl_float4));
        mcount += cell->neighbors[i]->count;
    }
    wu->m_count = mcount;
    return (wu);
}

typedef struct s_assemble_set
{
    t_tree **leaves;
    t_tree *root;
    unsigned int offset;
    unsigned int max;
    t_queue *todo_work;
    sem_t   *calc_thread_sem;
}               t_assemble;

void *assemble_thread(void *param)
{
    t_assemble *set = (t_assemble *)param;
    //printf("offset is %d\n", set->offset);
    for (int i = 0; i < set->max; i++)
    {
        set->leaves[set->offset + i]->neighbors = assemble_better(set->leaves[set->offset + i], set->root);
        //serialize to workunit
        t_workunit *wu = cell_to_wu(set->leaves[set->offset + i], set->offset + i);
        //push on queue
        queue_enqueue(&set->todo_work, queue_create_new(wu));
        //post on semaphore
        sem_post(set->calc_thread_sem);
    }
    //printf("thread with offset %d finished\n", set->offset);
    free(set);
    return (0);
}

void mt_assemble_neighborhoods(t_tree **leaves, t_tree *root, t_standalone *sim)
{
    int count = count_tree_array(leaves);
    int lpt = ceil((float)count / THREADCOUNT);
    //printf("lpt is %d\n", lpt);
    pthread_t *assemblers = calloc(THREADCOUNT, sizeof(pthread_t));
    for (int i = 0; i < THREADCOUNT; i++)
    {
        t_assemble *set = calloc(1, sizeof(t_assemble));
        set->leaves = leaves;
        set->root = root;
        set->offset = i * lpt;
        set->max = count - set->offset > lpt ? lpt : count - set->offset;
        set->todo_work = sim->todo_work;
        set->calc_thread_sem = sim->calc_thread_sem;
        pthread_create(&assemblers[i], NULL, assemble_thread, set);
    }
    //printf("threads made\n");
    for (int i = 0; i < THREADCOUNT; i++)
        pthread_join(assemblers[i], NULL);
    //printf("threads joined\n");
    free(assemblers);
}

static t_bounds bounds_from_bodies(cl_float4 *positions, int count)
{
    //at the start of making the tree, we need a box that bounds all the bodies.
    float xmin = 0, xmax = 0;
    float ymin = 0, ymax = 0;
    float zmin = 0, zmax = 0;

    for (int i = 0; i < count; i++)
    {
        //expand bounds if needed
        if (positions[i].x < xmin)
            xmin = positions[i].x;
        if (positions[i].x > xmax)
            xmax = positions[i].x;

        if (positions[i].y < ymin)
            ymin = positions[i].y;
        if (positions[i].y > ymax)
            ymax = positions[i].y;

        if (positions[i].z < zmin)
            zmin = positions[i].z;
        if (positions[i].z > zmax)
            zmax = positions[i].z;
    }

    //bounds must be a cube.
    float min = xmin;
    float max = xmax;
    if (ymin < min)
        min = ymin;
    if (zmin < min)
        min = zmin;
    if (ymax > max)
        max = ymax;
    if (zmax > max)
        max = zmax;
    //printf("bounds from bodies were %f %f %f %f %f %f\n", min, max, min, max, min, max);
    return ((t_bounds){min, max, min, max, min, max});
}

/*
    Creation/division of octree
*/

#define xmid parent.xmax - (parent.xmax - parent.xmin) / 2
#define ymid parent.ymax - (parent.ymax - parent.ymin) / 2
#define zmid parent.zmax - (parent.zmax - parent.zmin) / 2

t_bounds bounds_from_code(t_bounds parent, unsigned int code)
{
    //generate the appropriate subbounds of the parent cell for this 3-bit string
    t_bounds bounds;
    if(code >> 2)
    {
        bounds.xmin = xmid;
        bounds.xmax = parent.xmax;
    }
    else
    {
        bounds.xmin = parent.xmin;
        bounds.xmax = xmid;
    }
    code %= 4;
    if (code >> 1)
    {
        bounds.ymin = ymid;
        bounds.ymax = parent.ymax;
    }
    else
    {
        bounds.ymin = parent.ymin;
        bounds.ymax = ymid;
    }
    code %= 2;
    if (code)
    {
        bounds.zmin = zmid;
        bounds.zmax = parent.zmax;
    }
    else
    {
        bounds.zmin = parent.zmin;
        bounds.zmax = zmid;
    }
    return (bounds);
}

unsigned int mask(uint64_t mort, int depth)
{
    mort = mort << (1 + 3 * depth);
    mort = mort >> 61;
    return (unsigned int)mort;
}

void new_split(t_tree *node)
{
    int depth = node_depth(node);

    int *counts = calloc(8, sizeof(int));
    for (int i = 0; i < node->count; i++)
        counts[mask(node->mortons[i], depth)]++;
    int *spots = calloc(8, sizeof(int));
    for (int i = 1; i < 8; i++)
        spots[i] = spots[i - 1] + counts[i - 1];
    //this should offer some pretty substantial speedup so i think it's ok to be aggressive with memory
    cl_float4 *pos_new = calloc(node->count, sizeof(cl_float4));
    cl_float4 *vel_new = calloc(node->count, sizeof(cl_float4));
    uint64_t *mort_new = calloc(node->count, sizeof(uint64_t));
    for (int i = 0; i < node->count; i++)
    {
        int new_ind = spots[mask(node->mortons[i], depth)]++;
        pos_new[new_ind] = node->positions[i];
        vel_new[new_ind] = node->velocities[i];
        mort_new[new_ind] = node->mortons[i]; //we could save some shifts by discarding the left parts of mortons, try later
    }
    memcpy(node->positions, pos_new, node->count * sizeof(cl_float4));
    free(pos_new);
    memcpy(node->velocities, vel_new, node->count * sizeof(cl_float4));
    free(vel_new);
    memcpy(node->mortons, mort_new, node->count * sizeof(uint64_t));
    free(mort_new);

    unsigned int offset = 0;
    node->children = (t_tree **)calloc(8, sizeof(t_tree *));
    for (int i = 0; i < 8; i++)
    {
        node->children[i] = (t_tree *)calloc(1, sizeof(t_tree));
        node->children[i]->positions = &(node->positions[offset]);
        node->children[i]->velocities = &(node->velocities[offset]);
        node->children[i]->mortons = &(node->mortons[offset]);
        node->children[i]->count = counts[i];
        node->children[i]->parent = node;
        node->children[i]->children = NULL;
        node->children[i]->bounds = bounds_from_code(node->bounds, i);
        offset += counts[i];
    }
    free(counts);
    free(spots);
}

void split(t_tree *node)
{
	//split this cell into 8 octants
    node->children = (t_tree **)calloc(8, sizeof(t_tree *));
    int depth = node_depth(node);
    unsigned int offset = 0;
    //printf("this cell has %d bodies\n", node->count);
    for (unsigned int i = 0; i < 8; i++)
    {
        node->children[i] = (t_tree *)calloc(1, sizeof(t_tree));
        node->children[i]->positions = &(node->positions[offset]);
        node->children[i]->velocities = &(node->velocities[offset]);
        node->children[i]->mortons = &(node->mortons[offset]);
        node->children[i]->count = 0;
        node->children[i]->parent = node;
        node->children[i]->children = NULL;
        node->children[i]->bounds = bounds_from_code(node->bounds, i);
		
		//scan through array for borders between 3-bit substring values for this depth
		//these are the dividing lines between octants
        unsigned int j = binary_border_search(node->mortons, offset, node->count, i, depth);
        offset += j;
        node->children[i]->count = j;
        //printf("child %d has %d bodies\n", i, j);
    }
}

void split_tree(t_tree *root)
{
	//we divide the tree until each leaf cell has < leaf threshold bodies.
	//we also have to halt at depth 21, but it is unlikely to need to divide that far.
    if (root->count < LEAF_THRESHOLD || node_depth(root) == 21)
    {
        root->as_single = make_as_single(root);
        return;
    }
    new_split(root);
    for (int i = 0; i < 8; i++)
        split_tree(root->children[i]);
    root->as_single = make_as_single(root);
}

void *split_thread(void *param)
{
    t_tree *root = (t_tree *)param;
    if (root->count < LEAF_THRESHOLD || node_depth(root) == 21)
    {
        root->as_single = make_as_single(root);
        return (0);
    }
    new_split(root);
    for (int i = 0; i < 8; i++)
        split_tree(root->children[i]);
    root->as_single = make_as_single(root);
    return(0);
}

void mt_split_tree(t_tree *root)
{
    pthread_t *split_threads = calloc(8, sizeof(pthread_t));
    new_split(root);
    for (int i = 0; i < 8; i++)
        pthread_create(&split_threads[i], NULL, split_thread, root->children[i]);
    for (int i = 0; i < 8; i++)
        pthread_join(split_threads[i], NULL);
    free(split_threads);
}

t_tree *make_tree(t_dataset *data)
{
	//determine bounding cube of particle set
    t_bounds root_bounds = bounds_from_bodies(data->positions, data->particle_cnt);

    uint64_t *mortons = calloc(data->particle_cnt, sizeof(uint64_t));
    for (int i = 0; i < data->particle_cnt; i++)
        mortons[i] = mortonize(data->positions[i], root_bounds);
    // ^^ could mthread this but meh for now ^^
    t_tree *root = new_tnode(data->positions, data->velocities, data->particle_cnt, NULL);
    root->bounds = root_bounds;
    root->mortons = mortons;
    //recursively divide the tree
    mt_split_tree(root);
    printf("made tree\n");
    free(mortons);
    return (root);
}

static void free_tree(t_tree *t)
{
    if (!t)
        return;
    if (t->children)
        for (int i = 0; i < 8; i++)
            free_tree(t->children[i]);
    free(t->neighbors);
    free(t->children);
    if (t->as_single)
    {
        free(t->as_single->positions);
        free(t->as_single->velocities);
        free(t->as_single);
    }
    free(t);
}

void    divide_dataset(t_standalone *sim)
{
    printf("divide start, particle count %ld\n", sim->dataset->particle_cnt);
    static t_tree *t;

    if (t != NULL)
        free_tree(t);
    t = make_tree(sim->dataset);
    sim->cells = enumerate_better(t);
    sim->cell_count = count_tree_array(sim->cells);
    printf("counted leaves, %d\n", sim->cell_count);
    printf("starting assemble\n");
    mt_assemble_neighborhoods(sim->cells, t, sim);
    printf("divide end (assembled)\n");
}
