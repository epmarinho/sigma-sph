/**@<sph-treecode_binsort.c>::**/
#include <sph-treecode_defaults.h>
double          key;
index_t         id;
counter_t       l;

void qnn_binsort(size_t n)
{
    descriptor_t    i;
    BINTREE        *root = mkbnode(qnnlist.distance[0], qnnlist.list[0]);

    l = 0;
    for (i = 1; i < n; i++) {
        key = qnnlist.distance[i];
        id = qnnlist.list[i];
        qnn_binsort_put(root);
    }
    qnn_binsort_get(root);
}

BINTREE        *root;
BINTREE        *mkbnode(double key, int id)
{
    BINTREE        *t = malloc(sizeof(BINTREE));

    t->key = key;
    t->id = id;
    return t;
}

void qnn_binsort_put(BINTREE * node)
{
    if (key > node->key) {
        if (node->right == NULL)
            node->right = mkbnode(key, id);
        else
            qnn_binsort_put(node->right);
    } else {
        if (node->left == NULL)
            node->left = mkbnode(key, id);
        else
            qnn_binsort_put(node->left);
    }
}

int             l;
void qnn_binsort_get(BINTREE * t)
{
    if (t->left)
        qnn_binsort_get(t->left);

    qnnlist.distance[l] = t->key;
    qnnlist.list[l] = t->id;
    l++;

    if (t->right)
        qnn_binsort_get(t->right);
}
