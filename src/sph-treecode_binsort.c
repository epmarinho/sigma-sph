/**@<sph-treecode_binsort.c>::**/
#include <sph-treecode_defaults.h>
double          key;
index_t         id;
counter_t       l;

extern void     bintree_release(BINTREE * t);

void qnn_binsort(size_t n)
{
    fprintf(stderr, "# entering qnn_binsort\n");
    descriptor_t    i;
    BINTREE        *root = mkbnode(qnnlist.distance[0], qnnlist.list[0]);

    l = 0;
    for (i = 1; i < n; i++) {
        key = qnnlist.distance[i];
        id = qnnlist.list[i];
        qnn_binsort_put(root);
    }
    qnn_binsort_get(root);
    bintree_release(root);
}

void qnn_sort_bruteforce(void)
{
    //     fprintf(stderr, "# entering qnn_sort_bruteforce\n");
    for (descriptor_t i = 0; i < qnnlist.nn - 1; i++) {
        for (descriptor_t j = i + 1; j < qnnlist.nn; j++) {
            if (qnnlist.distance[i] > qnnlist.distance[j]) {
                double          tmp = qnnlist.distance[j];

                qnnlist.distance[j] = qnnlist.distance[i];
                qnnlist.distance[i] = tmp;
                index_t         pmt = qnnlist.list[j];

                qnnlist.list[j] = qnnlist.list[i];
                qnnlist.list[i] = pmt;
            }
        }
    }
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

void bintree_release(BINTREE * t)
{
    if (t != NULL) {
        bintree_release(t->left);
        bintree_release(t->right);
        free(t);
    }
}

void qnn_swap(int i, int j)
{
    double tmp = qnnlist.distance[i];
    qnnlist.distance[i] = qnnlist.distance[j];
    qnnlist.distance[j] = tmp;
    int pmt = qnnlist.list[i];
    qnnlist.list[i] = qnnlist.list[j];
    qnnlist.list[j] = pmt;
}

void qnn_quicksort(int first, int last)
{
    int             i,
                    j,
                    pivot/*,
                    pmt;
    double          tmp*/;

    if (first < last) {
        pivot = first;
        i = first;
        j = last;
        while (i < j) {
            while (qnnlist.distance[i] <= qnnlist.distance[pivot] && i < last)
                i++;
            while (qnnlist.distance[j] > qnnlist.distance[pivot])
                j--;
            if (i < j) {
                qnn_swap(i,j);

            }
        }
        qnn_swap(pivot, j);
        qnn_quicksort(first, j - 1);
        qnn_quicksort(j + 1, last);
    }
}

// Quick sort in C

// function to swap elements
// void swap(int *a, int *b) {
//     int t = *a;
//     *a = *b;
//     *b = t;
// }

// function to find the partition position
int qnn_partition(int low, int high) {

    // select the rightmost element as pivot
    double pivot = qnnlist.distance[high];

    // pointer for greater element
    int i = low;

    // traverse each element of the array
    // compare them with the pivot
    for (int j = low; j < high; j++) {
        if (qnnlist.distance[j] <= pivot) {

            // if element smaller than pivot is found
            // swap it with the greater element pointed by i
            // swap element at i with element at j
            qnn_swap(i, j);
            i++;
        }
    }

    // swap the pivot element with the greater element at i
    qnn_swap(i, high);

    // return the partition point
    return i;
}

void qnn_quickSort(int low, int high) {
    if (low < high) {

        // find the pivot element such that
        // elements smaller than pivot are on left of pivot
        // elements greater than pivot are on right of pivot
        int p = qnn_partition(low, high);

        // recursive call on the left of pivot
        qnn_quickSort(low, p - 1);

        // recursive call on the right of pivot
        qnn_quickSort(p + 1, high);
    }
}

// function to print array elements
// void printArray(int array[], int size) {
//     for (int i = 0; i < size; ++i) {
//         printf("%d  ", array[i]);
//     }
//     printf("\n");
// }

// main function
// int main() {
//     int data[] = {8, 7, 2, 1, 0, 9, 6};
//
//     int n = sizeof(data) / sizeof(data[0]);
//
//     printf("Unsorted Array\n");
//     printArray(data, n);
//
//     // perform quicksort on data
//     quickSort(data, 0, n - 1);
//
//     printf("Sorted array in ascending order: \n");
//     printArray(data, n);
// }
