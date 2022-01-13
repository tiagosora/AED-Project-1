/*  Merkle-Hellman Cryptosystem using C to solve subset sums.
    Copyright (C) <2021>  <Tiago Sora>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/

    Contact the author by email : tiagogcarvalho@ua.pt.
*/

// AED, November 2021.
// Solution of the first practical assignement (subset sum problem).
// Place your student numbers and names here.
// The following authors have the same rights upon the code.
// Ana Raquel Paradinha   102491 
// Paulo Pinto            103234
// Tiago Carvalho         104142

#if __STDC_VERSION__ < 199901L
# error "This code must must be compiled in c99 mode or later (-std=c99)" // to handle the unsigned long long data type
#endif
#ifndef STUDENT_H_FILE
# define STUDENT_H_FILE "104142.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include "elapsed_time.h"
#include <math.h>
#include STUDENT_H_FILE

int ctr = 0;
int brute_force(int n, integer_t *p, integer_t sum, int next_index, integer_t partial_sum, int *b){
    if(sum == partial_sum){
        for(int i = next_index; i < n; i++){
            b[i] = 0;
        }
        return 1;
    }
    if(next_index == n){return 0;}
    b[next_index] = 0;
    if(brute_force(n, p, sum, next_index + 1, partial_sum, b) == 1){return 1;}
    b[next_index] = 1; 
    return brute_force(n, p, sum, next_index + 1, partial_sum + p[next_index], b);
}
int branch_and_bound(int n, integer_t *p, integer_t sum, int next_index, integer_t partial_sum, int *b){
    if(sum == partial_sum){
        for(int i = next_index; i < n; i++){
            b[i] = 0;
        }   
        return 1;
    }
    if(next_index == n) return 0;
    if(partial_sum > sum) return 0; 
    else {
        integer_t remaining_sum = 0;
        for(int j = next_index; j < n; j++){
            remaining_sum += p[j];
        }
        if((partial_sum + remaining_sum) < sum){
            return 0;
        } else {
            b[next_index] = 0;
            int result = branch_and_bound(n, p, sum, next_index+1, partial_sum, b);
            if(result == 1){return 1;}
            b[next_index] = 1;
            return branch_and_bound(n, p, sum,next_index+1, partial_sum+p[next_index], b);
        }
    }
}
int meet_in_the_middle(int na, int nb, integer_t *a, integer_t *b, integer_t sum, int *index) {
    int i = 0;
    int j = nb-1;
    integer_t sum_of_sums;
    do {
        sum_of_sums = a[i] + b[j];
        if (((i == na)) | ((j < 0))) {
            return EXIT_FAILURE;
        }
        if (sum_of_sums < sum) {
            i++;
        }
        else if (sum_of_sums > sum) {
            j--;
        }
    } while(sum_of_sums != sum);
    index[0] = i;
    index[1] = j;
    return EXIT_SUCCESS;
}
void subsetSums(integer_t *p, int index, int n_1, integer_t sum, integer_t *subsums) {
    if (index > n_1) {
        subsums[ctr] = sum;
        ctr++;
        return;
    }
    subsetSums(p, index + 1, n_1, sum + p[index],subsums);
    subsetSums(p, index + 1, n_1, sum,subsums);
}
int compare(const void *s1, const void *s2){
    if (*(integer_t*)s1 < *(integer_t*)s2 ) return -1;
    if (*(integer_t*)s1 > *(integer_t*)s2 ) return 1;
    return 0;
}
int main(void){
    int choice;
    // ----- THE USER MAY COMMENT OR UNCOMMENT THE FOLLOWING LINES DEPENDING ON THE INTENDED METHOD ----- //
    // choice = 1; // Use the "brute force" method
    // choice = 2; // Use the "branch and bound" method
    choice = 3; // Use the "meet in the middle" method
    // -------------------------------------------------------------------------------------------------- //
    FILE *stderr;
    FILE *times_file;
    stderr = fopen("Solution_104142.txt", "a"); 
    if (choice == 1){times_file = fopen("Time_104142.txt", "a");}
    if (choice == 2){times_file = fopen("Time2_104142.txt", "a");}
    if (choice == 3){times_file = fopen("Time3_104142.txt", "a");}
    fprintf(stderr,"Program configuration:\n");
    fprintf(stderr,"  min_n ....... %d\n",min_n);
    fprintf(stderr,"  max_n ....... %d\n",max_n);
    fprintf(stderr,"  n_sums ...... %d\n",n_sums);
    fprintf(stderr,"  n_problems .. %d\n",n_problems);
    fprintf(stderr,"  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    for(int i = 0; i < n_problems; i++) {
        int n = all_subset_sum_problems[i].n;
        integer_t *p = all_subset_sum_problems[i].p;
        fprintf(stderr, "n -> %d\n",n);
        double tempo = cpu_time();
        if (choice == 1){ 
            for(int j = 0 ; j < n_sums; j++){
                int *b = (int *) malloc(n * sizeof(int));
                integer_t sum = all_subset_sum_problems[i].sums[j];
                brute_force(n, p, sum, 0, 0, b);
                for(int index = 0; index < n; index++){
                    fprintf(stderr, "%d", b[index]);
                }
                fprintf(stderr, "\n");
                free(b);
            }
        }
        else if (choice == 2){ 
            for(int j = 0; j < n_sums; j++){
                int *b = (int *) malloc(n * sizeof(int));
                integer_t sum = all_subset_sum_problems[i].sums[j];
                branch_and_bound(n, p, sum, 0, 0, b);
                for(int index = 0; index < n; index++){
                    fprintf(stderr, "%d", b[index]);
                }
                fprintf(stderr, "\n");
                free(b);
            }
        }
        else if (choice == 3){ 
            int n1 = n/2;
            int n2 = n - n1;
            integer_t *p1 = (integer_t *) malloc(n1 * sizeof(integer_t));
            integer_t *p2 = (integer_t *) malloc(n2 * sizeof(integer_t));
            for (int k = 0; k < n1; k++){
                p1[k] = p[k];
                p2[k] = p[n1 + k];
            }
            p2[n2-1]=p[n-1];
            int na = pow(2,n1);
            integer_t *a = (integer_t *) malloc(na * sizeof(integer_t));
            subsetSums(p1, 0, n1-1, 0, a);
            ctr=0;
            int nb = pow(2,n2);
            integer_t *b = (integer_t *) malloc(nb * sizeof(integer_t));
            subsetSums(p2, 0, n2-1, 0, b);
            ctr=0;
            qsort(a, na, sizeof(integer_t), compare);
            qsort(b, nb, sizeof(integer_t), compare);
            for(int j = 0; j < n_sums; j++){
                integer_t sum = all_subset_sum_problems[i].sums[j];
                int index[2] = {0,0};
                meet_in_the_middle(na, nb, a, b, sum, index);
                int *bA = (int *) malloc(n1 * sizeof(int));
                int *bB = (int *) malloc(n2 * sizeof(int));
                branch_and_bound(n1, p1, a[index[0]], 0, 0, bA);
                branch_and_bound(n2, p2, b[index[1]], 0, 0, bB);
                int *bin = (int *) malloc(n * sizeof(int));
                for (int k = 0; k < n1; k++){
                    bin[k] = bA[k];
                    bin[n1+k] = bB[k];
                }
                bin[n-1] = bB[n2-1];
                for(int index = 0; index < n; index++){
                    fprintf(stderr, "%d", bin[index]);
                }
                fprintf(stderr, "\n");
                free(bA);
                free(bB);
                free(bin);
            }
            free(p1);
            free(p2);
            free(a);
            free(b);
        }
        tempo = cpu_time() - tempo;
        fprintf(times_file, "%15.17f\n", tempo);
        printf("Method %d : %d done in -> %17.15f seconds \n", choice, n, tempo);
    }
    return EXIT_SUCCESS;
}
