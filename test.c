#include <stdio.h>

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n)) 


int main() {
    int n = 27;
    int p = 8;

    printf("n = %d, p = %d:\n", n, p);
    for(int i = 0; i < p; i++) {
        printf("[%d, %d] (%d)\n", BLOCK_LOW(i, p, n), BLOCK_HIGH(i, p, n), BLOCK_SIZE(i, p, n));
    }

    /*
    for(int i = 0; i < 11; i++) {
        int peer = (i&(-1)) + 1 - 2*(i%2);
        printf("%d talks to %d\n", i, peer);
    }
    */

}
