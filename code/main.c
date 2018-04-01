#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cfd.h"

int main(int argc, char *argv[]){	
	clock_t begin = clock();
	
	clock_t end = clock();
	printf("%lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
}
