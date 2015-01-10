#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int a,b,c;
double d;

int product(int x, int y);

int main(void)
{
 int count=1;
 srand48(time (0));
 while (count<=100000)
 {
    d = drand48();
    printf ("%.20f\n", d);
    count++;
 }
 return 0;
}
