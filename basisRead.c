#include<stdio.h>
#include<stdlib.h>

int main() {
  FILE* fp;
  int buffLen = 100;
  char buff[buffLen], atomSymbol[2], L[2], chArr[10][2];
  double ex, co;
  int orbLen, orbCount = 0; 
  int readFlag = 0, readEx = 0, readCo = 0;

  fp = fopen("H-STO-3G.txt", "r");

  //while(fgets(buff, buffLength, fp)) {
  //    printf("%s\n", buff);
  //    printf("nextline\n");
  //}
  
  // read atom symbol
  fscanf(fp, "%2s", atomSymbol);

  // next char is 0
  fscanf(fp, "%2s", buff);

  // read orbital heading
  readFlag = 1;


  while (fscanf(fp, " %99s", buff) == 1) {
    if (!strcmp("0",buff)) {
      readFlag = 1;
    }
    else if(!strcmp("1.00",buff)) {
      printf("found 1.00\n");
    }
    puts(buff);
  }
  fclose(fp);
}
