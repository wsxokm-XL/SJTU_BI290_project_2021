# SJTU_BI290_project_2021

SJTU BI290 project of 10 group



project:实现在一组DNA序列中寻找转座子的精确匹配的算法。

contributors: Cheng Liu 、zhixuan-t 、 wsxokm-XL



---

版本：1.1 交互式单序列比对版本



特性：一对一的比对DNA序列



源代码：

```c
//1.1.c
#include<stdio.h>
#include<string.h>
#define MAXN 100000
int rank[MAXN], sa[MAXN], tax[MAXN], tp[MAXN], Height[MAXN];
//sa[i] ：排名为i的后缀的位置
//rank[i]：从第i个位置开始的后缀的排名
//tax[i]：i号元素出现了多少次，用于基数排序
//tp[i]:第二关键字排名为i的后缀的位置,用于基数排序
void GetHeight(int N,int rank[],int sa[],char s[]);
void SuffixSort(int N,char s[]);
int main() {
char s1[MAXN],s2[MAXN];
char s[MAXN];
printf("s1:\n");
scanf("%s", s1);
int lenofs1=strlen(s1);
printf("s2:\n");
scanf("%s", s2);
int lenofs2=strlen(s2);
strcat(s+1,s1);
strcat(s+1,"@");
strcat(s+1,s2);
int N; //N：字符串s的长度
N = strlen(s+1);
//printf("s:\n%s\n",s+1);
SuffixSort(N,s);
GetHeight(N,rank,sa,s);
int i,link,target; //判断是否分别来自s1,s2
for(i=1;i<=N;i++)
if(rank[i]==1) {link=i; break;}
for(i=1;i<=N;i++)
if(rank[i]==rank[link+1]+1) {target=i; break;}
int j=0;
if(Height[rank[target]]==lenofs2) printf("The location of the firstcharacter is %d\n",target);
else {j++;printf("Not found!\n");}
while(j==0&&rank[target]<N)
{
    for(i=1;i<=N;i++)
    if(rank[i]==rank[target]+1) {target=i;break;}
    if(Height[rank[target]]>=lenofs2) printf("The location of the firstcharacter is %d\n",target);
    else j++;
}
return 0;
}
void GetHeight(int N,int rank[],int sa[],char s[]) {
int i, j, k = 0;
printf("Height[rank[i]]: ");
for( i = 1; i <= N; i++) {
if(k) k--;
j = sa[rank[i] - 1];
while(s[i + k] == s[j + k]) k++;
Height[rank[i]] = k;}
for( i = 1; i <= N; i++)
if(rank[i]==1) Height[rank[i]] = -1;
for( i = 1; i <= N; i++)
printf("%d ",Height[rank[i]]);
printf("\n");
}
void SuffixSort(int N,char s[]) {
int i,w,p; //p：排名，用于基数排序
//w:当前倍增长度
int M = 69; //M ：字符集的大小，即排名的个数，用于基数排序。
for ( i = 1; i <= N; i++) {rank[i] = s[i] - '0' + 1, tp[i] = i;}
//以下4行为基数排序
for (i = 0; i <= M; i++) tax[i] = 0;
for (i = 1; i <= N; i++) tax[rank[i]]++;
for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
for ( w = 1, p = 0; p < N; M = p, w <= 1) {
p = 0;
for (i = 1; i <= w; i++) tp[++p] = N - w + i;
for (i = 1; i <= N; i++) if (sa[i] > w) tp[++p] = sa[i] - w;
for (i = 0; i <= M; i++) tax[i] = 0;
for (i = 1; i <= N; i++) tax[rank[i]]++;
for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
memcpy(tp,rank,sizeof(rank));
rank[sa[1]] = p = 1;
for (i = 2; i <= N; i++)
rank[sa[i]] = (tp[sa[i - 1]] == tp[sa[i]] && tp[sa[i - 1] + w] ==tp[sa[i] + w]) ? p : ++p;
}
printf("sa[i]: ");
for ( i = 1; i <= N; i++)
printf("%d ", sa[i]);
printf("\n");
printf("rank[i]: ");
for ( i = 1; i <= N; i++)
printf("%d ", rank[i]);
printf("\n");
}

```



---

版本：1.2 交互式多序列比对版本



特性：提供选择，输入要比对的序列个数，可以同时比对多个序列



源代码：

```c
//1.2.c
#include<stdio.h>
#include<string.h>
#define MAXN 100000
#define MAX 5000
int rank[MAXN], sa[MAXN], tax[MAXN], tp[MAXN], Height[MAXN];
//sa[i] ：排名为i的后缀的位置
//rank[i]：从第i个位置开始的后缀的排名
//tax[i]：i号元素出现了多少次，用于基数排序
//tp[i]:第二关键字排名为i的后缀的位置,用于基数排序
void GetHeight(int N,int rank[],int sa[],char s[]);
void SuffixSort(int N,char s[]);
int main() {
char s[MAXN];//后缀数组字符串
int amount,i; //amount 要比对的序列个数
char ss[MAXN]; int lenof;//主序列（如基因组序列）
char s0[amount][MAX];int lenofs[amount]; //比对的序列

printf("Please input the amount of sequence:");
scanf("%d",&amount);

printf("s:\n");
scanf("%s",ss);
lenof=strlen(ss);
strncat(s+1,ss,lenof);

for(i=0;i<amount;i++){
    printf("s%d:\n",i+1);
    scanf("%s",s0[i]);//printf("s%d:YES\n",i+1);
    lenofs[i]=strlen(s0[i]);
    //printf("s%d:YES\n",i+1);
    int ii;
    for(ii=0;ii<=i;ii++) strcat(s+1,"@");
    strncat(s+1,s0[i],lenofs[i]);//printf("%s\n",s+1);
}
//printf("s:\n%s\n",s+1);
int N; //N：字符串s的长度
N = strlen(s+1);//printf("N=%d\n",N);

SuffixSort(N,s);
GetHeight(N,rank,sa,s);

int link=lenof+1,target; //判断是否分别来自s1,s2
//for(i=1;i<=N;i++)
for(i=0;i<amount;i++)
{
    printf("s%d:\n",i+1);
    target=sa[rank[link+1]+1];//printf("link=%d,target=%d,lenofs=%d,Height[rank[target]]=%d\n",link,target,lenofs[i],Height[rank[target]]);
    int j=0,jj=0;
    //if(Height[rank[target]]==lenofs[i]) printf("The location of the firstcharacter is %d\n",target);
    //else {j++;printf("Not found!\n");}
    while(j==0)
    {
        if(Height[rank[target]]>=lenofs[i]) {if(target<=lenof) {printf("The location of the firstcharacter is %d\n",target);jj++;}}
        else j++;
        target=sa[rank[target]+1];
    }
    if(jj==0) printf("Not found!\n");
    link=link+lenofs[i]+i+2;
}
return 0;
}
void GetHeight(int N,int rank[],int sa[],char s[]) {
int i, j, k = 0;
//printf("Height[rank[i]]: ");
for( i = 1; i <= N; i++) {
if(k) k--;
j = sa[rank[i] - 1];
while(s[i + k] == s[j + k]) k++;
Height[rank[i]] = k;}
for( i = 1; i <= N; i++)
if(rank[i]==1) Height[rank[i]] = -1;
/*for( i = 1; i <= N; i++)
printf("%d ",Height[rank[i]]);
printf("\n");*/
}
void SuffixSort(int N,char s[]) {
int i,w,p; //p：排名，用于基数排序
//w:当前倍增长度
int M = 69; //M ：字符集的大小，即排名的个数，用于基数排序。
for ( i = 1; i <= N; i++) {rank[i] = s[i] - '0' + 1, tp[i] = i;}
//以下4行为基数排序
for (i = 0; i <= M; i++) tax[i] = 0;
for (i = 1; i <= N; i++) tax[rank[i]]++;
for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
for ( w = 1, p = 0; p < N; M = p, w <= 1) {
p = 0;
for (i = 1; i <= w; i++) tp[++p] = N - w + i;
for (i = 1; i <= N; i++) if (sa[i] > w) tp[++p] = sa[i] - w;
for (i = 0; i <= M; i++) tax[i] = 0;
for (i = 1; i <= N; i++) tax[rank[i]]++;
for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
memcpy(tp,rank,sizeof(rank));
rank[sa[1]] = p = 1;
for (i = 2; i <= N; i++)
rank[sa[i]] = (tp[sa[i - 1]] == tp[sa[i]] && tp[sa[i - 1] + w] ==tp[sa[i] + w]) ? p : ++p;
}
/*printf("sa[i]: ");
for ( i = 1; i <= N; i++)
printf("%d ", sa[i]);
printf("\n");
printf("rank[i]: ");
for ( i = 1; i <= N; i++)
printf("%d ", rank[i]);
printf("\n");*/
}
```



---

版本：1.3 命令行输入版本



特性：

1. 提供了两个参数：

   ​	[-a] : 可选参数，使用参数将输出建立的后缀数组和LCP



​			[-n amount]（required） ： 必须携带参数 [-n amount] ，amount为要比对的DNA序列个数



​			示例：`-n 2 seq.fastq seq1.fastq seq2.fastq -a`

​			seq.fastq为基因组序列，seq1.fastq和seq2.fastq为要比对的序列



2. 以fastq文件作为输入，结果输出文件“result.txt”



源代码：

```c
//1.3.c
#include<stdio.h>
#include<string.h>
#include<unistd.h>
#include<stdlib.h>
#define MAXN 100000

void usage(){
	printf("Usage():\t[-a] [-n amount](required)\n");
} 
int rank[MAXN], sa[MAXN], tax[MAXN], tp[MAXN], Height[MAXN];
//sa[i] ：排名为i的后缀的位置
//rank[i]：从第i个位置开始的后缀的排名
//tax[i]：i号元素出现了多少次，用于基数排序
//tp[i]:第二关键字排名为i的后缀的位置,用于基数排序
void GetHeight(int N,int rank[],int sa[],char s[]);
void SuffixSort(int N,char s[]);
int main(int argc,char *argv[]) 
{
	char s[MAXN];//后缀数组字符串
	int amount,i; //amount 要比对的序列个数
	int a=0,b=0;//bool a ,输出形式;bool b,是否传入参数n

	int o;
	const char *optstring="an:";

	while((o=getopt(argc,argv,optstring))!=-1)
	{
		
		switch(o){
			case 'a':
				a=1;
				break;
			case 'n':
				b=1;
				amount=atoi(optarg);
				break;
		}
	}

	if(b==0)
	{
		usage();
		return 0;
	}
	
	char ss[MAXN]; int lenof;//主序列（如基因组序列）
	char s0[amount];int lenofs[amount]; //比对的序列


	FILE *fp;char name[amount+1][50];
	if((fp=fopen(argv[optind],"r"))==NULL)
	{
		printf("error when open file %s",argv[optind]);
	}
	else
	{
		fgets(name[0],50,fp);
		name[0][strlen(name[0])-1]='\0';
		fgets(ss,MAXN,fp);
		lenof=strlen(ss)-1;ss[lenof]='\0';
		strncat(s+1,ss,lenof);
	}
	fclose(fp);


	for(i=0;i<amount;i++){
		if((fp=fopen(argv[optind+i+1],"r"))==NULL)
		{
			printf("error when open file %s",argv[optind]);
		}
		else
		{
			fgets(name[i+1],50,fp);
			name[i+1][strlen(name[i+1])-1]='\0';
			fgets(s0,MAXN,fp);
			lenofs[i]=strlen(s0)-1;s0[lenofs[i]]='\0';
			int ii;
    		for(ii=0;ii<=i;ii++) strcat(s+1,"@");
			strncat(s+1,s0,lenofs[i]);
		}
		fclose(fp);
	}

	int N; //N：字符串s的长度
	N = strlen(s+1);

	SuffixSort(N,s);
	GetHeight(N,rank,sa,s);

	int link=lenof+1,target; //判断是否分别来自s1,s2

	fp=fopen("result.txt","w+");
	fprintf(fp,"Result:\n\n");

	if (a==1)
	{
		fprintf(fp,"sa[i]: ");
		for ( i = 1; i <= N; i++)
		fprintf(fp,"%d ", sa[i]);
		fprintf(fp,"\n");
		fprintf(fp,"rank[i]: ");
		for ( i = 1; i <= N; i++)
		fprintf(fp,"%d ", rank[i]);
		fprintf(fp,"\n");
		fprintf(fp,"Height[rank[i]]: ");
		for( i = 1; i <= N; i++)
		fprintf(fp,"%d ",Height[rank[i]]);
		fprintf(fp,"\n\n");
	}

	for(i=0;i<amount;i++)
	{
    	fprintf(fp,"%s in %s :\n",name[i+1],name[0]);//
    	target=sa[rank[link+1]+1];
    	int j=0,jj=0;
    
    	while(j==0)
    	{
        	if(Height[rank[target]]>=lenofs[i]) {if(target<=lenof) {fprintf(fp,"The location of the firstcharacter is %d\n",target);jj++;}}
        	else j++;
        	target=sa[rank[target]+1];
    	}
    	if(jj==0) printf("Not found!\n");
    	link=link+lenofs[i]+i+2;
    	fprintf(fp,"\n");
	}
	fclose(fp);
return 0;
}
void GetHeight(int N,int rank[],int sa[],char s[]) 
{
	int i, j, k = 0;

	for( i = 1; i <= N; i++) {
	if(k) k--;
	j = sa[rank[i] - 1];
	while(s[i + k] == s[j + k]) k++;
	Height[rank[i]] = k;}
	for( i = 1; i <= N; i++)
	if(rank[i]==1) Height[rank[i]] = -1;
}

void SuffixSort(int N,char s[]) 
{
	int i,w,p; //p：排名，用于基数排序 //w:当前倍增长度
	int M = 69; //M ：字符集的大小，即排名的个数，用于基数排序。
	for ( i = 1; i <= N; i++) {rank[i] = s[i] - '0' + 1, tp[i] = i;}
	//以下4行为基数排序
	for (i = 0; i <= M; i++) tax[i] = 0;
	for (i = 1; i <= N; i++) tax[rank[i]]++;
	for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
	for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
	for ( w = 1, p = 0; p < N; M = p, w <= 1) {
	p = 0;
	for (i = 1; i <= w; i++) tp[++p] = N - w + i;
	for (i = 1; i <= N; i++) if (sa[i] > w) tp[++p] = sa[i] - w;
	for (i = 0; i <= M; i++) tax[i] = 0;
	for (i = 1; i <= N; i++) tax[rank[i]]++;
	for (i = 1; i <= M; i++) tax[i] += tax[i - 1];
	for (i = N; i >= 1; i--) sa[ tax[rank[tp[i]]]-- ] = tp[i];
	memcpy(tp,rank,sizeof(rank));
	rank[sa[1]] = p = 1;
	for (i = 2; i <= N; i++)
	rank[sa[i]] = (tp[sa[i - 1]] == tp[sa[i]] && tp[sa[i - 1] + w] ==tp[sa[i] + w]) ? p : ++p;
	}
}

```

