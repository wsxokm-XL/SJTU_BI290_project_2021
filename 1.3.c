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
