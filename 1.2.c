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