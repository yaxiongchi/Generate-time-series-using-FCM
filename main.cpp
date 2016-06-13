#include <stdio.h>
#include<stdlib.h>
#include <time.h>
#include <math.h>

# define Nt 10//观测序列的时刻数
# define Nn 5 //定义权重边长，基因的个数
# define Ns 1
# define WeightDensity 0.4
# define pi 3.1415926
# define r  -5  //f(x)的参数
FILE *fp1,*fp2, *fp11, *fp22;

double w_lowBound = -1;//权重取值范围下限
double w_uppBound = 1;//权重取值范围上限
double s_lowBound = 0;//状态取值范围下限
double s_uppBound = 1;//状态取值范围上限
//long   rnd_uni_init;//64位的有符号整形
double G[Ns][Nn][Nt]={0};//输入，观测值
double S[Ns][Nn][Nt]={0};//响应
double W[Nn][Nn]={0};//
double love; //产生的随机数将赋予love

struct Individual	//定义染色体个体结构体
{ 
	double weight[Nn][Nn];	//定义变量
	double error ;     //误差率函数
};

/************************************************函数声明*********************************************/
void   GenerateInitialState();
double TransFunc(double pop0);

void GenerateInitialState()//初始化state
{
	 srand((unsigned)time(NULL));
	
	 // generate the weight matrix
	 for(int k=0; k < (Nn*Nn*WeightDensity); k++)//
	 {
		 love = (double)rand()/RAND_MAX;
		 double weightvalue = w_lowBound + love * (w_uppBound - w_lowBound);
	     if(fabs(weightvalue)<0.05)
			 weightvalue = 0.0;
		 
		 int xx = rand()%Nn;
		 int yy = rand()%Nn;
		 if(W[xx][yy] == 0)
		 {
		    W[xx][yy] = weightvalue;
		 }
		 else
			 k--;
	 }
	 
	 /*每个response的第一列都是随机数*/
	 for(int i = 0; i < Ns; i++)
	 {
		for(int j=0; j < Nn; j++)
		{
			love = (double)rand()/RAND_MAX;
			G[i][j][0] = s_lowBound + love * (s_uppBound - s_lowBound); //s_uppBound是状态值
			S[i][j][0] = G[i][j][0];
		}
	 }
	 
	 //double temp;
	 for(int s = 0; s < Ns; s++)
	 {
		 for(int t=0; t<(Nt-1); t++)//
		 { 
			 for(int n=0; n<Nn; n++)//基因有Nn=5个
			 {
				 double temp=0.0;
				 for(int m=0; m<Nn; m++)
				 {
					 temp += W[n][m]*G[s][m][t];
				 }
				 G[s][n][t+1] = TransFunc(temp);//第s个响应，基因n，t+1时刻的估计状态
				 S[s][n][t+1] = G[s][n][t+1];
		     }
		 }
	 }
}

double TransFunc(double pop0)// s函数
{
       double divider = exp(r*pop0);
	   double f = 1/(divider+1);
	   return f;
} 

void main()
{
      fp1 = fopen("TimeSequence.txt", "a+");//建立一个只写文件
	  fp11 = fopen("TimeSequence.tsv", "a+");
	  fp2 = fopen("Network.txt", "a+");//建立一个只写文件
	  fp22 = fopen("Network.tsv", "a+");
	  GenerateInitialState();//初始化state
	  
	  for(int s = 0; s < Ns; s++)
	  {
		  //fprintf(fp1, "/******************response %d******************/", s+1);
		  //fprintf(fp11, "/******************response %d******************/", s+1);
		  //fprintf(fp1, "\n");
		  //fprintf(fp11, "\n");
		  for(int n = 0; n < Nn; n++)//基因有Nn=5个
		  { 
			  //fprintf(fp1, "gene%d  ", n);
			  for(int t = 0; t < Nt; t++)//时刻有Nt=11
			  {
				  fprintf(fp1, "%f ",G[s][n][t]);//把f1写入文件
				  fprintf(fp11, "%f ",G[s][n][t]);
			  }
			  fprintf(fp1, "\n");//把f1写入文件
			  fprintf(fp11, "\n");
		  }
		  fprintf(fp1, "\n");
		  fprintf(fp11, "\n");
	  }

      for(int m = 0; m < Nn; m++)
	  {
		 //fprintf(fp2, "gene%d  ", m);
		 for(int n = 0; n < Nn; n++)
		 {
		     fprintf(fp2, "%f ",W[m][n]);//
			 fprintf(fp22, "%f ",W[m][n]);
		 }
		 fprintf(fp2, "\n");// 
		 fprintf(fp22, "\n");
	  }
      fclose(fp1);
	  fclose(fp11);
	  fclose(fp2);
	  fclose(fp22);
}