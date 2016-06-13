#include <stdio.h>
#include<stdlib.h>
#include <time.h>
#include <math.h>

# define Nt 10//�۲����е�ʱ����
# define Nn 5 //����Ȩ�ر߳�������ĸ���
# define Ns 1
# define WeightDensity 0.4
# define pi 3.1415926
# define r  -5  //f(x)�Ĳ���
FILE *fp1,*fp2, *fp11, *fp22;

double w_lowBound = -1;//Ȩ��ȡֵ��Χ����
double w_uppBound = 1;//Ȩ��ȡֵ��Χ����
double s_lowBound = 0;//״̬ȡֵ��Χ����
double s_uppBound = 1;//״̬ȡֵ��Χ����
//long   rnd_uni_init;//64λ���з�������
double G[Ns][Nn][Nt]={0};//���룬�۲�ֵ
double S[Ns][Nn][Nt]={0};//��Ӧ
double W[Nn][Nn]={0};//
double love; //�����������������love

struct Individual	//����Ⱦɫ�����ṹ��
{ 
	double weight[Nn][Nn];	//�������
	double error ;     //����ʺ���
};

/************************************************��������*********************************************/
void   GenerateInitialState();
double TransFunc(double pop0);

void GenerateInitialState()//��ʼ��state
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
	 
	 /*ÿ��response�ĵ�һ�ж��������*/
	 for(int i = 0; i < Ns; i++)
	 {
		for(int j=0; j < Nn; j++)
		{
			love = (double)rand()/RAND_MAX;
			G[i][j][0] = s_lowBound + love * (s_uppBound - s_lowBound); //s_uppBound��״ֵ̬
			S[i][j][0] = G[i][j][0];
		}
	 }
	 
	 //double temp;
	 for(int s = 0; s < Ns; s++)
	 {
		 for(int t=0; t<(Nt-1); t++)//
		 { 
			 for(int n=0; n<Nn; n++)//������Nn=5��
			 {
				 double temp=0.0;
				 for(int m=0; m<Nn; m++)
				 {
					 temp += W[n][m]*G[s][m][t];
				 }
				 G[s][n][t+1] = TransFunc(temp);//��s����Ӧ������n��t+1ʱ�̵Ĺ���״̬
				 S[s][n][t+1] = G[s][n][t+1];
		     }
		 }
	 }
}

double TransFunc(double pop0)// s����
{
       double divider = exp(r*pop0);
	   double f = 1/(divider+1);
	   return f;
} 

void main()
{
      fp1 = fopen("TimeSequence.txt", "a+");//����һ��ֻд�ļ�
	  fp11 = fopen("TimeSequence.tsv", "a+");
	  fp2 = fopen("Network.txt", "a+");//����һ��ֻд�ļ�
	  fp22 = fopen("Network.tsv", "a+");
	  GenerateInitialState();//��ʼ��state
	  
	  for(int s = 0; s < Ns; s++)
	  {
		  //fprintf(fp1, "/******************response %d******************/", s+1);
		  //fprintf(fp11, "/******************response %d******************/", s+1);
		  //fprintf(fp1, "\n");
		  //fprintf(fp11, "\n");
		  for(int n = 0; n < Nn; n++)//������Nn=5��
		  { 
			  //fprintf(fp1, "gene%d  ", n);
			  for(int t = 0; t < Nt; t++)//ʱ����Nt=11
			  {
				  fprintf(fp1, "%f ",G[s][n][t]);//��f1д���ļ�
				  fprintf(fp11, "%f ",G[s][n][t]);
			  }
			  fprintf(fp1, "\n");//��f1д���ļ�
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