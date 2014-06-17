/*����Ϊ��� R=C6_CON/a��a,b��ʱ���������� C2 �� C6 ��ָ�������ʱ��ı仯*/
/*����ͷ�ļ�*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /*ʹ�õ�ǰʱ����Ϊ���������*/
#include <math.h>
#define MAT 50           /*���������л��е���������������������߳�*/
#define MC_Time 20000       /*���ؿ�����ʱ�� MCS����������ʵʱ��*/
#define Label_NUM MAT*MAT/3+1  //C6���ӱ�Ƿ�

int MAT_SIZE = MAT*MAT;    //�����С
double C2_CON = 0.5;          //C2Ũ��
double C6_CON = 7.5;         //C6Ũ��
double C2_NUM,C6_NUM,Empty_NUM;

int main()
{  /*�������*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc,label[Label_NUM] = {0};

	int i,j,n,nn,x,y,s,k,w,cur,tag;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*ÿ��MCS��ϵͳ���淢���ĸ���Ӧ��*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*ϵͳ�ڸ���Ӧ����*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*ϵͳ�ڸ���Ӧ����*/
	float R1;int R2,R3;                     /*���ɾ��������*/
	int surf[MAT][MAT]={0};                   /*��ʼ���������*/
	/*��ÿ MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
	int tor[MC_Time];int ar[MC_Time];int br[MC_Time];int ax[MC_Time];int bx[MC_Time];int at[MC_Time];int bt[MC_Time];
	int C2[MC_Time];int C6[MC_Time];int empty[MC_Time];     /*�洢ÿһ��ѭ��ǰ�ı���״̬*/
	/*�����ָ�룬���������� TXT ��ʽ���*/
	FILE *fp1,*fp2;
	srand(time(NULL));                      /*��ʼ�������*/
	/*��ѭ�� MC_Time �Σ�PP ��ֵΪ���ؿ���ʱ�� MCS*/
	Empty_NUM=MAT_SIZE,C2_NUM=0,C6_NUM=0;
	int count=1;
	for(n=0;n<MC_Time;n++)
	{
	    printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*ÿ��ѭ��ǰ������Ӧ������*/
		empty[n]=Empty_NUM;C2[n]=C2_NUM;C6[n]=C6_NUM;       /*ȡÿ MCS �ĵ�һ��ѭ����Ϊ��� MCS �ڵı���״̬���*/
		/*һ�� MCS������ MAT_SIZE ��ѭ��*/
		for(nn=0;nn<MAT_SIZE;nn++)
		{
			pc = 0;
			s=0,k=0,tag=0;
			/*���� label ������δʹ�õ� tag*/
			for(i=1;i<Label_NUM;i++)
			{
				if(label[i]==0) {tag = i;break;}
			}

			/*��������ؼ���ʽ*/
			aa=10*C2_CON*0.16*0.443*Empty_NUM/MAT_SIZE;                                          /* C2 ��������*/
			ba=509.36*C6_CON*0.16*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE)*(0.443*Empty_NUM/MAT_SIZE);    /* C6 ��������*/
			ag=0.05*0.443*C2_NUM/MAT_SIZE;                                               /* C2 �Ѹ�����*/
			bg=0.2*0.443*C6_NUM/(MAT_SIZE*3);                                            /* C6 �Ѹ�����*/
			ad=0.0386*0.443*C2_NUM/MAT_SIZE;                                             /* C2 ��������*/
			bd=0.0657*0.443*C6_NUM/(MAT_SIZE*3);                                         /* C6 ��������*/
			ab=7*0.443*0.443*C2_NUM*C6_NUM/(MAT_SIZE*MAT_SIZE*3);                     /*˫���ӷ�Ӧ����*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*������*/

			/*��Ӧ��Ӧ����*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*����[0,1]���������*/
			R1=(float)rand()/RAND_MAX;
			//R3=(float)rand()/RAND_MAX;

			/*�������Ӧ�����㷨��Ϊ�����߲���ʵ��*/
			/*1����������С������ C2,��-1 ��ʾ*/
			if(R1<=paa)
			{
				if(Empty_NUM<1) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);  /*�������һ�����λ��p,q��*/
				surf[x][y]=-1;nax=nax+1;    /*�ڴ˿�λ������ C2��C2 �������� 1 */
				C2_NUM += 1;
				Empty_NUM -= 1;
			}

			/*2������С���ӷ����� C6������������ʾ*/
			else if(R1<=paa+pba)
			{
				if(Empty_NUM<3) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=0);    /*�������һ�����λ*/

				x1=x-1,y1=y; //��
				x2=x,y2=y-1; //��
				x3=x+1,y3=y; //��
				x4=x,y4=y+1; //��
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]==0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				int flag = 0;
				if(pc==4)
				{
					if(coord[0]!=coord[2]&&coord[1]!=coord[3])
					{
						surf[coord[0]][coord[1]]=tag;
						surf[coord[2]][coord[3]]=tag;
						flag = 1;
					}
				}
				if(pc==6)
				{
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					// if(R3<=(float)1/3){surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					// else if(R3<=(float)2/3){surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					// else{surf[coord[0]][coord[1]]=tag;surf[coord[4]][coord[5]]=tag;}
					do{
						R3=rand()%3;
						if(R3==0)
						{
							if(coord[0]!=coord[2]&&coord[1]!=coord[3])
							{
								surf[coord[0]][coord[1]]=tag;
								surf[coord[2]][coord[3]]=tag;
								flag = 1;
							}
						}
						else if(R3==1)
						{
							if(coord[2]!=coord[4]&&coord[3]!=coord[5])
							{
								surf[coord[2]][coord[3]]=tag;
								surf[coord[4]][coord[5]]=tag;
								flag = 1;
							}
						}
						else
						{
							if(coord[0]!=coord[4]&&coord[1]!=coord[5])
							{
								surf[coord[0]][coord[1]]=tag;
								surf[coord[4]][coord[5]]=tag;
								flag = 1;
							}
						}
					}while(!flag);    /*�����������ڿ�λ*/
				}
				if(pc==8)
				{
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� tag ��ʾ��ͬʱ C6 �������� 1*/
					R3=rand()%4;
					if(R3==0) {surf[coord[0]][coord[1]]=tag;surf[coord[2]][coord[3]]=tag;}
					else if(R3==1) {surf[coord[2]][coord[3]]=tag;surf[coord[4]][coord[5]]=tag;}
					else if(R3==2) {surf[coord[4]][coord[5]]=tag;surf[coord[6]][coord[7]]=tag;}
					else{surf[coord[0]][coord[1]]=tag;surf[coord[6]][coord[7]]=tag;}
					flag = 1;
				}
				if(pc>=2&&flag)
				{	
					surf[x][y]=tag;nbx=nbx+1;
					C6_NUM += 3;
					Empty_NUM -= 3;
					label[tag] = 1;
					flag = 0;
				}
			}

			/*3������С������ C2 �Ѹ�*/
			else if(R1<=paa+pba+pag)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);  /*���ѡ��һ C2 ����λ*/
				surf[x][y]=0;nat=nat+1;  /* C2 �Ѹ���ͬʱ C2 �Ѹ����� 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*4������С������ C2 ����*/
			else if(R1<=paa+pba+pag+pad)
			{
				if(C2_NUM<1) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*���ѡ��һ C2 ����λ*/
				surf[x][y]=0;nad=nad+1;  /* C2 �Ѹ���ͬʱ C2 �������� 1*/
				C2_NUM -= 1;
				Empty_NUM += 1;
			}

			/*5��С���ӷ����� C6 �Ѹ�*/
			else if(R1<=paa+pba+pag+pad+pbg)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);         /*���ѡ��һ C6 ����λ*/
				cur=surf[x][y];nbt=nbt+1;    /* C6 �Ѹ����� 1*/
				/* C6 �Ѹ����� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur){surf[i][j]=0;}
					}
				}
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*6��С���ӷ����� C6 ����*/
			else if(R1<=paa+pba+pag+pad+pbg+pbd)
			{
				if(C6_NUM<3) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]<1);      /*���ѡ��һ C6 ����λ*/
				cur=surf[x][y];nbd=nbd+1;  /* C6 �������� 1*/
				/* C6 �������� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				
				
				for(i=0;i<MAT;i++)
				{
					for(j=0;j<MAT;j++)
					{
						if(surf[i][j]==cur){surf[i][j]=0;}
					}
				}
				
				C6_NUM -= 3;
				Empty_NUM += 3;
				label[cur] = 0;
			}

			/*7������С������ C2 ����ΧС���ӷ����� C6 ��Ӧ*/
			else
			{
				if(C6_NUM<3||C2_NUM<1) continue;
				do{
					R2=rand()%(MAT_SIZE)+1;
					x=(R2-1)/MAT;
					y=(R2-1)%MAT;
				}while(surf[x][y]!=-1);    /*���ѡ��һ C2 ����λ*/

				x1=x-1,y1=y;
				x2=x,y2=y-1;
				x3=x+1,y3=y;
				x4=x,y4=y+1;
				if(x1>=0&&x1<MAT&&y1>=0&&y1<MAT&&surf[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<MAT&&y2>=0&&y2<MAT&&surf[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<MAT&&y3>=0&&y3<MAT&&surf[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<MAT&&y4>=0&&y4<MAT&&surf[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*�������һ�� C6 ����λ�������䷢��˫���ӷ�Ӧ*/
					cur=surf[coord[0]][coord[1]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
				}
				if(pc==4)
				{/*�������λ���� C6 ����*/
					/*���ѡ������һ��λ*/
					R3=rand()%2;
					if(R3==0) cur=surf[coord[0]][coord[1]];
					else cur=surf[coord[2]][coord[3]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
				}
				if(pc==6)
				{/*�������λ���� C6 ����*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					R3=rand()%3;
					if(R3==0) cur=surf[coord[0]][coord[1]];
					else if(R3==1) cur=surf[coord[2]][coord[3]];
					else cur=surf[coord[4]][coord[5]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
				}
				if(pc==8)
				{/*�������λ���� C6 ����*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					R3=rand()%4;
					if(R3==0) cur=surf[coord[0]][coord[1]];
					else if(R3==1) cur=surf[coord[2]][coord[3]];
					else if(R3==2) cur=surf[coord[4]][coord[5]];
					else cur=surf[coord[6]][coord[7]];
					for(i=0;i<MAT;i++)
					{
						for(j=0;j<MAT;j++)
						{
							if(surf[i][j]==cur){surf[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
				}
				if(pc >= 2)
				{
					surf[x][y]=0;    /*���� C2 ����λת��Ϊ��λ*/
					nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
					C6_NUM -= 3;
					C2_NUM -= 1;
					Empty_NUM += 4;
					label[cur] = 0;
				}
			}
		}
		printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf\n",C2_NUM,C6_NUM,Empty_NUM);

/*���� MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
		tor[n]=nr;  //tor �洢˫���ӷ�Ӧ��
		ar[n]=nad;  //ar �洢С���ӳ�����
		br[n]=nbd;  //br �洢�����������
		ax[n]=nax;  //ax �洢С����������
		bx[n]=nbx;  //bx �洢�����������
		at[n]=nat;  //at �洢С�����Ѹ���
		bt[n]=nbt;  //bt �洢������Ѹ���
	}

/*���½����յı�������״̬�� TXT ��ʽ���*/
 	fp1=fopen("./20000.txt","w");
	for(i=0;i<MAT;i++)
	{
		for(j=0;j<MAT;j++)
		{
			fprintf(fp1,"%d ",surf[i][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);
	printf("C2_NUM:%.lf C6_NUM:%.lf Empty_NUM:%.lf",C2_NUM,C6_NUM,Empty_NUM);
/*���·ֱ�����ģ�������ÿ�� MCS �ڵĸ���Ӧ���� TXT ��ʽ���*/
	fp2=fopen("./result.txt","w");
	for(i=0;i<MC_Time;i++)
	{
		fprintf(fp2,"%d %d %d %d %d %d %d\n",tor[i],ar[i],br[i],ax[i],bx[i],at[i],bt[i]);
	}
	fclose(fp2);
	return 0;
}
