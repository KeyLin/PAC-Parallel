/*����Ϊ��� R=b/a��a,b��ʱ���������� C2 �� C6 ��ָ�������ʱ��ı仯*/
/*����ͷ�ļ�*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /*ʹ�õ�ǰʱ����Ϊ���������*/
#include <math.h>
#define m 50           /*���������л��е���������������������߳�*/
#define pp 100       /*���ؿ�����ʱ�� MCS����������ʵʱ��*/
#define a 0.5
#define b 0.05
int main()
{  /*�������*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc;
	
	int i,j,n,nn,p,q,n1,n2,n3,s,k,w,o,l;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*ÿ��MCS��ϵͳ���淢���ĸ���Ӧ��*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*ϵͳ�ڸ���Ӧ����*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*ϵͳ�ڸ���Ӧ����*/
	float r1,r3;int r2;                     /*���ɾ��������*/
	int sol[m][m]={0};                   /*��ʼ���������*/
	/*��ÿ MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
	int tor[pp];int ar[pp];int br[pp];int ax[pp];int bx[pp];int at[pp];int bt[pp];
	int small[pp];int large[pp];int empty[pp];     /*�洢ÿһ��ѭ��ǰ�ı���״̬*/
	/*�����ָ�룬���������� TXT ��ʽ���*/
	FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fp11,*fp12,*fp13,*fp14;
	srand(time(NULL));                      /*��ʼ�������*/
	/*��ѭ�� PP �Σ�PP ��ֵΪ���ؿ���ʱ�� MCS*/
	int count=1;
	for(n=0;n<pp;n++)
	{
	    printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*ÿ��ѭ��ǰ������Ӧ������*/
		/*һ�� MCS������ m*m ��ѭ��*/
		for(nn=0;nn<m*m;nn++)
		{
			pc = 0;
			n1=0,n2=0;n3=0,s=0,k=0,l=1;
			/*ɨ�������������棬������Ŀ�λ����Ϊ n1������С������ C2 ����Ϊ n2��С���ӷ����� C6 ����Ϊ n3*/
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				{
					if(sol[i][j]==0){n1=n1+1;}  /*ͳ�ƻ�������λ��*/
					else if(sol[i][j]==-1){n2=n2+1;}/*ͳ�ƻ����������С������ C2 ��*/
					else{n3=n3+1;}              /*ͳ�ƻ������С���ӷ����� C6 ��*/
				}
			}
			/*ȡÿ MCS �ĵ�һ��ѭ����Ϊ��� MCS �ڵı���״̬���*/
			if(nn==0)
			{
				empty[n]=n1;small[n]=n2;large[n]=n3;
			}
			/*������������ C6 �Բ�ͬ����������ʾ����������*/
			/*���� sol �����в����ڵ���С������ l*/
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				{
					if(sol[i][j]==l) {l=l+1;i=0;}
				}
			}

			/*��������ؼ���ʽ*/
			aa=10*a*0.16*0.443*n1/2500;                                          /* C2 ��������*/
			ba=509.36*b*0.16*(0.443*n1/2500)*(0.443*n1/2500)*(0.443*n1/2500);    /* C6 ��������*/
			ag=0.05*0.443*n2/2500;                                               /* C2 �Ѹ�����*/
			bg=0.2*0.443*n3/(2500*3);                                            /* C6 �Ѹ�����*/
			ad=0.0386*0.443*n2/2500;                                             /* C2 ��������*/
			bd=0.0657*0.443*n3/(2500*3);                                         /* C6 ��������*/
			ab=7*0.443*0.443*n2*n3/(2500*2500*3);                                /*˫���ӷ�Ӧ����*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*������*/

			/*��Ӧ��Ӧ����*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*����[0,1]���������*/
			r1=(float)rand()/RAND_MAX;
			r3=(float)rand()/RAND_MAX;

			/*�������Ӧ�����㷨��Ϊ�����߲���ʵ��*/
			/*1����������С������ C2,��-1 ��ʾ*/
			if(r1<=paa)
			{
				if( n1<1 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;
					q=(r2-1)%m;
				}while(sol[p][q]!=0);  /*�������һ�����λ��p,q��*/
				sol[p][q]=-1;nax=nax+1;    /*�ڴ˿�λ������ C2��C2 �������� 1 */
			}

			/*2������С���ӷ����� C6������������ʾ*/
			else if(r1<=paa+pba)
			{
				if( n1<3 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=0);    /*�������һ�����λ*/

				x1=p-1,y1=q;
				x2=p,y2=q-1;
				x3=p+1,y3=q;
				x4=p,y4=q+1;
				if(x1>=0&&x1<m&&y1>=0&&y1<m&&sol[x1][y1]==0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<m&&y2>=0&&y2<m&&sol[x2][y2]==0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<m&&y3>=0&&y3<m&&sol[x3][y3]==0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<m&&y4>=0&&y4<m&&sol[x4][y4]==0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==4)
				{
					sol[p][q]=l;
					sol[coord[0]][coord[1]]=l;
					sol[coord[2]][coord[3]]=l;
					nbx=nbx+1;
				}
				if(pc==6)
				{
					sol[p][q]=l;nbx=nbx+1;
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� l ��ʾ��ͬʱ C6 �������� 1*/
					if(r3<=(float)1/3){sol[coord[0]][coord[1]]=l;sol[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/3){sol[coord[2]][coord[3]]=l;sol[coord[4]][coord[5]]=l;}
					else{sol[coord[0]][coord[1]]=l;sol[coord[4]][coord[5]]=l;}
				}
				if(pc==8)
				{
					sol[p][q]=l;nbx=nbx+1;
					/*�ڴ�����λλ�����ѡ������λ��˵��γ�����λ���� C6��ͬ��������C6 �� l ��ʾ��ͬʱ C6 �������� 1*/
					if(r3<=(float)1/4){sol[coord[0]][coord[1]]=l;sol[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/4){sol[coord[2]][coord[3]]=l;sol[coord[4]][coord[5]]=l;}
					else if(r3<=(float)3/4){sol[coord[4]][coord[5]]=l;sol[coord[6]][coord[7]]=l;}
					else{sol[coord[0]][coord[1]]=l;sol[coord[6]][coord[7]]=l;}
				}
			}

			/*3������С������ C2 �Ѹ�*/
			else if(r1<=paa+pba+pag)
			{
				if( n2<1 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);  /*���ѡ��һ C2 ����λ*/
				sol[p][q]=0;nat=nat+1;  /* C2 �Ѹ���ͬʱ C2 �Ѹ����� 1*/
			}

			/*4������С������ C2 ����*/
			else if(r1<=paa+pba+pag+pad)
			{
				if( n2<1 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);    /*���ѡ��һ C2 ����λ*/
				sol[p][q]=0;nad=nad+1;  /* C2 �Ѹ���ͬʱ C2 �������� 1*/
			}

			/*5��С���ӷ����� C6 �Ѹ�*/
			else if(r1<=paa+pba+pag+pad+pbg)
			{
				if( n3<3 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]<1);         /*���ѡ��һ C6 ����λ*/
				o=sol[p][q];sol[p][q]=0;nbt=nbt+1;    /* C6 �Ѹ����� 1*/
				/* C6 �Ѹ����� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				for(i=0;i<m;i++)
				{
					for(j=0;j<m;j++)
					{
						if(sol[i][j]==o){sol[i][j]=0;}
					}
				}
			}

			/*6��С���ӷ����� C6 ����*/
			else if(r1<=paa+pba+pag+pad+pbg+pbd)
			{
				if( n3<3 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]<1);      /*���ѡ��һ C6 ����λ*/
				o=sol[p][q];nbd=nbd+1;  /* C6 �������� 1*/
				/* C6 �������� C6 ռ�ݵ���������λ�ϵ�ֵ����Ϊ 0*/
				for(i=0;i<m;i++)
				{
					for(j=0;j<m;j++)
					{
					if(sol[i][j]==o){sol[i][j]=0;}
					}
				}
			}

			/*7������С������ C2 ����ΧС���ӷ����� C6 ��Ӧ*/
			else
			{
				if( n2<1 || n3<3 ) continue;
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);    /*���ѡ��һ C2 ����λ*/

				x1=p-1,y1=q;
				x2=p,y2=q-1;
				x3=p+1,y3=q;
				x4=p,y4=q+1;
				if(x1>=0&&x1<m&&y1>=0&&y1<m&&sol[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<m&&y2>=0&&y2<m&&sol[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<m&&y3>=0&&y3<m&&sol[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<m&&y4>=0&&y4<m&&sol[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*�������һ�� C6 ����λ�������䷢��˫���ӷ�Ӧ*/
					sol[p][q]=0;    /*���� C2 ����λת��Ϊ��λ*/
					o=sol[coord[0]][coord[1]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
					nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
				}
				if(pc==4)
				{/*�������λ���� C6 ����*/
					sol[p][q]=0;    /*���� C2 ����λת��Ϊ��λ*/
					/*���ѡ������һ��λ*/
					if(r3<0.5) o=sol[coord[0]][coord[1]];
					else o=sol[coord[2]][coord[3]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
					nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
				}
				if(pc==6)
				{/*�������λ���� C6 ����*/
					sol[p][q]=0;    /*���� C2 ����λת��Ϊ��λ*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					if(r3<=(float)1/3) o=sol[coord[0]][coord[1]];
					else if(r3<=(float)2/3) o=sol[coord[2]][coord[3]];
					else o=sol[coord[4]][coord[5]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
					nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
				}
				if(pc==8)
				{/*�������λ���� C6 ����*/
					sol[p][q]=0;    /*���� C2 ����λת��Ϊ��λ*/
					/*�ڴ�����λλ�����ѡ��һ��λ*/
					if(r3<=(float)1/4) o=sol[coord[0]][coord[1]];
					else if(r3<=(float)2/4) o=sol[coord[2]][coord[3]];
					else if(r3<=(float)3/4) o=sol[coord[4]][coord[5]];
					else o=sol[coord[6]][coord[7]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*���� C6 �漰��������λ����ת��Ϊ��λ*/
					nr=nr+1;  /*����˫���ӷ�Ӧ���� 1*/
				}
			}
		}
/*���� MCS ��ϵͳ����Ӧ���洢�ڸ�������*/
		tor[n]=nr;  //tor �洢˫���ӷ�Ӧ��
		ar[n]=nad;  //ar �洢С���ӳ�����
		br[n]=nbd;  //br �洢�����������
		ax[n]=nax;  //ax �洢С����������
		bx[n]=nbx;  //bx �洢�����������
		at[n]=nat;  //at �洢С�����Ѹ���
		bt[n]=nbt;  //bt �洢������Ѹ���

// /* /*���·ֱ𽫼����м�ʱ�̵ȵı�������״̬�� TXT ��ʽ���*/
// 		if(n==499)
// 		{
// 			fp3=fopen("./500.txt","w");
// 			for(i=0;i<m;i++)
// 			{
// 				for(j=0;j<m;j++)
// 				fprintf(fp3,"%d ",sol[i][j]);
// 				fprintf(fp3,"\n");
// 			}
// 			fclose(fp3);
// 		}
// //�� �� */
	}

/*���½����յı�������״̬�� TXT ��ʽ���*/
	printf("Empty_NUM:%d C2_NUM:%d C6_NUM:%d",n1,n2,n3);
 	fp1=fopen("./20000.txt","w");
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			fprintf(fp1,"%d ",sol[i][j]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);
/*���·ֱ�����ģ�������ÿ�� MCS �ڵĸ���Ӧ���� TXT ��ʽ���*/
	fp2=fopen("./result.txt","w");
	for(i=0;i<pp;i++)
	{
		fprintf(fp2,"%d %d %d %d %d %d %d\n",tor[i],ar[i],br[i],ax[i],bx[i],at[i],bt[i]);
	}
	fclose(fp2);
	return 0;
}
