/*以下为求解 R=b/a（a,b）时，表面吸附 C2 和 C6 组分覆盖率随时间的变化*/
/*定义头文件*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /*使用当前时间作为随机数种子*/
#include <math.h>
#define m 100           /*方形网格行或列的网格数，代表吸附表面边长*/
#define pp 10000       /*蒙特卡罗总时间 MCS，正比于真实时间*/
#define a 1
#define b 6
int main()
{  /*定义变量*/
	int x1,y1,x2,y2,x3,y3,x4,y4,coord[8],pc;
	
	int i,j,n,nn,p,q,n1,n2,n3,s,k,w,o,l;
	int nr,nad,nbd,nax,nbx,nat,nbt;           /*每个MCS内系统表面发生的各反应数*/
	float aa,ba,ag,bg,ad,bd,ab,all;            /*系统内各反应速率*/
	float paa,pba,pag,pbg,pbd,pad,pab;      /*系统内各反应概率*/
	float r1,r3;int r2;                     /*生成均匀随机数*/
	int sol[m][m]={0};                   /*初始化基体表面*/
	/*将每 MCS 内系统各反应数存储在各数组内*/
	int tor[pp];int ar[pp];int br[pp];int ax[pp];int bx[pp];int at[pp];int bt[pp];
	int small[pp];int large[pp];int empty[pp];     /*存储每一次循环前的表面状态*/
	/*定义各指针，将计算结果以 TXT 格式输出*/
	FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fp11,*fp12,*fp13,*fp14;
	srand(time(NULL));                      /*初始化随机数*/
	/*总循环 PP 次，PP 的值为蒙特卡罗时间 MCS*/
	int count=1;
	for(n=0;n<pp;n++)
	{
	    printf("runing...%d\n",count);count++;

		nr=0;nad=0;nbd=0;nax=0;nbx=0;nat=0;nbt=0;    /*每次循环前将各反应数归零*/
		/*一个 MCS，包含 m*m 次循环*/
		for(nn=0;nn<m*m;nn++)
		{
			pc = 0;
			n1=0,n2=0;n3=0,s=0,k=0,l=1;
			/*扫描整个吸附表面，将表面的空位数记为 n1，线性小分子烃 C2 数记为 n2，小分子芳香烃 C6 数记为 n3*/
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				{
					if(sol[i][j]==0){n1=n1+1;}  /*统计基体表面空位数*/
					else if(sol[i][j]==-1){n2=n2+1;}/*统计基体表面线性小分子烃 C2 数*/
					else{n3=n3+1;}              /*统计基体表面小分子芳香烃 C6 数*/
				}
			}
			/*取每 MCS 的第一次循环作为这次 MCS 内的表面状态输出*/
			if(nn==0)
			{
				empty[n]=n1;small[n]=n2;large[n]=n3;
			}
			/*将表面吸附的 C6 以不同的正整数表示并加以区别*/
			/*查找 sol 数组中不存在的最小正整数 l*/
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				{
					if(sol[i][j]>=l) l=l+1;
				}
			}

			/*以下是相关计算式*/
			aa=10*a*0.16*0.443*n1/2500;                                          /* C2 吸附速率*/
			ba=509.36*b*0.16*(0.443*n1/2500)*(0.443*n1/2500)*(0.443*n1/2500);    /* C6 吸附速率*/
			ag=0.05*0.443*n2/2500;                                               /* C2 脱附速率*/
			bg=0.2*0.443*n3/(2500*3);                                            /* C6 脱附速率*/
			ad=0.0386*0.443*n2/2500;                                             /* C2 沉积速率*/
			bd=0.0657*0.443*n3/(2500*3);                                         /* C6 沉积速率*/
			ab=7*0.443*0.443*n2*n3/(2500*2500*3);                                /*双分子反应速率*/
			all=aa+ba+ag+bg+ad+bd+ab;                                            /*总速率*/

			/*对应反应概率*/
			paa=aa/all;pba=ba/all;pag=ag/all;pbg=bg/all;pad=ad/all;pbd=bd/all;pab=ab/all;

			/*产生[0,1]均匀随机数*/
			r1=(float)rand()/RAND_MAX;
			r3=(float)rand()/RAND_MAX;

			/*具体各方应过程算法分为以下七步来实现*/
			/*1、吸附线性小分子烃 C2,以-1 表示*/
			if(r1<=paa)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;
					q=(r2-1)%m;
				}while(sol[p][q]!=0);  /*随机生成一表面空位（p,q）*/
				sol[p][q]=-1;nax=nax+1;    /*在此空位上吸附 C2，C2 吸附数加 1 */
			}

			/*2、吸附小分子芳香烃 C6，以正整数表示*/
			else if(r1<=paa+pba)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=0);    /*随机生成一表面空位*/

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
					/*在此三邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 l 表示，同时 C6 吸附数加 1*/
					if(r3<=(float)1/3){sol[coord[0]][coord[1]]=l;sol[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/3){sol[coord[2]][coord[3]]=l;sol[coord[4]][coord[5]]=l;}
					else{sol[coord[0]][coord[1]]=l;sol[coord[4]][coord[5]]=l;}
				}
				if(pc==8)
				{
					sol[p][q]=l;nbx=nbx+1;
					/*在此三邻位位上随机选择两邻位与此点形成三空位吸附 C6，同样吸附的C6 以 l 表示，同时 C6 吸附数加 1*/
					if(r3<=(float)1/4){sol[coord[0]][coord[1]]=l;sol[coord[2]][coord[3]]=l;}
					else if(r3<=(float)2/4){sol[coord[2]][coord[3]]=l;sol[coord[4]][coord[5]]=l;}
					else if(r3<=(float)3/4){sol[coord[4]][coord[5]]=l;sol[coord[6]][coord[7]]=l;}
					else{sol[coord[0]][coord[1]]=l;sol[coord[6]][coord[7]]=l;}
				}
			}

			/*3、线性小分子烃 C2 脱附*/
			else if(r1<=paa+pba+pag)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);  /*随机选择一 C2 吸附位*/
				sol[p][q]=0;nat=nat+1;  /* C2 脱附，同时 C2 脱附数加 1*/
			}

			/*4、线性小分子烃 C2 沉积*/
			else if(r1<=paa+pba+pag+pad)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);    /*随机选择一 C2 吸附位*/
				sol[p][q]=0;nad=nad+1;  /* C2 脱附，同时 C2 沉积数加 1*/
			}

			/*5、小分子芳香烃 C6 脱附*/
			else if(r1<=paa+pba+pag+pad+pbg)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]<1);         /*随机选择一 C6 吸附位*/
				o=sol[p][q];sol[p][q]=0;nbt=nbt+1;    /* C6 脱附数加 1*/
				/* C6 脱附：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<m;i++)
				{
					for(j=0;j<m;j++)
					{
						if(sol[i][j]==o){sol[i][j]=0;}
					}
				}
			}

			/*6、小分子芳香烃 C6 沉积*/
			else if(r1<=paa+pba+pag+pad+pbg+pbd)
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]<1);      /*随机选择一 C6 吸附位*/
				o=sol[p][q];nbd=nbd+1;  /* C6 沉积数加 1*/
				/* C6 沉积：将 C6 占据的三个吸附位上的值更新为 0*/
				for(i=0;i<m;i++)
				{
					for(j=0;j<m;j++)
					{
					if(sol[i][j]==o){sol[i][j]=0;}
					}
				}
			}

			/*7、线性小分子烃 C2 与周围小分子芳香烃 C6 反应*/
			else
			{
				do{
					r2=rand()%(m*m)+1;
					p=(r2-1)/m;q=(r2-1)%m;
				}while(sol[p][q]!=-1);    /*随机选择一 C2 吸附位*/

				x1=p-1,y1=q;
				x2=p,y2=q-1;
				x3=p+1,y3=q;
				x4=p,y4=q+1;
				if(x1>=0&&x1<m&&y1>=0&&y1<m&&sol[x1][y1]>0){coord[pc] = x1,coord[pc+1] = y1;pc += 2;}
				if(x2>=0&&x2<m&&y2>=0&&y2<m&&sol[x2][y2]>0){coord[pc] = x2,coord[pc+1] = y2;pc += 2;}
				if(x3>=0&&x3<m&&y3>=0&&y3<m&&sol[x3][y3]>0){coord[pc] = x3,coord[pc+1] = y3;pc += 2;}
				if(x4>=0&&x4<m&&y4>=0&&y4<m&&sol[x4][y4]>0){coord[pc] = x4,coord[pc+1] = y4;pc += 2;}

				if(pc==2)
				{/*如果存在一个 C6 吸附位，则与其发生双分子反应*/
					sol[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					o=sol[coord[0]][coord[1]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==4)
				{/*如果两邻位都被 C6 吸附*/
					sol[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*随机选择其中一邻位*/
					if(r3<0.5) o=sol[coord[0]][coord[1]];
					else o=sol[coord[2]][coord[3]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==6)
				{/*如果三邻位都被 C6 吸附*/
					sol[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*在此三邻位位上随机选择一邻位*/
					if(r3<=(float)1/3) o=sol[coord[0]][coord[1]];
					else if(r3<=(float)2/3) o=sol[coord[2]][coord[3]];
					else o=sol[coord[4]][coord[5]];
					for(i=0;i<m;i++)
					{
						for(j=0;j<m;j++)
						{
							if(sol[i][j]==o){sol[i][j]=0;}
						}
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
				if(pc==8)
				{/*如果四邻位都被 C6 吸附*/
					sol[p][q]=0;    /*将此 C2 吸附位转变为空位*/
					/*在此四邻位位上随机选择一邻位*/
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
					}    /*将此 C6 涉及的三吸附位重新转变为空位*/
					nr=nr+1;  /*表面双分子反应数加 1*/
				}
			}
		}
/*将此 MCS 内系统各反应数存储在各数组内*/
		tor[n]=nr;  //tor 存储双分子反应数
		ar[n]=nad;  //ar 存储小分子沉积数
		br[n]=nbd;  //br 存储大分子吸附数
		ax[n]=nax;  //ax 存储小分子吸附数
		bx[n]=nbx;  //bx 存储大分子吸附数
		at[n]=nat;  //at 存储小分子脱附数
		bt[n]=nbt;  //bt 存储大分子脱附数

/* /*以下分别将计算中间时刻等的表面吸附状态以 TXT 格式输出*/
		if(n==499)
		{
			fp3=fopen("./500.txt","w");
			for(i=0;i<m;i++)
			{
				for(j=0;j<m;j++)
				fprintf(fp3,"%d ",sol[i][j]);
				fprintf(fp3,"\n");
			}
			fclose(fp3);
		}
//… … */
	}

/*以下将最终的表面吸附状态以 TXT 格式输出*/
	printf("Empty:%d C2:%d C6:%",empty[pp-1],small[pp-1],large[pp-1]);
/* 	//fp1=fopen("./20000.txt","w");
	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			if(sol[i][j]==-1) 
		//fprintf(fp1,"%d ",sol[i][j]);
		//fprintf(fp1,"\n");
		}
	}
	//fclose(fp1); */
/*以下分别将整个模拟过程中每个 MCS 内的各反应数以 TXT 格式输出*/
	fp2=fopen("./result.txt","w");
	for(i=0;i<pp;i++)
	{
		fprintf(fp2,"%d %d %d %d %d %d %d\n",tor[i],ar[i],br[i],ax[i],bx[i],at[i],bt[i]);
	}
	fclose(fp2);
	return 0;
}
