#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

// �����С100*100
#define WIDTH 100
#define HEIGHT 100
#define SIZE 10000

// ��ѧ��Ӧ���ʵļ��㷽ʽ
#define CALC_AA(x) 1*10*0.16*0.443*x/SIZE
#define CALC_BA(x) 6*509.36*0.16*pow((0.443*x/SIZE), 3)
#define CALC_AG(y) 0.05*0.443*y/SIZE
#define CALC_BG(z) 0.2*0.443*z/(SIZE*3)
#define CALC_AD(y) 0.0386*0.443*y/SIZE;
#define CALC_BD(z) 0.0657*0.443*z/(SIZE*3)
#define CALC_AB(y,z) 7*0.443*0.443*y*z/(SIZE*SIZE*3)

// �������
int surface[100][100];
// ��������Ӹ���
//extern int C2_NUM, C6_NUM, EMPTY_NUM;
// ���ڱ���˫���ӷ�ӦʱC2������
int double_react_x, double_react_y;

// ��������C6�����������ƥ��ģʽ
int C6_PATTERN[][4] = {
		{1,0,0,-1}, {1,-1,1,0}, {1,-1,0,-1},// ����
		{1,0,0,1}, {1,1,0,1}, {1,1,0,1}, // ����
		{-1,0,0,-1}, {-1,-1,0,-1}, {-1,1,-1,0}, // ����
		{1,-1,0,-1}, {-1,-1,0,-1}, {0,1,-1,0}, // ����
		{0,-1,0,1}, // ��
		{1,0,-1,0} // ��
		};

// ˫���ӷ�Ӧʱ���ڼ��C2��Χ��ûC6
int C2_C6_REACT_PATTERN[][2] = {
	{0,1},{1,0},{-1,0},{0,-1}};

// ���ɷ�Ӧ����
void genSurface(int C2_NUM, int C6_NUM){
	// ��ʼ���������
	memset(surface, 0, sizeof(int) * SIZE);

	int i, x, y, p, t, j;
	// ����C6�ĵ�
	for(i = 0; i < C6_NUM; i+=3){
		x = rand() % (WIDTH);
		y = rand() % (HEIGHT);
		p = rand() % (12);
		t = rand() % (100);

		while(1){
			if(x + C6_PATTERN[p][0] < 0 ||
				x + C6_PATTERN[p][0] >= WIDTH ||
				y + C6_PATTERN[p][1] < 0 ||
				y + C6_PATTERN[p][1] >= HEIGHT ||
				x + C6_PATTERN[p][2] < 0 ||
				x + C6_PATTERN[p][2] >= WIDTH ||
				y + C6_PATTERN[p][3] < 0 ||
				y + C6_PATTERN[p][3] >= HEIGHT ||
				!surface[x][y] ||
				!surface[x + C6_PATTERN[p][0]][y + C6_PATTERN[p][1]] ||
				!surface[x + C6_PATTERN[p][2]][y + C6_PATTERN[p][3]]){
				x = rand() % (WIDTH);
				y = rand() % (HEIGHT);
				p = rand() % (12);
			}else{
				break;
			}
		}
		surface[x + C6_PATTERN[p][2]][y + C6_PATTERN[p][3]] = t;
		surface[x + C6_PATTERN[p][0]][y + C6_PATTERN[p][1]] = t;
		surface[x][y] = t;
	}
	if(C2_NUM != 0){
		t = rand() % (C2_NUM);
	}

	// ����C2�ĵ�
	for(i = 0; i < C2_NUM; i++){
		x = rand() % (WIDTH);
		y = rand() % (HEIGHT);
		while(1){
			if(surface[x][y]){
				x = rand() % (WIDTH);
				y = rand() % (HEIGHT);
			}else{
				break;
			}
		}
		surface[x][y] = -1;
		if(t == i){
			double_react_x = x;
			double_react_y = y;
		}
	}
}

// �ж�C6�Ƿ��ܽ���������Ӧ
int canC6_Adsorb(int C2_NUM, int C6_NUM){
	int i, x, y;
	x = rand() % (WIDTH);
	y = rand() % (HEIGHT);
	genSurface(C2_NUM, C6_NUM);

	for(i = 0; i < 12; i++){
		if( (x + C6_PATTERN[i][0]) >=0 &&
			(x + C6_PATTERN[i][0]) < WIDTH &&
			(x + C6_PATTERN[i][2]) >= 0 &&
			(x + C6_PATTERN[i][2]) < WIDTH &&
			(y + C6_PATTERN[i][1]) >=0 &&
			(y + C6_PATTERN[i][1]) < HEIGHT &&
			(y + C6_PATTERN[i][3]) >= 0 &&
			(y + C6_PATTERN[i][3]) < HEIGHT){
			if(surface[x][y] == 0 &&
					surface[x + C6_PATTERN[i][0]][y + C6_PATTERN[i][1]] == 0 &&
					surface[x + C6_PATTERN[i][2]][y + C6_PATTERN[i][3]] == 0){
				return 1;
			}
		}
	}
	return 0;
}

// �ж�C2��C6�Ƿ��ܷ���˫���ӷ�Ӧ
int canC2_C6_React(int C2_NUM, int C6_NUM){
	int i;

	genSurface(C2_NUM, C6_NUM);
	for(i = 0; i < 4; i++){
		if(surface[double_react_x + C2_C6_REACT_PATTERN[i][0]][double_react_y + C2_C6_REACT_PATTERN[i][1]]){
			return 1;
		}
	}
	return 0;
}


int tag = 0;
void alert(int t, int C2_NUM, int C6_NUM, int EMPTY_NUM){
	if(0 != (C2_NUM + C6_NUM + EMPTY_NUM) && tag == 0){
		cout<<"fuck:"<<t<<endl;
		tag = 1;
		cout<<"C2:"<<C2_NUM<<" C6_NUM:"<<C6_NUM<<" EMPTY_NUM:"<<EMPTY_NUM<<" Total:"<<(C2_NUM + C6_NUM + EMPTY_NUM)<<endl;
	}
}

int main(int argc, char **argv) {
	// ѭ�����
	int temp_x, temp_y;
	// ����ʱ����
	clock_t start, finish, whole_start, whole_end;

	int C2_NUM = 0, C6_NUM = 0, EMPTY_NUM = SIZE;

	// ��Ӧ�����Լ���Ӧ���
	float all, paa, pba, pag, pbg, pad, pbd, pab;
	float AA, BA, AG, BG, AD, BD, AB;
	// ��������
	float r1;
	srand(time(NULL));

	// ��ʼȫ�ּ�ʱ
	whole_start = clock();
	// ��10��ʵ��ÿ����10000�η�Ӧ
	for(temp_x = 0; temp_x < 100; temp_x++){
		// ��ʼ��C2 C6�������
		C2_NUM = 0; C6_NUM = 0; EMPTY_NUM = SIZE;
		// ��ʼ�ֲ���ʱ
		start = clock();

		for(temp_y = 0; temp_y < 20000; temp_y++){
			// ���㵱ǰ��Ӧ����
			AA = CALC_AA(EMPTY_NUM);
			BA = CALC_BA(EMPTY_NUM);
			AG = CALC_AG(C2_NUM);
			BG = CALC_BG(C6_NUM);
			AD = CALC_AD(C2_NUM);
			BD = CALC_BD(C6_NUM);
			AB = CALC_AB(C2_NUM,C6_NUM);
			all=AA+BA+AG+BG+AD+BD+AB;
			paa=AA/all;
			pba=BA/all;
			pag=AG/all;
			pbg=BG/all;
			pad=AD/all;
			pbd=BD/all;
			pab=AB/all;
			// ��������������жϷ����ĸ���Ӧ
			r1 = ((float)rand()/RAND_MAX);
			//cout<<r1<<endl;
			// ������ӦC2����
			if(r1<=paa){
				if(EMPTY_NUM > 0){
					C2_NUM++;
					EMPTY_NUM--;
				}
				//alert(1, C2_NUM, C6_NUM, EMPTY_NUM);
			}
			// ������ӦC6����
			else if(r1<=paa+pba){
				//genSurface(C2_NUM, C6_NUM);
				if(canC6_Adsorb(C2_NUM, C6_NUM)){
					if(EMPTY_NUM > 2){
						C6_NUM = C6_NUM + 3;
						EMPTY_NUM = EMPTY_NUM - 3;
					}
				}
				//alert(2, C2_NUM, C6_NUM, EMPTY_NUM);
			}
			// C2������
			else if(r1<=paa+pba+pag){
				if(C2_NUM > 0){
					C2_NUM--;
					EMPTY_NUM++;
				}
				//alert(3, C2_NUM, C6_NUM, EMPTY_NUM);
			}
			// C2����
			else if(r1<=paa+pba+pag+pad){
				if(C2_NUM > 0){
					C2_NUM--;
					EMPTY_NUM++;
				}
				//alert(4, C2_NUM, C6_NUM, EMPTY_NUM);
			}
			// C6�Ѹ�
			else if(r1<=paa+pba+pag+pad+pbg){
				if(C6_NUM > 2){
					C6_NUM-=3;
					EMPTY_NUM+=3;
				}
				//alert(5, C2_NUM, C6_NUM, EMPTY_NUM);
			}
			// C6����
			else if(r1<=paa+pba+pag+pad+pbg+pbd){
				if(C6_NUM > 2){
					C6_NUM-=3;
					EMPTY_NUM+=3;
				}
				//alert(6, C2_NUM, C6_NUM, EMPTY_NUM);
			}else{// ˫���ӷ�Ӧ
				//genSurface(C2_NUM, C6_NUM);
				if(canC2_C6_React(C2_NUM, C6_NUM)){
					if(C6_NUM > 2 && C2_NUM > 0){
						C2_NUM--;
						C6_NUM-=3;
						EMPTY_NUM+=4;
					}
				}
				//alert(7, C2_NUM, C6_NUM, EMPTY_NUM);
			}
		}
		// �����ֲ���ʱ
		finish = clock();
		// �����Ӧ������Ӹ���
		cout<<"C2:"<<C2_NUM<<" C6_NUM:"<<C6_NUM<<" EMPTY_NUM:"<<EMPTY_NUM<<endl;
		// ����ֲ���ʱ
		cout<<"Using time:"<<(double)(finish - start) / CLOCKS_PER_SEC<<endl;
	}
	// ����ȫ�ּ�ʱ
	whole_end = clock();
	// ���������ʱ
	cout<<endl<<"####################"<<endl<<"Whole using time:"<<(double)(whole_end - whole_start) / CLOCKS_PER_SEC<<endl;

	/*
	������Ҫ������й�
	*/

	cout<<"C2:"<<C2_NUM<<" C6_NUM:"<<C6_NUM<<" EMPTY_NUM:"<<EMPTY_NUM<<endl;

	// ������������ľ���
	genSurface(C2_NUM, C6_NUM);
	int i, j;
	// �����result�ļ�����֮����python��ʾ
	ofstream f ("result.txt");
	for(i = 0; i < WIDTH; i++){
		for(j = 0; j < HEIGHT; j++){
			f<<surface[i][j];
			if(j != (HEIGHT-1)){
				f<<" ";
			}else{
				f<<endl;
			}
		}
	}
	f.close();

	// ����
	getchar();
	return 1;
}
