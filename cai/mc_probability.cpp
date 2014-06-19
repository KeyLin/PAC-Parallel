#include <iostream>
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
using namespace std;

// 矩阵大小100*100
#define WIDTH 50
#define HEIGHT 50
//#define C2_CON 1 
//#define C6_CON 6
#define MCTime 5000
#define BOARD 2*(WIDTH+HEIGHT)-4
double  SIZE = WIDTH*HEIGHT;

// 化学反应速率的计算方式
#define CALC_AA(x,y) y*10*0.16*0.443*x/SIZE
#define CALC_BA(x,y) y*509.36*0.16*pow((0.443*x/SIZE), 3)
#define CALC_AG(y) 0.05*0.443*y/SIZE
#define CALC_BG(z) 0.2*0.443*z/(SIZE*3)
#define CALC_AD(y) 0.0386*0.443*y/SIZE;
#define CALC_BD(z) 0.0657*0.443*z/(SIZE*3)
#define CALC_AB(y,z) 7*0.443*0.443*y*z/(SIZE*SIZE*3)

// 表面矩阵
int surface[WIDTH][HEIGHT];
// 表面的粒子个数
//extern int C2_NUM, C6_NUM, EMPTY_NUM;
// 用于保存双分子反应时C2的坐标
int double_react_x, double_react_y;

// 用来生成C6另外两个点的匹配模式
int C6_PATTERN[][4] = {
        {1,0,0,-1}, {1,-1,1,0}, {1,-1,0,-1},// 右上
        {1,0,0,1}, {1,1,0,1}, {1,1,0,1}, // 右下
        {-1,0,0,-1}, {-1,-1,0,-1}, {-1,1,-1,0}, // 左下
        {1,-1,0,-1}, {-1,-1,0,-1}, {0,1,-1,0}, // 左上
        {0,-1,0,1}, // 竖
        {1,0,-1,0} // 横
        };

// 双分子反应时用于检测C2周围有没C6
int C2_C6_REACT_PATTERN[][2] = {
    {0,1},{1,0},{-1,0},{0,-1}};

// 生成反应表面
void genSurface(int C2_NUM, int C6_NUM){
    // 初始化表面矩阵
    memset(surface, 0, sizeof(int) * SIZE);
    
    int i, x, y, p, t, j;
    // 生成C6的点
    for(i = 0; i < C6_NUM; i+=3){
        x = rand() % (WIDTH);
        y = rand() % (HEIGHT);
        t = rand() % (100);

        for(p = 0; p < 13; p++){
            if(x + C6_PATTERN[p][0] < 0 ||
                x + C6_PATTERN[p][0] >= WIDTH ||
                y + C6_PATTERN[p][1] < 0 ||
                y + C6_PATTERN[p][1] >= HEIGHT ||
                x + C6_PATTERN[p][2] < 0 ||
                x + C6_PATTERN[p][2] >= WIDTH ||
                y + C6_PATTERN[p][3] < 0 ||
                y + C6_PATTERN[p][3] >= HEIGHT){
                if(!surface[x][y] &&
                !surface[x + C6_PATTERN[p][0]][y + C6_PATTERN[p][1]] &&
                !surface[x + C6_PATTERN[p][2]][y + C6_PATTERN[p][3]]){
                    if(p == 12){
                        i -= 3;
                    }
                }
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

    // 生成C2的点
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

// 判断C6是否能进行吸附反应
int canC6_Adsorb(int C2_NUM, int C6_NUM){
    int i, x, y;
    x = rand() % (WIDTH);
    y = rand() % (HEIGHT);
    genSurface(C2_NUM, C6_NUM);
    
    for(i = 0; i < 13; i++){
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

// 判断C2和C6是否能发生双分子反应
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

bool C6_Adsorb(int EMPTY_NUM)
{
    double p0,p1,p2,p3,p4;
    int r1 = rand()/(int)(SIZE);
    double r2 = ((double)rand()/RAND_MAX);
    double p = (double)EMPTY_NUM/SIZE;
    if(r1<4)
    {
        if(EMPTY_NUM >= 3)
        {
            p0 = (1-p)*(1-p);
            p1 = 2*p*(1-p);
            p2 = p*p;
        }
    }
    else if(r1 < BOARD)
    {
        if(EMPTY_NUM == 3)
        {
            p0 = (1-p)*(1-p);
            p1 = 2*p*(1-p);
            p2 = p*p;
        }
        else if(EMPTY_NUM >= 4)
        {
            p0 = (1-p)*(1-p)*(1-p);
            p1 = 3*p*(1-p)*(1-p);
            p2 = 3*p*p*(1-p);
            p3 = p*p*p;
        }
    }
    else
    {
        if(EMPTY_NUM == 3)
        {
            p0 = (1-p)*(1-p);
            p1 = 2*p*(1-p);
            p2 = p*p;
        }
        else if(EMPTY_NUM == 4)
        {
            p0 = (1-p)*(1-p)*(1-p);
            p1 = 3*p*(1-p)*(1-p);
            p2 = 3*p*p*(1-p);
            p3 = p*p*p;
        }
        else if(EMPTY_NUM >=5 )
        {
            p0 = (1-p)*(1-p)*(1-p)*(1-p);
            p1 = 4*p*(1-p)*(1-p)*(1-p);
            p2 = 6*p*p*(1-p)*(1-p);
            p3 = 4*p*p*p*(1-p);
            p4 = p*p*p*p;
        }
        if(r2<p0+p1) return false;
        else return true;
    }
}

bool C2C6_Depos(int C6_NUM)
{
    bool flag = true;
    double p0,p1,p2,p3,p4;
    int r1 = rand()/(int)(SIZE);
    double r2 = ((double)rand()/RAND_MAX);
    double NUM = SIZE - C6_NUM;
    double p = NUM/SIZE;
    if(r1<4)
    {
        if(NUM >= 2)
        {
            p0 = (1-p)*(1-p);
            p1 = 2*p*(1-p);
            p2 = p*p;
            if(r2<p2) flag = false;
        }
    }
    else if(r1 < BOARD)
    {
        if(NUM >= 3)
        {
            p0 = (1-p)*(1-p)*(1-p);
            p1 = 3*p*(1-p)*(1-p);
            p2 = 3*p*p*(1-p);
            p3 = p*p*p;
            if(r2<p3) flag = false;
        }
    }
    else
    {
        if(NUM >= 4)
        {
            p0 = (1-p)*(1-p)*(1-p)*(1-p);
            p1 = 4*p*(1-p)*(1-p)*(1-p);
            p2 = 6*p*p*(1-p)*(1-p);
            p3 = 4*p*p*p*(1-p);
            p4 = p*p*p*p;
            if(r2<p4) flag = false;
        }
    }
    return flag;
}

int main(int argc, char **argv) 
{
    double data[]={0.1,0.1,0.1,1,0.1,15,0.3,0.1,0.3,4.2,0.3,4.4,0.3,15,0.5,0.1,0.5,11.6,0.5,11.8,0.5,15};
    double C2_CON,C6_CON;
    // 循环标记
    int temp_x, temp_y;
    // 运行时间标记
    for(int i=0; i<22; i+=2)
    {
        C2_CON = data[i];
        C6_CON = data[i]*data[i+1];
        cout<<"TEST"<<i/2+1<<":"<<"C2_CON="<<C2_CON<<" R="<<data[i+1];

    clock_t start, finish, whole_start, whole_end;

    int C2_NUM = 0, C6_NUM = 0, EMPTY_NUM = SIZE;
    
    // 反应速率以及反应标记
    float all, paa, pba, pag, pbg, pad, pbd, pab;
    float AA, BA, AG, BG, AD, BD, AB;
    // 随机数标记
    float r1;
    srand(time(NULL));
    // 初始化C2 C6数量标记
    C2_NUM = 0; C6_NUM = 0; EMPTY_NUM = SIZE;

    // 开始全局计时
    whole_start = clock();
    // 做 MCTime 时间做 SIZE 次反应
    for(temp_x = 0; temp_x < MCTime; temp_x++){
        
        // 开始局部计时
        start = clock();
        
        for(temp_y = 0; temp_y < SIZE; temp_y++){
            // 计算当前反应速率
            AA = CALC_AA(EMPTY_NUM,C2_CON);
            BA = CALC_BA(EMPTY_NUM,C6_CON);
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
            // 生成随机数用于判断发生哪个反应
            r1 = ((float)rand()/RAND_MAX);
            //cout<<r1<<endl;
            // 发生反应C2吸附
            if(r1<=paa){
                if(EMPTY_NUM > 0){
                    C2_NUM++;
                    EMPTY_NUM--;
                }
                //alert(1, C2_NUM, C6_NUM, EMPTY_NUM);
            }
            // 发生反应C6吸附
            else if(r1<=paa+pba){
                //genSurface(C2_NUM, C6_NUM);
                if(EMPTY_NUM >= 3){
                    if(C6_Adsorb){
                        C6_NUM = C6_NUM + 3;
                        EMPTY_NUM = EMPTY_NUM - 3;
                    }
                }
                //alert(2, C2_NUM, C6_NUM, EMPTY_NUM);
            }
            // C2脱吸附
            else if(r1<=paa+pba+pag){
                if(C2_NUM > 0){
                    C2_NUM--;
                    EMPTY_NUM++;
                }
                //alert(3, C2_NUM, C6_NUM, EMPTY_NUM);
            }
            // C2沉积
            else if(r1<=paa+pba+pag+pad){
                if(C2_NUM > 0){
                    C2_NUM--;
                    EMPTY_NUM++;
                }
                //alert(4, C2_NUM, C6_NUM, EMPTY_NUM);
            }
            // C6脱付
            else if(r1<=paa+pba+pag+pad+pbg){
                if(C6_NUM > 2){
                    C6_NUM-=3;
                    EMPTY_NUM+=3;
                }
                //alert(5, C2_NUM, C6_NUM, EMPTY_NUM);
            }
            // C6沉积
            else if(r1<=paa+pba+pag+pad+pbg+pbd){
                if(C6_NUM > 2){
                    C6_NUM-=3;
                    EMPTY_NUM+=3;
                }
                //alert(6, C2_NUM, C6_NUM, EMPTY_NUM);
            }else{// 双分子反应
                //genSurface(C2_NUM, C6_NUM);
                if(C6_NUM >= 3 && C2_NUM >= 1){
                    if(C2C6_Depos){
                        C2_NUM--;
                        C6_NUM-=3;
                        EMPTY_NUM+=4;
                    }
                }
                //alert(7, C2_NUM, C6_NUM, EMPTY_NUM);
            }
        }
        // 结束局部计时
        finish = clock();
        // 输出反应后各分子个数
        //cout<<"C2:"<<C2_NUM<<" C6_NUM:"<<C6_NUM<<" EMPTY_NUM:"<<EMPTY_NUM<<endl;
        // 输出局部用时
        //cout<<"Using time:"<<(double)(finish - start) / CLOCKS_PER_SEC<<endl;
    }
    // 结束全局计时
    whole_end = clock();
    // 输出整体用时
    cout<<endl<<"####################"<<endl<<"Whole using time:"<<(double)(whole_end - whole_start) / CLOCKS_PER_SEC<<endl;
    
    /*
    下面主要与输出有关
    */

    cout<<"C2:"<<C2_NUM<<" C6_NUM:"<<C6_NUM<<" EMPTY_NUM:"<<EMPTY_NUM<<endl<<endl;
    }

    // // 生成用于输出的矩阵
    // genSurface(C2_NUM, C6_NUM);
    // int i, j;
    // // 输出到result文件夹中之后用python显示
    // ofstream f ("result.txt");
    // for(i = 0; i < WIDTH; i++){
    //     for(j = 0; j < HEIGHT; j++){
    //         f<<surface[i][j];
    //         if(j != (HEIGHT-1)){
    //             f<<" ";
    //         }else{
    //             f<<endl;
    //         }
    //     }
    // }
    // f.close();
    
    // 结束
    getchar();
    return 1;
}