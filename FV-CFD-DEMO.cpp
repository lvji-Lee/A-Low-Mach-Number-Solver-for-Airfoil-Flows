// FV-CFD-DEMO.cpp : 定义控制台应用程序的入口点。
//
#pragma  warning(disable:4996)
#include "stdafx.h"
#include <cmath>

// int _tmain(int argc, _TCHAR* argv[])
// {
// 	return 0;
// }

// FV-CFD-DEMO.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <stdio.h>
#include <stdlib.h>

#define _Grid_IMax 400 // I 方向最大点数
#define _Grid_JMax 400 // J 方向最大点数

#define _GHOST_NUM 2   // 虚拟单元层数
#define pi 3.14// Π值


//定义远场边界条件，前面加_则默认为全局变量
double _alpha = 0;        // 迎角
double _gamma = 1.4;        // 比热比
double _M_infty = 2.0;      // 远场马赫数
double _rho_infty = 1;    //远场的无量纲密度    

//定义时间推进的相关初始参数
double _deltaT = 0.000001;
int _time_level;

// 二维平面上点/向量的结构体
typedef struct Point2D
{
	double x;
	double y;
} Pnt2D;

// 二维流动守恒变量的结构体
typedef struct ConserVar2D
{
	double rho;
	double rhoE;
	double p;
	//压强P可由其他量定义
	double u;       //i方向速度
	double v;       //j方向速度
	double rhou;
	double rhov;
} CV2D;

//边界条件的定义
typedef struct BonCon2D
{

	int I0;        //i的起始值
	int IM;        //i的最大值
	int J0;        //j的起始值
	int JM;        //j的最大值
	int type;      //边界条件类型

}BC2D;

BC2D _domain_boundaries[100];          //储存边界条件的数组结构体

// 节点数目
int _node_number[3];


//边界条件数目
int _bc_number;


// 存储网格节点坐标的变量
Pnt2D _grid_node[_Grid_JMax][_Grid_IMax];
// 存储单元边法向向量的变量，是一个三维结构体数组，第一个下标表示方向；
// 若第一个下标为0，则表示i方向单元边界上的法向向量,第一个下标为1：表示j方向单元边界上的法向向量
//并且对于数组的每一个单元来说，储存的是一个含有x值和y值的结构体，或者说看成向量
Pnt2D _cell_normal[2][_Grid_JMax][_Grid_IMax];
// 存储单元面积的变量
double _cell_volume[_Grid_JMax - 1][_Grid_IMax - 1];
// 存储每一个单元上局部时间步长的变量
double _cell_time_step[_Grid_JMax][_Grid_IMax];
// 存储单元上解的均值的变量，第一个下标是时间层下标,0：n时刻，1：n+1时刻
CV2D _Solutions[2][_Grid_JMax + 2 * _GHOST_NUM][_Grid_IMax + 2 * _GHOST_NUM];
// 存储单元边界上通量的变量，第一个下标为0，表示i方向单元边界上的通量,第一个下标为1：表示j方向单元边界上的通量
CV2D _Fluxes[2][_Grid_JMax + 1][_Grid_IMax + 1];
// 存储单元上的残差
CV2D _RHS[_Grid_JMax][_Grid_IMax];


//取得n时刻（j，i）单元上的守恒变量值
#define _SOL(n,j,i)  _Solutions[(n)% 2][(j) + _GHOST_NUM][(i) + _GHOST_NUM]     //定义了一个考虑了虚拟单元之后的宏,相当于作了一个坐标系的平移变换
//使得下标ij能够直接表示真实单元的坐标（不需考虑为了表示虚拟单元而产生负值的问题）
//例如两层虚拟单元，那么就从i=2，j=2开始成为真实单元    所有地方全部都用_SOL，不用_Solutions!!!!!!!!!!!!!!!!!


CV2D _residual;     //残差的范数

//取得n时刻（j,i）节点上的坐标值
#define _GRD(j,i) _grid_node[j][i]

///取得（j,i）单元上的面积值
#define _VOL(dir,j,i) _cell_volume[j][i]

//取得dir方向编号为（j，i）的边界法向向量值
#define _VEC(dir,j,i) _cell_normal[dir][j][i]

//取得dir方向编号为（j，i）的边界通量值
#define _FLX(dir,j,i) _Fluxes[dir][j][i]

//取得（j,i）单元上的RHS值
#define _RHS(j,i) _RightHandSide[j][i]

//取得（j,i）单元上依据CFL条件计算得到的时间步长值
#define _DT(j,i) _cell_time_step[j][i]

//所有单元上全部时间步长中的最小值
double _min_delta_time = 0.000001;


// 读取参数文件
int read_parameter()
{
	return 0;
}

// 读取网格文件
int read_grid()
{
	FILE* file;
	int blk_num;

	file = fopen("N0012.dat", "r");
	if (file == NULL)
	{
		printf("Failed to open the file.\n");
		return 0;
	}

	if (fscanf(file, "%d", &blk_num) != EOF)
	{
		printf("BLK NUMBER : %d\n", blk_num);
	}

	if (fscanf(file, "%d", &_node_number[0]) != EOF)
	{
		printf("I NUMBER : %d\n", _node_number[0]);
	}

	if (fscanf(file, "%d", &_node_number[1]) != EOF)
	{
		printf("J NUMBER : %d\n", _node_number[1]);
	}
	int k;
	if (fscanf(file, "%d", &k) != EOF)
	{
		printf("K NUMBER : %d\n", k);
	}
	printf("%d %d\n", _node_number[0], _node_number[1]);
	for (int j = 0; j < _node_number[1]; ++j)
	{
		for (int i = 0; i < _node_number[0]; ++i)
		{
			double coord;
			if (fscanf(file, "%lf", &coord) != EOF)
			{
				_grid_node[j][i].x = coord;

			}

		}
	};

	for (int j = 0; j < _node_number[1]; ++j)
	{
		for (int i = 0; i < _node_number[0]; ++i)
		{
			double coord;
			if (fscanf(file, "%lf", &coord) != EOF)
			{
				_grid_node[j][i].y = coord;

			}

		}
	}
	printf("横向节点数为：%d,           纵向节点数为：%d\n", _node_number[0], _node_number[1]);
	fclose(file);
	return 0;
}


// 读取边界条件文件
int read_bc()
{
	FILE* fp0;
	fp0 = fopen("N0012.inp", "r");

	if (fp0 == NULL)
	{
		printf("Failed to open the boundary condition file.\n");
		return -1;  // 返回-1表示操作失败
	}

	// 跳过前面4行描述性信息
	char buffer[256];
	for (int i = 0; i < 4; ++i)
	{
		if (fgets(buffer, sizeof(buffer), fp0) == NULL)
		{
			printf("Failed to read the header information.\n");
			fclose(fp0);
			return -1;  // 返回-1表示操作失败
		}
	}

	// 读取文件中的边界条件数目
	if (fscanf(fp0, "%d", &_bc_number) != 1)
	{
		printf("Failed to read the number of boundary conditions.\n");
		fclose(fp0);
		return -1;  // 返回-1表示操作失败
	}

	// 读取每一个边界条件并存储到 _domain_boundaries 数组中
	for (int i = 0; i < _bc_number; ++i)
	{
		if (fscanf(fp0, "%d %d %d %d %d",
			&_domain_boundaries[i].I0,
			&_domain_boundaries[i].IM,
			&_domain_boundaries[i].J0,
			&_domain_boundaries[i].JM,
			&_domain_boundaries[i].type)
			!= 5)                        //注意内部边界条件的存在，会使得_bc_number发生变化
		{

			printf("Failed to read boundary condition %d.\n", i + 1);
			fclose(fp0);
			return -1;  // 返回-1表示操作失败
		}
		printf("%d %d %d %d %d\n", _domain_boundaries[i].I0, _domain_boundaries[i].IM,
			_domain_boundaries[i].J0, _domain_boundaries[i].JM, _domain_boundaries[i].type);
		if (_domain_boundaries[i].type == -1)
		{
			_bc_number = _bc_number + 1;
		}
	}

	fclose(fp0);
	return 0;  // 返回0表示操作成功
}

// 写出流场数据（以tecplot格式）
int write_field()
{
	int n = _time_level ;
	FILE* fp1;
	fp1 = fopen("field_juxing.dat", "w");
	if (fp1 == NULL)
	{
		printf("Failed to open the file.\n");
		return 0;
	}
	//fprintf(fp1, "Title =\"Field\"\n ");
	fprintf(fp1, "VARIABLES =\"X\",\"Y\",\"rho\",\"rhou\",\"rhov\",\"rhoE\",\"p\",\"Mach\"\n ");              //t通过块儿的格式输出流场数据
	fprintf(fp1, "ZONE T=ZONE1,  I = %d,J = %d,  DATAPACKING=BLOCK,VARLOCATION=([3-8]=CELLCENTERED) \n ", _node_number[0], _node_number[1]);
	//给出数据文件的文件头定义

	for (int j = 0; j < _node_number[1]; j++)
	{

		for (int i = 0; i < _node_number[0]; i++)
		{
			float corx = _grid_node[j][i].x;

			fprintf(fp1, "%f\n  ", corx);             //以tecplot形式存储节点坐标
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}

	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿

	for (int j = 0; j < _node_number[1]; j++)
	{

		for (int i = 0; i < _node_number[0]; i++)
		{
			float cory = _grid_node[j][i].y;

			fprintf(fp1, "%f\n  ", cory);             //以tecplot形式存储节点坐标
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}

	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	//存入密度：
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			fprintf(fp1, "%f  \n", _SOL(n, j, i).rho);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	//存入rhou：
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			fprintf(fp1, "%f \n ", _SOL(n, j, i).rhou);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	//存入rhov
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			fprintf(fp1, "%f\n  ", _SOL(n, j, i).rhov);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			fprintf(fp1, "%f \n ", _SOL(n, j, i).rhoE);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	//存储p
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			fprintf(fp1, "%f \n ", _SOL(n, j, i).p);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fprintf(fp1, "%\n\n\n "); //隔开不同的块儿
	//存马赫数Mach
	for (int j = 0; j < _node_number[1] - 1; j++)
	{

		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			double c = sqrt(_gamma * _SOL(n, j, i).p / _SOL(n, j, i).rho);
			double Ma = 1;// sqrt(_SOL(n, j, i).u * _SOL(n, j, i).u + _SOL(n, j, i).v * _SOL(n, j, i).v) / c;

			fprintf(fp1, "%f \n ", Ma);             //
			//printf(" %f %f\n", _grid_node[j][i].x, _grid_node[j][i].y);
		}
	}
	fclose(fp1);
	return 0;
}


//写出物面数据
int write_wall()
{
	int n = _time_level;
	// 1. 打开要写入数据的文件
	FILE* fp2;
	fp2 = fopen("wall_juxing.dat", "w");
	if (fp2 == NULL)
	{
		// 如果文件打开失败，输出错误信息并返回
		printf("Failed to open the wall data file.\n");
		return 0;
	}

	// 2. 写入 Tecplot 文件头
	fprintf(fp2, "TITLE = \"Wall Data\"\n");
	fprintf(fp2, "VARIABLES = \"X\", \"Y\", \"rho\", \"rhou\", \"rhov\", \"rhoE\", \"p\", \"Mach\"\n");

	// 3. 计算物面边界条件的总点数
	int total_wall_points = 0;
	for (int ib = 0; ib < _bc_number; ++ib)
	{
		if (_domain_boundaries[ib].type == 2) // Type 2: 物面边界条件
		{
			// 判断是沿 I 不变的方向还是 J不变的方向的边界
			if (_domain_boundaries[ib].I0 == _domain_boundaries[ib].IM)                 //i不变
			{
				total_wall_points = total_wall_points + (_domain_boundaries[ib].JM - _domain_boundaries[ib].J0 + 1) - 1;    //_domain_boundaries[ib].JM - _domain_boundaries[ib].J0 + 1 是一个表达式,
				// 用来计算当前边界条件在 J 方向上的网格点个数,如果是边界网格中心，则需要再减去1
			}
			else if (_domain_boundaries[ib].J0 == _domain_boundaries[ib].JM)               //j不变                     //此时无需考虑虚拟单元
			{
				total_wall_points = total_wall_points + (_domain_boundaries[ib].IM - _domain_boundaries[ib].I0 + 1) - 1;
			}
		}
	}

	// 4. 写入 Tecplot 数据区
	fprintf(fp2, "ZONE I=%d, J=1, F=POINT\n", total_wall_points - 1);                //默认是j不变的边界

	// 5. 遍历所有物面边界条件，写入数据
	for (int ib = 0; ib < _bc_number; ++ib)
	{
		if (_domain_boundaries[ib].type == 2) // Type 2: 物面边界条件
		{
			//printf("%d   %d\n", _domain_boundaries[ib].I0 - 1, _domain_boundaries[ib].J0 - 1);
			// 判断是沿 I 方向还是 J 方向的边界
			if (_domain_boundaries[ib].I0 == _domain_boundaries[ib].IM)
			{
				int i = _domain_boundaries[ib].I0 - 1;      //i不变的单元
				//printf("该处i值为：%d\n", _domain_boundaries[ib].I0 - 1);
				for (int j = _domain_boundaries[ib].J0 - 1; j != _domain_boundaries[ib].JM - 1; ++j)   //编号默认是从1开始的，因此开始要减去1
				{
					//printf("该处j为：%d\n", j);
					// 获取节点坐标和流场解
					double x1 = _grid_node[j][i].x;
					//printf("该处横坐标2为：%lf\n", _grid_node[j][i].x);
					double y1 = _grid_node[j][i].y;
					double x2 = _grid_node[j + 1][i].x;
					double y2 = _grid_node[j + 1][i].y;
					double xx = (x1 + x2) / 2;
					double yy = (y1 + y2) / 2;
					//printf("%lf    %d\n", _grid_node[3][0].x, _domain_boundaries[ib].JM-1 );
					   //计算边界的中心点处的坐标     
					//计算边界中点的值，计算方式为相邻两个单元的值相加除以2
					//printf("物面边界是I不变的\n");
					double rho = (_SOL(n, j, i).rho + _SOL(n, j, i - 1).rho) / 2;
					double u = (_SOL(n, j, i).u + _SOL(n, j, i - 1).u) / 2;
					double v = (_SOL(n, j, i).v + _SOL(n, j, i - 1).v) / 2;
					double rhou = rho * u;
					double rhov = rho * v;
					double rhoE = (_SOL(n, j, i).rhoE + _SOL(n, j, i - 1).rhoE) / 2;
					double p = (_SOL(n, j, i).p + _SOL(n, j, i - 1).p) / 2;
					// 将数据写入文件
					fprintf(fp2, "%f %f %f %f %f %f %f %f\n", xx, yy, rho, rhou, rhov, rhoE, p, _M_infty);
					//printf("各项数据依次为：%f %f %f %f %f %f %f %f\n", xx, yy, rho, u, v, rhoE, p, _M_infty);
				}
			}
			else if (_domain_boundaries[ib].J0 == _domain_boundaries[ib].JM)
			{
				int j = _domain_boundaries[ib].J0 - 1;
				for (int i = _domain_boundaries[ib].I0 - 1; i < _domain_boundaries[ib].IM - 1; ++i)        //由于边界条件里面的点编号默认是从1开始的，因此这里是要-1
				{
					// 获取节点坐标和流场解
					double x1 = _grid_node[j][i].x;
					double y1 = _grid_node[j][i].y;
					double x2 = _grid_node[j][i + 1].x;

					double y2 = _grid_node[j][i + 1].y;
					double x = (x1 + x2) / 2;
					double y = (y1 + y2) / 2;           //计算边界的中心点处的坐标
					//printf("物面边界是J不变的\n");
					//计算边界中点的值，计算方式为相邻两个单元的值相加除以2
					double rho = (_SOL(n, j - 1, i).rho + _SOL(n, j, i).rho) / 2;
					double rhou = (_SOL(n, j - 1, i).rhou + _SOL(n, j, i).rhou) / 2;
					double rhov = (_SOL(n, j - 1, i).rhov + _SOL(n, j, i).rhov) / 2;
					double rhoE = (_SOL(n, j - 1, i).rhoE + _SOL(n, j, i).rhoE) / 2;
					double p = (_SOL(n, j - 1, i).p + _SOL(n, j, i).p) / 2;
					//printf("该处虚拟网格密度为：%lf\n", _SOL(0, j - 1, i).rho);
					// 将数据写入文件
					fprintf(fp2, "%f %f %f %f %f %f %f %f\n", x, y, rho, rhou, rhov, rhoE, p, _M_infty);
				}
			}
		}

	}

	// 6. 关闭文件
	fclose(fp2);
	return 0;
}



// 计算metrics的值（即计算单元边界上的法向向量）

// 预定义一个计算二维向量的长度的函数(不要了)


// 计算metrics的值（即计算单元边界上的法向向量）
int compute_metrics()
{
	// 计算 i 方向单元边界上的法向向量
	for (int j = 0; j < _node_number[1] - 1; j++)
	{
		for (int i = 0; i < _node_number[0]; i++)
		{
			// 获取相邻节点的坐标
			Pnt2D p1 = _grid_node[j][i];
			Pnt2D p2 = _grid_node[j + 1][i];

			// 计算边的向量
			double vec_x = p2.x - p1.x;
			double vec_y = p2.y - p1.y;

			// 计算法向量（顺时针方向）
			double normal_x = vec_y;
			double normal_y = -vec_x;

			// 归一化法向量
			//double length = 1.0;                 //sqrt(normal_x * normal_x + normal_y * normal_y);无需归一化处理，要的就是有长度的法向量
			//normal_x /= length;
			//normal_y /= length;

			// 存储对应点（j,i）处的法向量
			_cell_normal[0][j][i] = { normal_x, normal_y };
		}
	}

	// 计算 j 方向单元边界上的法向向量
	for (int j = 0; j < _node_number[1]; j++)
	{
		for (int i = 0; i < _node_number[0] - 1;i++)
		{
			// 获取相邻节点的坐标
			Pnt2D p1 = _grid_node[j][i];
			Pnt2D p2 = _grid_node[j][i + 1];

			// 计算边的向量
			double vec_x = p2.x - p1.x;
			double vec_y = p2.y - p1.y;

			// 计算法向量（逆时针方向）
			double normal_x = -vec_y;
			double normal_y = vec_x;

			// 归一化法向量
			//double length = sqrt(normal_x * normal_x + normal_y * normal_y);
			//normal_x /= length;
			//normal_y /= length;

			// 存储法向量
			_cell_normal[1][j][i] = { normal_x, normal_y };
		}
	}
	
	return 0;
}


// 计算单元的面积值
int compute_volume()
{
	for (int j = 0; j < _node_number[1] - 1; j++)
	{
		for (int i = 0; i < _node_number[0] - 1; i++)
		{
			// 获取四个顶点的坐标
			Pnt2D p1 = _grid_node[j][i];
			Pnt2D p2 = _grid_node[j][i + 1];
			Pnt2D p3 = _grid_node[j + 1][i + 1];
			Pnt2D p4 = _grid_node[j + 1][i];

			// 使用多边形面积公式计算单元面积
			double area = 0.5 * fabs(
				p1.x * p2.y + p2.x * p3.y + p3.x * p4.y + p4.x * p1.y -
				(p2.x * p1.y + p3.x * p2.y + p4.x * p3.y + p1.x * p4.y)
			);

			_cell_volume[j][i] = area;
			//printf("面积值为%20.16e\n\n\n", _cell_volume[j][i]);
		}
	}
	return 0;
}

// 初始化
int initialize()
{
	compute_metrics();
	compute_volume();
	FILE* fp;
	//fp = fopen("111.txt", "w");
	double u_infty = _M_infty * cos(_alpha * pi / 180);  // 远场速度的计算
	double v_infty = _M_infty * sin(_alpha * pi / 180);
	double p_infty = 1.0 / _gamma;  // 浮点数除法
	double E_infty = p_infty / ((_rho_infty) * (_gamma - 1)) + (u_infty * u_infty+ v_infty* v_infty) / 2.0;  //  指数运算改为乘以自身,不然会报错？？

	for (int j = 0; j < _node_number[1]-1+4 ; j++)
	{
		for (int i = 0; i < _node_number[0]-1+4 ; i++)                     //不需要考虑虚拟单元，只初始化真实单元的值即可，因此用宏变量_SOL		
		{
			_SOL(0, j-2, i-2).p = p_infty;
			_SOL(0, j - 2, i - 2).rho = _rho_infty;
			_SOL(0, j - 2, i - 2).u = u_infty;
			_SOL(0, j - 2, i - 2).v = v_infty;
			_SOL(0, j - 2, i - 2).rhou = _rho_infty * u_infty;  // 初始动量 x 方向
			_SOL(0, j - 2, i - 2).rhov = _rho_infty * v_infty;  // 初始动量 y 方向
			_SOL(0, j - 2, i - 2).rhoE = _rho_infty * E_infty;  // 初始总能量
			
			//printf( "密度为：%lf    v速度为：%lf\n",_SOL(0, j, i).rho, _SOL(0, j, i).u);
			//printf("速度为：%lf    v速度为：%lf\n", _SOL(0, j, i).u, _SOL(0, j, i).v);
			//printf("u速度为：%lf    v速度为：%lf\n", _SOL(0, j, i).u, _SOL(0, j, i).v);
			//printf("u速度为：%lf    v速度为：%lf\n", _SOL(0, j, i).u, _SOL(0, j, i).v);
		}

	}
	//fclose(fp);
	return 0;
}

// 设置边界条件
//特别要注意的是，边界条件设置不是目的，而是要通过边界条件来把虚拟单元上的值给定，方便后续计算（根据边界条件给虚拟单元赋值）
//例如无粘物面边界条件要求速度关于边界镜像对称，则据此条件即可设置好对应的虚拟单元上的值
int set_boundary_conditions()                       //需要完成哪些设置：内部边界，物面边界，远场边界
{
	//FILE* fp;
	//fp = fopen("ceshi11.txt", "w");
	//要注意每次计算后边界条件都要更新一次，因此外面还要再嵌套一个关于时间推进的大循环（等写完内部内容后再写）
	int  n = _time_level;       //先假设为n时刻；

	//默认物面边界为i不变，内部边界为j不变
	//先设置物面边界条件
	int i = 0;
	for (int j = 0; j != _node_number[1] - 1; j++)
	{
		//下面通过边界条件为虚拟单元赋值：
		_SOL(n, j, i - 1).rho = _SOL(n, j, i).rho;
		_SOL(n, j, i - 2).rho = _SOL(n, j, i + 1).rho;     //密度
		_SOL(n, j, i - 1).p = _SOL(n, j, i).p;
		_SOL(n, j, i - 2).p = _SOL(n, j, i + 1).p;     //压强

		/*double vx = _grid_node[j + 1][i].x - _grid_node[j][i].x;
		double vy = _grid_node[j+1][i].y - _grid_node[j][i].y;
		double m = vx * vx + vy * vy;
		_SOL(0, j, i-1).u = (_SOL(0, j, i).u * vx * vx + 2 * vx * vy * _SOL(0, j, i).v - _SOL(0, j, i).u * vy * vy) / m;
		_SOL(0, j , i-1).v = (_SOL(0, j, i).v * vy * vy + 2 * vx * vy * _SOL(0, j, i).u - _SOL(0, j, i).v * vx * vx) / m;
		
		_SOL(0, j , i-2).u = (_SOL(0, j, i+1).u * vx * vx + 2 * vx * vy * _SOL(0, j , i+1).v - _SOL(0, j , i+1).u * vy * vy) / m;
		_SOL(0, j , i-2).v = (_SOL(0, j, i+1).v * vy * vy + 2 * vx * vy * _SOL(0, j , i+1).u - _SOL(0, j , i+1).v * vx * vx) / m;*/
		
		//printf("物面是I不变的边界！！！！\n");

		//速度的设置根据法向分量为零，即关于物面镜像对称来确定  but。。是需要设置速度还是直接设置rhou，rhov。。等量？
		//若边界为i不变，则法向量取i变化方向(即第一个下标为0)的第一组
		double len1= sqrt(_cell_normal[0][j][i].x * _cell_normal[0][j][i].x + _cell_normal[0][j][i].y * _cell_normal[0][j][i].y);
		double voc_nor_1 = _SOL(n, j, i).u * _cell_normal[0][j][i].x /len1+ _SOL(n, j, i).v * _cell_normal[0][j][i].y/len1;       //计算法向速度模值
		double voc_nor_1u = voc_nor_1 * _cell_normal[0][j][i].x/len1;
		double voc_nor_1v = voc_nor_1 * _cell_normal[0][j][i].y / len1;
		//printf("该处（%d,%d）的法向量为：(%20.16e,%20.16e)\n\n\n", i, j, _cell_normal[0][j][i].x, _cell_normal[0][j][i].y);
		_SOL(n, j, i - 1).u = _SOL(n, j, i).u - 2 * voc_nor_1u;
		_SOL(n, j, i - 1).v = _SOL(n, j, i).v - 2 * voc_nor_1v;
		//printf("该处虚拟网格的速度为：%20.16e", _SOL(n, j, i).u);
		double len2 = sqrt(_cell_normal[0][j][i].x * _cell_normal[0][j][i+1].x + _cell_normal[0][j][i].y * _cell_normal[0][j][i+1].y);
		double voc_nor_2 = _SOL(n, j, i).u * _cell_normal[0][j][i + 1].x/len2+ _SOL(n, j, i + 1).v * _cell_normal[0][j][i + 1].y/len2;       //计算法向速度模值
		double voc_nor_2u = voc_nor_2 * _cell_normal[0][j][i + 1].x/len2;
		double voc_nor_2v = voc_nor_2 * _cell_normal[0][j][i + 1].y/len2;
		_SOL(n, j, i - 2).v = _SOL(n, j, i + 1).v - 2 * voc_nor_1v;
		_SOL(n, j, i - 2).u = _SOL(n, j, i + 1).u - 2 * voc_nor_1u;
		//printf("该处的法向向量为：(%lf  ,  %lf)\n", _cell_normal[0][j][i].x, _cell_normal[0][j][i].y);
		//printf("该处的虚拟网格速度为：\nU=%lf\nV=%lf\n", _SOL(n, j, i - 1).u, _SOL(n, j, i - 1).v);
		//计算rhoE

		double E1 = _SOL(n, j, i - 1).p / ((_gamma - 1) * _SOL(n, j, i - 1).rho) + (_SOL(n, j, i - 1).u * _SOL(n, j, i - 1).u + _SOL(n, j, i - 1).v * _SOL(n, j, i - 1).v) / 2;
		double E2 = _SOL(n, j, i - 2).p / ((_gamma - 1) * _SOL(n, j, i - 2).rho) + (_SOL(n, j, i - 2).u * _SOL(n, j, i - 2).u + _SOL(n, j, i - 2).v * _SOL(n, j, i - 2).v) / 2;
		
		_SOL(n, j, i - 1).rhoE = _SOL(n, j, i - 1).rho * E1;
		_SOL(n, j, i - 2).rhoE = _SOL(n, j, i - 2).rho * E2;

		_SOL(n, j, i - 1).rhou = _SOL(n, j, i - 1).rho * _SOL(n, j, i - 1).u;
		_SOL(n, j, i - 2).rhou = _SOL(n, j, i - 2).rho * _SOL(n, j, i - 2).u;

		_SOL(n, j, i - 1).rhov = _SOL(n, j, i - 1).rho * _SOL(n, j, i - 1).v;
		_SOL(n, j, i - 2).rhov = _SOL(n, j, i - 2).rho * _SOL(n, j, i - 2).v;
		//fprintf(fp, "%lf     %lf         %d     %d\n\n\n", _SOL(n, j, i - 1).p, _SOL(n, j, i).p, j,i);
	}

	//设置内部周期性边界条件：
	int j0 = 0;
	for (int i = 0; i< _node_number[0] - 1; i++)
	{
	     //定义边界附近的单元，单元数要比节点数少1
		int jm = _node_number[1] - 2;       //同上
		_SOL(n, j0 - 1, i).rho = _SOL(n, jm, i).rho;
		_SOL(n, j0 - 2, i).rho = _SOL(n, jm - 1, i).rho;     //密度

		_SOL(n, j0 - 1, i).u = _SOL(n, jm, i).u;
		_SOL(n, j0 - 2, i).u = _SOL(n, jm - 1, i).u;       //rhou

		_SOL(n, j0 - 1, i).p = _SOL(n, jm, i).p;
		_SOL(n, j0 - 2, i).p = _SOL(n, jm - 1, i).p;       //rhou


		_SOL(n, j0 - 1, i).v = _SOL(n, jm, i).v;
		_SOL(n, j0 - 2, i).v = _SOL(n, jm - 1, i).v;       //rhou

		_SOL(n, j0 - 1, i).rhou = _SOL(n, jm, i).rhou;
		_SOL(n, j0 - 2, i).rhou = _SOL(n, jm - 1, i).rhou;       //rhou

		_SOL(n, j0 - 1, i).rhov = _SOL(n, jm, i).rhov;
		_SOL(n, j0 - 2, i).rhov = _SOL(n, jm - 1, i).rhov;       //rhov

		_SOL(n, j0 - 1, i).rhoE = _SOL(n, jm, i).rhoE;
		_SOL(n, j0 - 2, i).rhoE = _SOL(n, jm - 1, i).rhoE;       //rhoE

		_SOL(n, jm+1, i).rho = _SOL(n, j0, i).rho;
		_SOL(n, jm+ 2, i).rho = _SOL(n, j0+1, i).rho;     //密度

		_SOL(n, jm + 1, i).u = _SOL(n, j0, i).u;
		_SOL(n, jm + 2, i).u = _SOL(n, j0 + 1, i).u;       

		_SOL(n, jm + 1, i).p = _SOL(n, j0, i).p;
		_SOL(n, jm + 2, i).p = _SOL(n, j0 + 1, i).p;


		_SOL(n, jm + 1, i).v = _SOL(n, j0, i).v;
		_SOL(n, jm + 2, i).v = _SOL(n, j0 + 1, i).v;       

		_SOL(n, jm + 1, i).rhou = _SOL(n, j0, i).rhou;
		_SOL(n, jm + 2, i).rhou = _SOL(n, j0 + 1, i).rhou;      //rhou

		_SOL(n, jm + 1, i).rhov = _SOL(n, j0, i).rhov;
		_SOL(n, jm + 2, i).rhov = _SOL(n, j0 + 1, i).rhov;        //rhov

		_SOL(n, jm + 1, i).rhoE = _SOL(n, j0, i).rhoE;
		_SOL(n, jm + 2, i).rhoE = _SOL(n, j0 + 1, i).rhoE;        //rhoE
		//printf("i 和 j 方向的节点数分别为：%d  %d\n\n\n", _node_number[0], _node_number[1]);
		
		//printf("两种算法的差别：：%20.16e,    %20.16e\n\n\n", _SOL(n, jm + 1, i).v, _SOL(n, jm + 1, i).rhov / _SOL(n, jm + 1, i).rho);
	}
	
	//设置远场边界条件,默认也是i不变的区域
	 i = _node_number[0] - 1; //远场起始位置
		for (int j = 0; j != _node_number[1] - 1;j++)
		{
			double u_infty = _M_infty * cos(_alpha * pi / 180);  // 远场速度的计算
			double v_infty = _M_infty * sin(_alpha * pi / 180);
			double p_infty = 1.0 / _gamma;  // 浮点数除法
			double E_infty = p_infty / ((_rho_infty) * (_gamma - 1)) + (u_infty * u_infty + v_infty * v_infty) / 2.0;  //  指数运算改为乘以自身,不然会报错？？
			_SOL(n, j, i).p = p_infty;
			_SOL(n, j, i).rho = _rho_infty;
			_SOL(n, j, i).u = u_infty;
			_SOL(n, j, i).v = v_infty;
			_SOL(n, j, i).rhou = _rho_infty * u_infty;  // 初始动量 x 方向
			_SOL(n, j, i).rhov = _rho_infty * v_infty;  // 初始动量 y 方向
			_SOL(n, j, i).rhoE = _rho_infty * E_infty;  // 初始总能量

			_SOL(n, j, i+1).p = p_infty;
			_SOL(n, j, i+1).rho = _rho_infty;
			_SOL(n, j, i+1).u = u_infty;
			_SOL(n, j, i+1).v = v_infty;
			_SOL(n, j, i+1).rhou = _rho_infty * u_infty;  // 初始动量 x 方向
			_SOL(n, j, i+1).rhov = _rho_infty * v_infty;  // 初始动量 y 方向
			_SOL(n, j, i+1).rhoE = _rho_infty * E_infty;  // 初始总能量



		}
	
	//结束全部边界条件的设置工作
  //fclose(fp);
	return 0;
}

// 计算时间步长
int compute_time_step()
{
	return 0;
}

// 计算准确的通量，
// flx 存储计算得到的通量
// Q 存储守恒变量值
// n 单元边界法向量
// 此函数返回逆变速度（流动速度与n的点积）
/*double getFlux：函数返回一个?double?类型的值，这是一个逆变速度（即流体速度与法向量的点积）。
CV2D* flx：一个指向?CV2D?结构的指针，用于存储计算得到的通量。
CV2D const* Q：一个指向常量?CV2D?结构的指针，存储了守恒变量值，如密度和动量。
Pnt2D const* n：一个指向常量?Pnt2D?结构的指针，表示边界法向量。*/
double getFlux(CV2D* flx, CV2D const* Q, Pnt2D const* n)
{
	double const u = Q->rhou / Q->rho;  // 计算x方向速度分量
	double const v = Q->rhov / Q->rho;  // 计算y方向速度分量
	double const U = u * n->x + v * n->y;  // 逆变速度，即速度与法向量的点积
	double const p = Q->p;  // 计算压力

	//计算每个单元中心的通量大小，主要目的是为了在下一个函数中进行调用
	flx->rho = Q->rho * U;  // 计算质量通量
	flx->rhou = Q->rhou * U + p * n->x;  // 计算x方向动量通量
	flx->rhov = Q->rhov * U + p * n->y;  // 计算y方向动量通量
	flx->rhoE = (Q->rhoE + p) * U;  // 计算能量通量

	return U;  // 返回逆变速度(单元中心的逆变速度)
}

// 用L-F格式计算迎风通量。
// 这一步要计算的是单元边界上的通量，基于上面已经得到的单元中心通量的值，采用迎风格式计算
// flx 存储计算得到的通量
// ql 存储左单元（边界法向量规定从左单元指向右单元）的守恒变量值
// qr 存储右单元（边界法向量规定从左单元指向右单元）的守恒变量值
// n 单元边界法向量，从左单元指向右单元，大小为单元边界的长度
int LFflux(CV2D* flx, CV2D* ql, CV2D* qr, Pnt2D* n)
{
	CV2D FL;  // 由左(下)单元值计算的通量值
	CV2D FR;  // 由右（上）单元值计算的通量值
	CV2D Qc;  // 左右单元中间边上的变量值

	double UL = fabs(getFlux(&FL, ql, n));  // 调用getFlux函数进行计算，并把计算过程中得到的值赋到FL，并返回左单元计算得到的逆变速度值
	double UR = fabs(getFlux(&FR, qr, n));  // 调用getFlux函数进行计算，并把计算过程中得到的值赋到FR，并返回左单元计算得到的逆变速度值
	double len = sqrt(n->x * n->x + n->y * n->y);  // 单元边界的长度
	double c1 = sqrt(_gamma * ql->p / ql->rho);    //分别计算左右单元上的声速
	double c2 = sqrt(_gamma * qr->p / qr->rho);
	double c = (c1 + c2) / 2.0;             //边界的声速
	//下面是计算单元边界上的值：
	Qc.rho = (ql->rho + qr->rho) / 2.0;  // 计算左右单元密度的平均值
	Qc.rhou = (ql->rhou + qr->rhou) / 2.0;  // 计算左右单元x方向动量的平均值
	Qc.rhov = (ql->rhov + qr->rhov) / 2.0;  // 计算左右单元y方向动量的平均值
	Qc.rhoE = (ql->rhoE + qr->rhoE) / 2.0;  // 计算左右单元能量的平均值
	//Qc.u = (ql->u + qr->u) / 2.0;
	//Qc.v = (ql->v + qr->v) / 2.0;


	double UaC = fabs(Qc.rhou / Qc.rho * n->x + Qc.rhov / Qc.rho * n->y) + c * len;
	// 曲线坐标系上方程Jacobian矩阵的特征值的最大值
	//printf("第%d次迭代的数据：\n\n", _time_level+1);
	//printf("该处密度为：%20.16e\n", Qc.rho);
	//printf("该处逆变速度为：%20.16e\n", UaC);
	// Lax-Friedrich格式计算通量。
	flx->rho = (FL.rho + FR.rho) / 2.0 - UaC / 2.0 * (qr->rho - ql->rho);
	flx->rhou = (FL.rhou + FR.rhou) / 2.0 - UaC / 2.0 * (qr->rhou - ql->rhou);
	flx->rhov = (FL.rhov + FR.rhov) / 2.0 - UaC / 2.0 * (qr->rhov - ql->rhov);
	flx->rhoE = (FL.rhoE + FR.rhoE) / 2.0 - UaC / 2.0 * (qr->rhoE - ql->rhoE);
	//printf("该处的通量为：%20.16e\n\n\n", flx->rhou);
	return 0;
}


// 计算通量
int compute_flux()
{

	int const n = _time_level;
	CV2D QL;
	CV2D QR;

	//首先调用之前的LFflux函数，用以计算单元边界的通量，并存储到_FLX数组中（具体可见上述函数需要调用的参数）
	//再代入公式进行结合时间推进的计算
	// 
	// i 方向
	//printf("第%d次迭代的数据：\n\n", _time_level + 1);
	for (int j = 0; j != _node_number[1]-1 ; ++j)
	{
		for (int i = 0; i != _node_number[0]; ++i)
		{
			//RECON(&QL, &QR, &_SOL(N, j, i - 2), &_SOL(N, j, i - 1), &_SOL(N, j, i), &_SOL(N, j, i + 1));         //二阶格式时候会用到，一阶格式暂不需要
			//LFflux(&_FLX(0, j, i), &QL, &QR, &_VEC(0, j, i));
			LFflux(&_FLX(0, j, i), &_SOL(n, j, i - 1), &_SOL(n, j, i), &_VEC(0, j, i));
			
			//printf("该处(%d,%d)I的通量为：%20.16e\n\n\n",i,j, _FLX(0,j,i).rho);
		}
	}

	// j 方向
	for (int j = 0; j != _node_number[1]; ++j)
	{
		for (int i = 0; i != _node_number[0]-1 ; ++i)
		{
			//RECON(&QL, &QR, &_SOL(N, j - 2, i), &_SOL(N, j - 1, i), &_SOL(N, j, i), &_SOL(N, j + 1, i));
			//LFflux(&_FLX(1, j, i), &QL, &QR, &_VEC(1, j, i));
			LFflux(&_FLX(1, j, i), &_SOL(n, j - 1, i), &_SOL(n, j, i), &_VEC(1, j, i));
			//printf("该处(%d,%d)I的通量为：%20.16e\n\n\n", i, j, _FLX(0, j, i).rho);
			//printf("该处(%d,%d)J的通量为：%20.16e\n\n\n", i, j, _FLX(1, j, i).rho);
		}
	}

	return 0;
}

// 计算残差
int compute_residual()
{
	return 0;
}

// 更新下一步的解
int update_solution()
{

	int n = _time_level;
	double dt = _deltaT;
	//计算得到各单元边界的通量后，带入时间推进公式进行Q的计算
	for (int j = 0; j != _node_number[1] - 1; j++)
	{
		for (int i = 0; i != _node_number[0] - 1; i++)
		{
			//printf("%lf\n", _cell_volume[j][i]);

			_SOL(n + 1, j, i).rho = _SOL(n, j, i).rho + dt / _cell_volume[j][i] * (_FLX(0, j, i ).rho - _FLX(0, j, i+1).rho + _FLX(1, j , i).rho - _FLX(1, j+1, i).rho);
			_SOL(n + 1, j, i).rhou = _SOL(n, j, i).rhou + dt / _cell_volume[j][i] * (_FLX(0, j, i).rhou - _FLX(0, j, i+1 ).rhou + _FLX(1, j, i).rhou - _FLX(1, j+1 , i).rhou);
			_SOL(n + 1, j, i).rhov = _SOL(n, j, i).rhov + dt / _cell_volume[j][i] * (_FLX(0, j, i).rhov - _FLX(0, j, i+1 ).rhov + _FLX(1, j, i).rhov - _FLX(1, j+1 , i).rhov);
			_SOL(n + 1, j, i).rhoE = _SOL(n, j, i).rhoE + dt / _cell_volume[j][i] * (_FLX(0, j, i).rhoE - _FLX(0, j, i+1 ).rhoE + _FLX(1, j, i).rhoE - _FLX(1, j+1 , i).rhoE);
			//printf("面积值为：%e20.16\n\n\n", _cell_volume[j][i]);
			//更新速度等其他量：
			_SOL(n + 1, j, i).u = _SOL(n + 1, j, i).rhou / _SOL(n + 1, j, i).rho;
			_SOL(n + 1, j, i).v = _SOL(n + 1, j, i).rhov / _SOL(n + 1, j, i).rho;
			double uv = _SOL(n + 1, j, i).u * _SOL(n + 1, j, i).u + _SOL(n + 1, j, i).v * _SOL(n + 1, j, i).v;
			_SOL(n + 1, j, i).p = (_gamma - 1) * (_SOL(n + 1, j, i).rhoE - _SOL(n + 1, j, i).rho * uv / 2);

			/*printf("\n第%d次迭代:\n该时刻(%d,%d)单元左通量为：%20.16e\n", n + 1, i, j, _FLX(0, j, i).rhou);
			printf("该时刻(%d,%d)单元右通量为：%20.16e\n\n\n", i, j, _FLX(0, j, i + 1).rhou);
			printf("该时刻(%d,%d)单元下通量为：%20.16e\n", i, j, _FLX(1, j, i).rhou);
			printf("该时刻(%d,%d)单元上通量为：%20.16e\n\n\n", i, j, _FLX(1, j+1, i ).rhou);
			printf("此处密度为：%lf\n\n\n\n\n", _SOL(n + 1, j, i).rho);*/

			/*printf("I方向法向量为：(%lf  ,  %lf)\n", _VEC(0, j, i).x, _VEC(0, j, i).y);
			printf("I方向法向量为：(%lf  ,  %lf)\n", _VEC(0, j, i+1).x, _VEC(0, j, i+1).y);
			printf("J方向法向量为：(%lf  ,  %lf)\n", _VEC(1, j, i).x, _VEC(1, j, i).y);
			printf("J方向法向量为：(%lf  ,  %lf)\n\n\n\n", _VEC(1, j+1, i).x, _VEC(1, j+1, i).y);*/

			//printf("该时刻(%d,%d)单元的密度为：%lf", i,j,_SOL(n + 1, j, i).rho);
			/*_SOL(n, j, i).rho = _SOL(n + 1, j, i).rho;
			_SOL(n, j, i).rhou = _SOL(n+1, j, i).rhou;
			_SOL(n, j, i).rhov = _SOL(n+1, j, i).rhov;
			_SOL(n, j, i).rhoE = _SOL(n+1, j, i).rhoE;
			_SOL(n, j, i).u = _SOL(n+1, j, i).u;
			_SOL(n, j, i).v = _SOL(n+1, j, i).v;
			_SOL(n, j, i).p = _SOL(n+1, j, i).p;*/
			
		}
	}

	return 0;
}

// 检测残差是否小于设定值
int check_residual()
{
	return 0;
}

// 时间推进
int iteration()
{
	_time_level = 0;
	int it = 0;
	for (it = 0; it !=10000; it++)            //初始设置为10000，现在可以先修改一下
	{
		_time_level = it;
		set_boundary_conditions();
		compute_flux();
		//compute_residual();
		update_solution();

		if (it % 9999 == 0)                    //初始设的值为100，可以改,相当于每100步出一个结果
		{
			write_wall();
			write_field();
		}
	}
	return 0;
}

int main(int argc, char** argv)
{
	printf("Hello,CFD!");

	//read_parameter();
	read_grid();
	read_bc();
	initialize();
	//set_boundary_conditions();
	compute_metrics;              //
	compute_volume();              //
    iteration();
	//set_boundary_conditions();       //后续设置好时间推进后这里要删除哦
	//write_wall();
	//write_field();


}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
