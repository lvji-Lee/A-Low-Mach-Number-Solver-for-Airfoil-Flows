// FV-CFD-DEMO.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


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
	double rhou;
	double rhov;
	double rhoE;
} CV2D;


// 节点数目
int _node_number[3];

// 存储网格节点坐标的变量
Pnt2D _grid_node[_Grid_JMax][_Grid_IMax];
// 存储单元边法向向量的变量，第一个下标为0，表示i方向单元边界上的法向向量,每一个下标为1：表示j方向单元边界上的法向向量
Pnt2D _cell_normal[2][_Grid_JMax][_Grid_IMax];
// 存储单元面积的变量
double _cell_volume[_Grid_JMax - 1][_Grid_IMax - 1];
// 存储每一个单元上局部时间步长的变量
double _cell_time_step[_Grid_JMax][_Grid_IMax];
// 存储单元上解的均值的变量，第一个下标是时间层下标,0：n时刻，1：n+1时刻
CV2D _Solutions[2][_Grid_JMax + 2 * _GHOST_NUM][_Grid_IMax + 2 * _GHOST_NUM];
// 存储单元边界上通量的变量，第一个下标为0，表示i方向单元边界上的通量,每一个下标为1：表示j方向单元边界上的通量
CV2D _Fluxes[2][_Grid_JMax + 1][_Grid_IMax + 1];
// 存储单元上的残差
CV2D _RHS[_Grid_JMax][_Grid_IMax];

#define _SOL(n,j,i)  _Solutions[n % 2][j - _GHOST_NUM][i - _GHOST_NUM]


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

	for (int j=0; j<_node_number[1]; j++)
	{
		for (int i = 0; i < _node_number[0]; i++)
		{
			double coord;
			if (fscanf(file, "%lf", &coord) != EOF)
			{
				_grid_node[j][i].x = coord;
				printf(" %f\n", _grid_node[j][i].x);
			}
		}
	}
	fclose(file);
	return 0;
}

// 读取边界条件文件
int read_bc()
{
	return 0;
}

// 写出流场数据（以tecplot格式）
int write_field()
{
	return 0;
}

// 写出物面数据（以tecplot格式）
int write_wall()
{
	return 0;
}

// 计算metrics的值（即计算单元边界上的法向向量）
int compute_metrics()
{
	return 0;
}

// 计算单元的面积值
int compute_volume()
{
	return 0;
}

// 初始化
int initialize()
{
	compute_metrics();
	compute_volume();

	return 0;
}

// 设置边界条件
int set_boundary_conditions()
{
	return 0;
}

// 计算时间步长
int compute_time_step()
{
	return 0;
}

// 计算通量
int compute_flux()
{
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
	int it = 0;
	for (it = 0; it != 10000; ++it)
	{
		set_boundary_conditions();
		compute_flux();
		compute_residual();
		update_solution();

		if (it % 100 == 0)
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

	read_parameter();
	read_grid();
	read_bc();
	initialize();
	iteration();
	write_wall();
	write_field();


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
