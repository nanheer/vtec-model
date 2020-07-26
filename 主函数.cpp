#include"����.h"
#include"stdafx.h"
#include<iostream>
#include<string>
#include"stdlib.h"
#include <io.h>
#include <Eigen/Dense>
using Eigen::MatrixXd;
int _tmain(int argc, _TCHAR* argv[])
{
	pobs po = new obs;
	psp3 pn = new all_sate_ephem;
	pio ion = new ionex;
	pdcb pd = new dcb;
	psh ps = new sh_file;
	//pvt pt = new result;
	double **ctable = new double*[25]; //������  
	for (int i = 0; i < 25; i++)
		ctable[i] = new double[17]; //������ 
	double coe[(15+1)*(15+1)][25];
	vector<obs> sta_obs;
	vector<string> files;//���е���dcb�ļ����ҵõ��Ĳ�վ
	vector<obs> obs_value;
	obs obs_data;
	string file_path = "E:\\GNSS\\����\\o�ļ�\\2018152";
	string strn = "E:\\GNSS\\����\\igs20035.sp3";
	string stri = "E:\\GNSS\\����\\codg1520.18i";
	string strd = "E:\\GNSS\\����\\CAS0MGXRAP_20181520000_01D_01D_DCB.BSX";
	string strs = "E:\\GNSS\\ion\\CODE\\2018\\COD20035.ION";
	string stra = "E:\\GNSS\\����\\a0.dat";
	string strc = "E:\\����\\ϵ��.dat";
	readsp3file(strn, pn);
	read_ionex(stri, ion);
	read_sh(strs, ps);
	read_dcb(strd, pd);
	read_a0(stra, ctable);
	read_coe(strc, coe);
	getFiles(file_path, pd,files);
	int n = files.size();
	//pobs obs_value = new obs[n];
	for (int i = 0; i <20; i++)
	{
		obs_data.obsdata.clear();
		if (readofile_vtec(files[i], &obs_data)){
			obs_value.push_back(obs_data);
		}
	}
	cout << "o�ļ�������" << obs_value.size()<< endl;
	vtec(obs_value,ctable,coe, pn,ion,pd,ps);
	delete po;
	delete pn;
	delete ion;
	delete pd;
	delete ps;
	for (int i = 0; i < 25; i++)
		delete[] ctable[i];
	delete[] ctable;
	return 0;
}