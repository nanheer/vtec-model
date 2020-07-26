#include"函数.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<string>
#include <Eigen/Dense>
#include <io.h>//文件夹操作
#include <algorithm>//字符串操作
using Eigen::MatrixXd;
using namespace std;

//读取sp3文件
void readsp3file(string strn, psp3 sp3file)
{
	s_sp3_ephe onesp3;
	ifstream nfile(strn, ios::in);
	if (!nfile){
		cout << "sp3文件打开错误" << endl;
		exit(0);
	}
	cout << "开始读sp3文件" << endl;
	for (int i = 0; i < 32; i++){//prn初始化
		if (i < 9){
			sp3file->all_ephem[i].prn = "G0" + to_string(i + 1);
		}
		else{sp3file->all_ephem[i].prn = "G" + to_string(i + 1);}
	}
	string str1;
	getline(nfile, str1);
	while (str1.substr(0, 1) != "*")
	{
		getline(nfile, str1);
	}
	while (!nfile.eof())
	{
		if (str1.substr(0, 1) == "*"){
			onesp3.utime_n.year = atoi(str1.substr(3, 4).c_str());
			onesp3.utime_n.month = atoi(str1.substr(8, 2).c_str());
			onesp3.utime_n.day = atoi(str1.substr(11, 2).c_str());
			onesp3.utime_n.hour = atoi(str1.substr(14, 2).c_str());
			onesp3.utime_n.minute = atoi(str1.substr(17, 2).c_str());
			onesp3.utime_n.second = atof(str1.substr(20, 10).c_str());
		}
		else if(str1.substr(0, 2) == "PG"){
			onesp3.x = atof(str1.substr(5, 13).c_str())*1000.0;
			onesp3.y = atof(str1.substr(19, 13).c_str())*1000.0;
			onesp3.z = atof(str1.substr(33, 13).c_str())*1000.0;
			for (int i = 0; i < 32; i++){
				if (sp3file->all_ephem[i].prn == str1.substr(1, 3).c_str()){
					sp3file->all_ephem[i].sate_ephem.push_back(onesp3);
				}
			}
		}
		getline(nfile, str1);
	}
	nfile.close();
	cout << "sp3文件读取完毕" << endl;
}

//读取BDGIM 模型的非发播系数及预报周期
void read_a0(string path, double** ctable){
	ifstream ifile(path, ios::in);
	if (!ifile){
		cout << "a0文件打开错误" << endl;
		exit(0);
	}
	string str;
	int row = 0;
	getline(ifile, str);
	do{
		for (int i = 0; i < 17; i++){
			ctable[row][i] = atof(str.substr(1 + i * 16, 10).c_str())*pow(10,atof(str.substr(12 + i * 16, 4).c_str()));
		}
		row++;
		do{
			getline(ifile, str);
		} while (str.length() == 0 && !ifile.eof());
	} while (!ifile.eof());
	ifile.close();
	cout << "a0文件读取完毕" << endl;
}
//读取球谐函数的预报周期对应的系数
void read_coe(string path, double coe[][25]) {
	cout << "开始读取系数文件" << endl;
	ifstream ifile(path, ios::in);
	if (!ifile) {
		cout << "系数文件打开错误" << endl;
		exit(0);
	}
	string str;
	getline(ifile, str);
	for (int i = 0; i < (15+1)*(15+1); i++) 
	{
		getline(ifile, str);
		for (int j = 0; j < 25; j++)
	    {
			coe[i][j] = atof(str.substr(10 + j * 15, 14).c_str());
			//cout <<std::right<<setw(15)<< coe[i][j];
		}
		//cout << endl;
	} 
	ifile.close();
	cout << "系数文件读取完毕" << endl;
}
//读取o文件
bool readofile_vtec(string stro, pobs obsfile)//按照相位平滑伪距计算方式读取，即同一颗卫星放一起
{
	cout << "打开文件：" << stro << endl;
	obs_value val;//一个时间一颗卫星的观测数据
	od one_time_val;//一个时间所有卫星的观测数据
	//GPS
	int pos[4] = { -1, -1, -1, -1 };//测距码在所给数据中的位置
	double value;//对应的数值
	int prn_num;//有的测站给出的卫星prn号是G 1的形式,统一改为G01的形式
	string  prn_name;
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << stro<<"o文件打开错误" << endl;
		return false;
		//exit(0);
	}
	string str1, str_ins;//str_ins临时数据，用于str2和str3的pushback函数使用
	vector<string> str2, str3;//str2存储卫星prn，str3存储一颗卫星的所有数据
	getline(ofile, str1);
	//cout << str1 << endl;
	int a = 0;//所有的测距码类型
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "开始读文件头"<< endl;
		//cout <<str1<< endl;
		//GPS
		if (str1.substr(60, 19) == "# / TYPES OF OBSERV")
		{
			a = atoi(str1.substr(4, 2).c_str());
			int b = 0;//防止一行数据写不完所有的测距码类型
			//cout << "a的值：" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(10 + 6 * (i - b), 2) == "C1") pos[0] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "P2") pos[1] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L1") pos[2] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L2") pos[3] = i;
				if ((i+1)%9 ==0) getline(ofile, str1), b = i+1;//一行最多写9种类型的数据
			}
			if (pos[0] == -1)cout << "没有C1数据" << endl;
			if (pos[1] == -1)cout << "没有P2数据" << endl;
			if (pos[2] == -1)cout << "没有L1数据" << endl;
			if (pos[3] == -1)cout << "没有L2数据" << endl;
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		if (str1.substr(60, 11) == "MARKER NAME"){
			obsfile->obsheaddata.station = str1.substr(0, 4);
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	//cout << "读取o文件头成功" << endl;
	//getline(ofile, str1);
	//int cou=0;
	while (!ofile.eof())
	{
		//getline(ofile, str1);
		// << "str1:" << str1 << endl;
		//cout << "str1.length():" << str1.length() << endl;
		//if (str1.empty() == true && ofile.eof()) break;//最后一行为空行则为true，字符串为空时substr函数会报错
		do{
			getline(ofile, str1);
			//cout << "str1:" << str1 << endl;
		} while (str1.length() == 0 && !ofile.eof());//防止连续出现一个或者多个空行
		if (str1.empty() == true && ofile.eof()) break;//最后一行为空行则为true，字符串为空时substr函数会报错
		//cout << "str1:" << str1 << endl;
		//if (str1.substr(1, 2) != "17")cout << "str1" << str1<<endl;
		if (str1.substr(1, 2) == "18" && (str1.substr(32, 1) == "G" || str1.substr(32, 1) == "R" || str1.substr(32, 1) == "E" || str1.substr(32, 1) == "C"))// == "17"意味着读取的是17年的数据，待修改
		//if (str1.substr(32, 1) == "G" || str1.substr(32, 1) == "R" || str1.substr(32, 1) == "E")
		{
			//cou++;
			//先清空vector变量
			str2.clear();
			one_time_val.prn_list.clear();
			one_time_val.one_obs_data.clear();
			int sat = 0;
			int b = 0;//防止一行数据写不完所有的测距码类型
			one_time_val.utime_o.year = atoi(str1.substr(1, 2).c_str()) + 2000;
			one_time_val.utime_o.month = atoi(str1.substr(4, 2).c_str());
			one_time_val.utime_o.day = atoi(str1.substr(7, 2).c_str());
			one_time_val.utime_o.hour = atoi(str1.substr(10, 2).c_str());
			one_time_val.utime_o.minute = atoi(str1.substr(13, 2).c_str());
			one_time_val.utime_o.second = atof(str1.substr(16, 10).c_str());
			sat = atoi(str1.substr(30, 2).c_str());
			for (int k = 0; k < (sat-1) / 12; k++){//等于12 的整数倍时不换行，再大一才换，一共sat/12行数据，一行12个
				getline(ofile, str_ins);
				str2.push_back(str_ins);
			}
			for (int i = 0; i < sat; i++)
			{
				//先清空vector变量
				str3.clear();
				val.p1 = 0.0; val.p2 = 0.0; val.l1 = 0.0; val.l2 = 0.0;
				if (i % 12 == 0 && i != 0){ b = i; str1 = str2[i / 12 - 1]; };//换行
				//cout <<"i的值" <<i<< endl;
				for (int k = 0; k < (a-1) / 5 + 1; k++){//同理sat-1，一颗卫星的所有类型观测数据，一共a/5行数据，一行5种类型数据
					getline(ofile, str_ins);
					str3.push_back(str_ins);
				}
				if (str1.substr(32 + (i - b) * 3, 1) != "G") continue;//保证读取的数据为GPS
				prn_num = atoi(str1.substr(32 + (i - b) * 3+1, 2).c_str());//对于小于10的数字，有的测站给出的卫星prn号是G 1的形式,统一改为G01的形式
				if (prn_num < 9){
					prn_name = "G0" + to_string(prn_num);
				}
				else{ prn_name = "G" + to_string(prn_num); }
				//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
				for (int p = 0; p < 4; p++)
				{
					if (str3[pos[p] / 5].length() < (14 + 16 * (pos[p] % 5)))break;//保证读取的数据长度够长，因为atof函数的代入值不能为空
					value = atof(str3[pos[p] / 5].substr(1 + 16 * (pos[p] % 5), 14).c_str());//pos[p]种类型数据所在行数为pos[p]/5，列数为pos[p]%5
					switch (p)
					{
					case 0:val.p1 = value; break;
					case 1:val.p2 = value; break;
					case 2:val.l1 = value; break;
					case 3:val.l2 = value; break;
					}
				}
				if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l2 == 0) continue;//舍弃此颗卫星
				one_time_val.prn_list.push_back(prn_name);
				one_time_val.one_obs_data.push_back(val);
			}
			obsfile->obsdata.push_back(one_time_val);
		}
		//getline(ofile, str1);//放在最后一行防止出现读取空行用substr函数判断出错
	}
	ofile.close();
	cout << "读取o文件成功 " << endl;
	if (obsfile->obsdata.size() != 2880){ cout << "                                             历元数不等于2880   " << obsfile->obsdata.size() << endl; return false; }
	BLH gc;
	xyztoblh(&obsfile->obsheaddata.approx_coordinate, &gc);
	ofstream coo("E:\\GNSS\\数据\\coor.txt", ios::out);
	coo.setf(ios::fixed);
	coo.precision(4);
	coo << std::right << setw(15) << gc.latitude*180.0 / PI << std::right << setw(15) << gc.longitude*180.0 / PI << endl;
	coo.close();
	return true;
}

//读取电离层格网文件
void read_ionex(string stri, pio ionfile)//按照相位平滑伪距计算方式读取，即同一颗卫星放一起
{
	cout << "开始读ionex文件 " << endl;
	satedcb dcb_sate;
	stadcb dcb_sta;
	vtecmap vtec_val;
	//vtecmap val_rms;
	ifstream ofile(stri, ios::in);
	int lat = 0;//纬度个数
	int num;//一行数据包含的点数，第五行9个，其余16个，共73个
	if (!ofile){
		cout << "ionex文件打开错误" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "开始读o文件" << endl;
	string str1, str2;
	getline(ofile, str1);
	//cout << str1 << endl;
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		if (str1.substr(60, 16) == "PRN / BIAS / RMS")
		{
			//cout << "卫星dcb：" << str1 << endl;
			dcb_sate.prn = str1.substr(3, 3);
			dcb_sate.bias = atof(str1.substr(10, 6).c_str());
			dcb_sate.rms = atof(str1.substr(20, 6).c_str());
			ionfile->sate_dcb.push_back(dcb_sate);
		}
		if (str1.substr(60, 20) == "STATION / BIAS / RMS")
		{
			dcb_sta.station = str1.substr(6, 4);
			dcb_sta.bias = atof(str1.substr(30, 6).c_str());
			dcb_sta.rms = atof(str1.substr(40, 6).c_str());
			ionfile->sta_dcb.push_back(dcb_sta);
		}
		getline(ofile, str1);
	}
	cout << "读取ionex文件头成功" << endl;
	//cout <<str1<< endl;
	while (str1.substr(60, 11) != "END OF FILE")//getline函数遇到空行不再往下读
	{
		getline(ofile, str1);//跳过END OF TEC MAP或者END OF RMS MAP
		lat = 0;
		if (str1.substr(60, 16) == "START OF TEC MAP" || str1.substr(60, 16) == "START OF RMS MAP")//新的时间点对应的全天vtec
		{
			str2 = str1;//用于判断是START OF TEC MAP还是START OF RMS MAP
			getline(ofile, str1);
			if (str1.substr(60, 20) == "EPOCH OF CURRENT MAP"){
				//cout << str1 << endl;
				vtec_val.vtime.year = atoi(str1.substr(2, 4).c_str());
				vtec_val.vtime.month = atoi(str1.substr(10, 2).c_str());
				vtec_val.vtime.day = atoi(str1.substr(16, 2).c_str());
				vtec_val.vtime.hour = atoi(str1.substr(22, 2).c_str());
				vtec_val.vtime.minute = atoi(str1.substr(28, 2).c_str());
				vtec_val.vtime.second = atof(str1.substr(34, 10).c_str());
			}
			getline(ofile, str1);
			while (str1.substr(60, 20) == "LAT/LON1/LON2/DLON/H")//新的纬度对应的vtec
			{
				num = 16;
				for (int i = 0; i < 5; i++)
				{
					if (i == 4) num = 9;
					getline(ofile, str1);
					for (int j = 0; j < num; j++)
					{
						vtec_val.tec_values[lat][i * 16 + j] = atof(str1.substr(j * 5, 5).c_str());
					}
				}
				lat++;
				getline(ofile, str1);
			}
			if (str2.substr(60, 16) == "START OF TEC MAP"){
				ionfile->vtec_map.push_back(vtec_val);//vtec值
			}
			else{
				ionfile->vtec_rms.push_back(vtec_val);//vtec误差值
			}
			//getline(ofile, str1);//跳过END OF TEC MAP或者END OF RMS MAP
		}
	}
	ofile.close();
	cout << "读取ionex文件成功 " << endl;
}
//读取电离层球谐函数文件
//读取文件
void read_sh(string path, psh sh_file)
{
	ifstream  ifile(path, ios::in);
	one_coeff coe;
	one_line_coeff one_coe;
	if (!ifile){
		cout << "电离层文件打开错误" << endl;
		exit(0);
	}
	string str;
	getline(ifile, str);
	do{
		//数据头
		while (str.substr(0, 42) != "DEGREE  ORDER    VALUE (TECU)   RMS (TECU)")
		{
			if ((str.substr(0, 41) == "COORDINATES OF EARTH-CENTERED DIPOLE AXIS"))
			{
				getline(ifile, str);
				coe.pole_lat = atof(str.substr(50, 6).c_str());
				getline(ifile, str);
				coe.pole_lon= atof(str.substr(50, 6).c_str());
			}
			if ((str.substr(0, 18) == "PERIOD OF VALIDITY"))
			{
				getline(ifile, str);
				coe.jul.year = atoi(str.substr(49, 4).c_str());
				coe.jul.month = atoi(str.substr(54, 2).c_str());
				coe.jul.day = atoi(str.substr(57, 2).c_str());
				coe.jul.hour = atoi(str.substr(60, 2).c_str());
				coe.jul.minute = atoi(str.substr(63, 2).c_str());
				coe.jul.second = atof(str.substr(66, 2).c_str());
			}
			getline(ifile, str);
			if (str.empty() == true && ifile.eof()) break;
		}
		getline(ifile, str);
		coe.one_coe.clear();
		//数据体
		while (str.empty() != true)//每次更新数据最后会有一个空行
		{
			one_coe.degree = atoi(str.substr(2, 2).c_str());
			one_coe.order = atoi(str.substr(9, 3).c_str());
			one_coe.tec = atof(str.substr(18, 12).c_str());
			coe.one_coe.push_back(one_coe);
			getline(ifile, str);
		}
		sh_file->coeff.push_back(coe);
		//每次数据最后一行
		while (str.empty() == true){//防止最后出现多个空行
			getline(ifile, str);
			if (ifile.eof()) break;
		}
	} while (!ifile.eof());
	ifile.close();
	cout << "球谐函数系数读取完毕" <<endl;
}
//读取dcb文件
void read_dcb(string strd, pdcb dcbfile){
	cout << "开始读dcb文件 " << endl;
	dcb_f one_dcb;//一行数据
	ifstream ofile(strd, ios::in);
	if (!ofile){
		cout << "dcb文件打开错误" << endl;
		exit(0);
	}
	string str1;
	int n = 0;
	//cout << str1 << endl;
	while (!ofile.eof())
	{
		getline(ofile, str1);
		//cout << "dcb文件列数：" << str1.length() << endl;
		//cout << "str1.substr(1, 3)：" << str1.substr(1, 3) << endl;
		//cout << str1 << endl;
		if ((str1.substr(1, 3) == "DSB" || str1.substr(1, 3) == "DCB") && str1.length()>100)
		{
			//cout << "dcb文件行数：" << ++n << endl;
			one_dcb.prn = str1.substr(11, 3);
			one_dcb.station = str1.substr(15, 4);
			//cout << "str1.substr(11, 3):" << one_dcb.prn << "  str1.substr(15, 4):" << one_dcb.station << endl;
			one_dcb.obs1 = str1.substr(25, 3);
			one_dcb.obs2 = str1.substr(30, 3);
			one_dcb.bias = atof(str1.substr(82, 9).c_str());
			one_dcb.rms = atof(str1.substr(94, 9).c_str());
			dcbfile->dcb_val.push_back(one_dcb);
		}
		if (str1.substr(0, 8) == "%=ENDBIA")break;
	}
	ofile.close();
	cout << "读取dcb文件成功 " << endl;
	//cout << "dcb文件结果个数：" << dcbfile->dcb_val.size()<< endl;
	//for (int i = 0; i < dcbfile->dcb_val.size(); i++){
	//	cout << "prn" << dcbfile->dcb_val[i].prn << "station" << dcbfile->dcb_val[i].station << "obs1" << dcbfile->dcb_val[i].obs1 << "obs2" << dcbfile->dcb_val[i].obs2 << "dcb" << dcbfile->dcb_val[i].bias << endl;
//	}
}
void gpsttoutc(pgpst gt, ptc ut)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	gpsttojulianday(gt,ju);
	juliandaytoutc(ju, ut);
	free(ju);
}
void juliandaytobdt(pjulian ju, pbdt bt)
{
	double dt = ju->daynum + (ju->secondfrac + ju->secondnum) / _DAY_IN_SECOND - 2453736.500000;
	bt->week = int(dt / 7);
	bt->second = (dt - bt->week*7)*_DAY_IN_SECOND;
}
void utctobdt(ptc ut, pbdt bt)
{
	if (ut->year < 2006 || ut->month>12 || ut->month < 0 || ut->day>31 || ut->day < 0 || ut->hour>24 || ut->hour < 0 || ut->minute>60 || ut->minute < 0 || ut->second>60 || ut->second < 0)
	{
		cout << "时间大小有误"<<endl;
	}
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut,ju);
	juliandaytobdt(ju, bt);
	free(ju);
}

void utctogpst(ptc ut, pgpst gt)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut, ju);
	juliandaytogpst(ju, gt);
	free(ju);
}

void utctojulianday(ptc ut, pjulian ju)
{
	int		m;
	int		y;
	double	dhour;

	dhour = ut->hour + ut->minute / (double)_HOUR_IN_MINUTE
		+ ut->second / (double)_HOUR_IN_SECOND;

	if (ut->month <= 2) {
		y = ut->year - 1;
		m = ut->month + 12;
	}
	else {
		y = ut->year;
		m = ut->month;
	}

	ju->daynum = (long)(365.25*y) + (long)(30.6001*(m + 1))
		+ ut ->day + (long)(dhour / 24 + 1720981.5);
	ju->secondnum = ((ut->hour + 12) % _DAY_IN_HOUR)*_HOUR_IN_SECOND
		+ ut->minute*_MINUTE_IN_SECOND + (long)ut->second;
	ju->secondfrac = ut->second-(long)ut->second;
}

void juliandaytoutc(pjulian ju, ptc ut)
{
	int a, b, c, d, e;
	double JD;
	JD = ju->daynum + (ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;

	a = static_cast<int>(JD + 0.5);
	b = a + 1537;
	c = static_cast<int>((b - 122.1) / 365.25);
	d = static_cast<int>(365.25*c);
	e = static_cast<int>((b - d) / 30.6001);
	
	double day = b - d-(long)(30.6001*e) + JD + 0.5 - a;
	ut->day = int(day);
	ut->month = e - 1 - 12 * (int)(e / 14);
	ut->year = c - 4715 - (int)((7 + ut->month) / 10);

	ut->hour = int((day - ut->day)*24.0);
	ut->minute = (int)(((day - ut->day)*24.0 - ut->hour)*60.0);
	ut->second = ju->secondnum + ju->secondfrac - (int((ju->secondnum + ju->secondfrac)/60))*60.0;

}

void gpsttojulianday(pgpst gt, pjulian ju)
{
	double JD;
	JD = gt->weeknum * 7 + (gt->secondnum + gt->secondfrac) / _DAY_IN_SECOND + 2444244.5;
	ju->daynum= long(JD);
	
	ju->secondnum = long(gt->secondnum + (gt->weeknum * 7 + 2444244.5 - ju->daynum)*_DAY_IN_SECOND);
	ju->secondfrac = gt->secondfrac;
}


void juliandaytogpst(pjulian ju, pgpst gt)
{
	double JD;
	JD = ju->daynum +( ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;
	gt->weeknum = int((JD - 2444244.5) / 7);

	gt->secondnum = long((JD - 2444244.5 - gt->weeknum * 7)*_DAY_IN_SECOND);
	gt->secondfrac = ju->secondfrac;

}

//求儒略日差值
double deltjulianday(ptc u1, ptc u2)
{
	JULIANDAY j1, j2;
	utctojulianday(u1, &j1);
	utctojulianday(u2, &j2);
	double delt,d1,d2;
	d1 = j1.daynum + (j1.secondnum + j1.secondfrac) / _DAY_IN_SECOND;
	d2 = j2.daynum + (j2.secondnum + j2.secondfrac) / _DAY_IN_SECOND;
	//delt = (ju1->daynum - ju2->daynum)*_DAY_IN_SECOND + (ju1->secondnum - ju2->secondnum) + (ju1->secondfrac - ju2->secondfrac);
	delt = (d1 - d2)*_DAY_IN_SECOND;

	/*if (delt>302400)
		delt -= 604800;
	else if (delt<-302400)
		delt += 604800;
	else
		delt = delt;*/
	return delt;
}

//儒略日变换
void transjulian(pjulian ju1, double * dt, pjulian ju2)
{
	double JDold, JDnew;
	JDold = ju1->daynum + (ju1->secondnum + ju1->secondfrac) / _DAY_IN_SECOND;
	JDnew = JDold-*dt / _DAY_IN_SECOND;

	ju2->daynum = long(JDnew);
	ju2->secondnum = long((JDnew - long(JDnew))*_DAY_IN_SECOND);
	ju2->secondfrac = (JDnew - long(JDnew))*_DAY_IN_SECOND
		- long((JDnew - long(JDnew))*_DAY_IN_SECOND);
}






//坐标函数

//空间直角坐标系到大地坐标系
void xyztoblh(pxyz px, pblh pb)
{
	double pi = 4.0*atan(1.0);
	double E2 = 2.0*flattening - flattening * flattening;
	double E4 = E2*E2;
	double ALFA = (px->x*px->x + px->y*px->y + (1.0 - E2)*px->z*px->z) / (a*a);
	double BATA = (px->x*px->x + px->y*px->y - (1.0 - E2)*px->z*px->z) / (a*a);
	double Q = 1.0 + 13.50*E4*(ALFA*ALFA - BATA*BATA) / pow(ALFA - E4, 3);
	double A1 = -Q + sqrt(Q*Q - 1.0);
	double AL = (1.0 / 3.0)*log(-A1);
	AL = -exp(AL);
	double A2 = AL + 1.0 / AL;
	double A3 = AL - 1.0 / AL;
	double T23 = (ALFA + E4 / 2.0) / 3.0 - (ALFA - E4)*A2 / 12.0;
	double T32 = sqrt(T23 *T23 + ((ALFA - E4)*A3)*((ALFA - E4)*A3) / 48.0);
	double T1 = -E2*BATA / (4.0*T32);
	double DK = T1 + sqrt(T23 + T32) - (1.0 - E2 / 2.0);
	double EK = (1.0 + DK) / (1.0 + DK - E2);
	/*cout << "T1的值：" << T1 << endl;
	cout << "T32的值：" << T32 << endl;
	cout << "T23的值：" << T23 << endl;
	cout << "E2的值：" << E2 << endl;
	cout << "ALFA的值：" << ALFA << endl;
	cout << "AL的值：" << AL << endl;
	cout << "DK的值：" << DK << endl;
	cout << "EK的值：" << EK << endl;*/
	pb->height = (DK / (1.0 + DK))*sqrt(px->x*px->x + px->y*px->y + (EK*px->z)*(EK*px->z));
	double P = sqrt(px->x*px->x + px->y*px->y);
	pb->latitude = atan(EK*px->z / P);
	double COSFL = px->x / P;
	double SINFL = px->y / P;
	pb->longitude = asin(SINFL);
	if (SINFL > 0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;
	if (SINFL < 0.0&&COSFL>0.0) pb->longitude = 2.0*pi + pb->longitude;
	if (SINFL <0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;

	/*double e2;//第一偏心率的平方
	e2 = 2 * flattening - flattening*flattening;

	pb->longitude = atan(px->y / px->x);
	double W, N, N1 = 0, B, B1;
	B1 = atan(px->z / sqrt(px->x*px->x + px->y*px->y));
	while (1)
	{
		W = sqrt(1 - e2*sin(B1)*sin(B1));
		N1 = a / W;
		B = atan((px->z + N1*e2*sin(B1)) / sqrt(px->x*px->x + px->y*px->y));

		if (fabs(B - B1)<delta)
			break;
		else
			B1 = B;
	}

	pb->latitude = B;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));
	pb->height = sqrt(px->x*px->x + px->y*px->y) / cos(B) - N;*/
}

//大地坐标系到空间直角坐标系
void blhtoxyz(pblh pb, pxyz px)
{
	double e2;//第一偏心率的平方
	double N;//卯酉圈半径
	e2 = 2 * flattening - flattening*flattening;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));

	px->x = (N + pb->height)*cos(pb->latitude)*cos(pb->longitude);
	px->y = (N + pb->height)*cos(pb->latitude)*sin(pb->longitude);
	px->z = (N*(1 - e2) + pb->height)*sin(pb->latitude);
}

//笛卡尔坐标系到站心空间直角坐标系
void xyztoenu(pxyz pxcenter, pxyz px, penu pe)
{
	double dx, dy, dz;
	dx = px->x - pxcenter->x;
	dy = px->y - pxcenter->y;
	dz = px->z - pxcenter->z;

	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));

	xyztoblh(pxcenter, pd);

	pe->northing = -sin(pd->latitude)*cos(pd->longitude)*dx
		- sin(pd->latitude)*sin(pd->longitude)*dy
		+ cos(pd->latitude)*dz;
	pe->easting = -sin(pd->longitude)*dx
		+ cos(pd->longitude)*dy;
	pe->upping = cos(pd->latitude)*cos(pd->longitude)*dx
		+ cos(pd->latitude)*sin(pd->longitude)*dy
		+ sin(pd->latitude)*dz;
	free(pd);
}

//站心空间直角坐标系到笛卡尔坐标系
 void enutoxyz(pxyz pxcenter, penu pe, pxyz px)
{
	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));
	xyztoblh(pxcenter, pd);
	MatrixXd H(3, 3), DB(3, 1), DX(3, 1);
	DB(0, 0) = pe->northing;
	DB(1, 0) = pe->easting;
	DB(2, 0) = pe->upping;
	H(0, 0) = -sin(pd->latitude)*cos(pd->longitude);
	H(0, 1) = -sin(pd->latitude)*sin(pd->longitude);
	H(0, 2) = cos(pd->latitude);
	H(1, 0) = -sin(pd->longitude);
	H(1, 1) = cos(pd->longitude);
	H(1, 2) = 0;
	H(2, 0) = cos(pd->latitude)*cos(pd->longitude);
	H(2, 1) = cos(pd->latitude)*sin(pd->longitude);
	H(2, 2) = sin(pd->latitude);
	DX = (H.inverse())*DB;
	double dx, dy, dz;
	dx = DX(0,0 );
	dy = DX(1, 0);
	dz = DX(2, 0);
	px->x = pxcenter->x + dx;
	px->y = pxcenter->y + dy;
	px->z = pxcenter->z + dz;
	free(pd);
}

 //站心直角坐标系到站心极坐标系
 void enutoenupolar(penu pe, penupolar pep)
 {
	 pep->range = sqrt(pe->northing*pe->northing + pe->easting*pe->easting + pe->upping*pe->upping);
	 pep->azimuth = atan(pe->easting / pe->northing);
	 //atan2函数返回的范围为(-PI,PI]，当返回值大于零时，在1,2象限，小于零时在3,4象限
	 pep->azimuth = atan2(pe->easting / sqrt(pe->northing*pe->northing + pe->easting*pe->easting), pe->northing / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
	 if (pep->azimuth < 0.0)pep->azimuth += PI*2.0;
	 pep->elevation = atan(pe->upping / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
 }
 //求卫星高度角和方位角
 void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep){
	 penu pe;
	 pe = (penu)malloc(sizeof(ENU));
	 xyztoenu(pxcenter,  px, pe);
	 enutoenupolar(pe, pep);
	 //cout << "卫星方位角：" << pep->azimuth*180.0/PI << endl;
	 //cout << "卫星高度角：" << pep->elevation*180.0 / PI << endl;
	 //cout << "卫星距离：" << pep->range<< endl;
	 free(pe);
 }
 //求穿刺点地理坐标以及投影角
 void ipp_pos(pblh pb, penupolar pep,pblh pb1,double* MF){
	 double psi;//张角
	 psi = PI / 2.0 - pep->elevation - asin(ave_a / (ave_a + hion)*cos(pep->elevation));
	 pb1->latitude = asin(sin(pb->latitude)*cos(psi) + cos(pb->latitude)*sin(psi)*cos(pep->azimuth));
	 pb1->longitude = pb->longitude + atan(cos(pb->latitude)*sin(psi)*sin(pep->azimuth) / (cos(psi) - sin(pb->latitude)*sin(pb1->latitude)));
	 if (pb1->longitude < 0.0)pb1->longitude += PI*2.0;
	 if (pb1->longitude > PI*2.0)pb1->longitude -= PI*2.0;
	 pb1->height = hion;
	 *MF = 1.0 / sqrt(1.0 - (ave_a / (ave_a + hion)*cos(pep->elevation)*ave_a / (ave_a + hion)*cos(pep->elevation)));
 }
 //大地坐标转到日固地磁坐标
 void g2m(pblh pb,double mjd,double pole_lat,double pole_lon)
 {
	 double mag_lat, mag_lon,sun_lon;
	 //大地坐标转到地磁坐标
	 mag_lat = asin(sin(pole_lat)*sin(pb->latitude) + cos(pole_lat)*cos(pb->latitude)*cos(pb->longitude - pole_lon));
	 //mag_lon = atan(cos(pole_lat)*cos(pb->latitude)*sin(pb->longitude - pole_lon) / (sin(mag_lat)*sin(pole_lat) - sin(pb->latitude)));
	 //if (fabs(fabs(mag_lon) - PI / 2.0) < 1e-8)cout <<"接近临界值："<< mag_lon << endl;
	 mag_lon = atan2(cos(pb->latitude)*sin(pb->longitude - pole_lon) / cos(mag_lat), (sin(mag_lat)*sin(pole_lat) - sin(pb->latitude)) / (cos(mag_lat)*cos(pole_lat)));
	 pb->height = mag_lon;
	 //if (mag_lon < 0.0)mag_lon += PI*2.0;
	 //cout.precision(10);
	 sun_lon = PI*(1.0- 2.0 * (mjd - int(mjd)));//平太阳地理经度
	 //cout << "反正切：" << atan(sin(sun_lon - pole_lon) / sin(pole_lat) / cos(sun_lon - pole_lon)) << endl;
	 pb->latitude = mag_lat;
	 pb->sun_lon = atan2(sin(sun_lon - pole_lon) / sin(pole_lat), cos(sun_lon - pole_lon));
	 pb->longitude = mag_lon - atan2(sin(sun_lon - pole_lon)/ sin(pole_lat), cos(sun_lon - pole_lon));
	 if (pb->longitude < 0.0)pb->longitude += PI*2.0;
	 if (pb->longitude > PI*2.0)pb->longitude -= PI*2.0;
 }
 //由ionex文件插值计算任意点的vtec
 double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon){//lat和lon是大地经纬度
	 double p, q;//权
	 double v,v1, v2;//前后两个时间点的插值大小
	 //大地经纬度转为地心经纬度
	 double lat_geo, lon_geo;
	 lon_geo = lon;
	 if (lon_geo > 180.0)lon_geo -= 360.0;
	 lat_geo = atan((1 - flattening)*(1 - flattening)*tan(lat*PI/180.0))*180.0/PI;
	 if (fabs(lon_geo) > 180.0 || fabs(lat_geo) > 87.5 || deltjulianday(&vtime, &ionfile->vtec_map[0].vtime)<0.0){
		 cout << "穿刺点不在计算范围内" << endl;
		 return 0.0;
	 }
	 for (int i = 0; i < ionfile->vtec_map.size(); i++){
		 if (deltjulianday(&vtime, &ionfile->vtec_map[i].vtime) <= 0.0)//找到插值时间点后面的历元
		 {
			 for (int j = 0; j < 71; j++)
			 {
				 if (lat_geo > 87.5 - j*2.5)//找到插值点纬度后面的纬度节点
				 {
					 for (int k = 0; k < 73; k++)
					 {
						 if (lon_geo < -180.0 + k *5.0)//找到插值点经度后面的经度节点
						 {
							 p = fabs((lat_geo - (87.5-j*2.5))/2.5);
							 q = fabs((lon_geo -(-180.0+(k-1)*5.0)) / 5.0);
							//空间插值
							 //前一时间点
							 v1 = (1 - p)*(1 - q)*ionfile->vtec_map[i-1].tec_values[j][k - 1] +
								  (1 - p)*q*ionfile->vtec_map[i-1].tec_values[j][k] +
								   p*(1 - q)*ionfile->vtec_map[i-1].tec_values[j-1][k - 1] +
								   p*q*ionfile->vtec_map[i-1].tec_values[j-1][k];
							 //后一时间点
							 v2 = (1 - p)*(1 - q)*ionfile->vtec_map[i].tec_values[j][k - 1] +
								 (1 - p)*q*ionfile->vtec_map[i].tec_values[j][k] +
								 p*(1 - q)*ionfile->vtec_map[i].tec_values[j - 1][k - 1] +
								 p*q*ionfile->vtec_map[i].tec_values[j - 1][k];
							 //时间插值
							 v = (fabs(deltjulianday(&vtime, &ionfile->vtec_map[i - 1].vtime))*v2 +
								 fabs(deltjulianday(&vtime, &ionfile->vtec_map[i].vtime))*v1) / fabs(deltjulianday(&ionfile->vtec_map[i].vtime, &ionfile->vtec_map[i-1].vtime));
							 return v;
						 }
					 }
				 }
			 }
		 }
	 }
	 cout << "没能找到插值时间点" << endl;
	 return 0.0;
 }
 //计算归化勒让德多项式田谐项
 void geo_legendre(double lat, double**pg,int n)
 {
	 double sin_lat = sin(lat);
	 double tmp,tmp1,tmp2;
	 pg[0][0] = sqrt(3.0 * (1.0 - sin_lat*sin_lat));
	 //计算所有扇谐项
	for (int i = 1; i < n; i++){
		tmp = sqrt((2.0 * i + 3.0) / (2.0 * i + 2.0)*(1.0 - sin_lat*sin_lat));
		pg[i][i] = tmp*pg[i - 1][i - 1];
	}
	//计算除第一列以外的田谐项
	for (int j = 1; j < n - 1; j++){
		for (int i = j + 1; i < n; i++){
			tmp1 =sqrt((2.0 * i + 3.0)*(2.0 * i + 1.0)/(i+j+2)/(i-j))*sin_lat;
			tmp2 =sqrt((2.0 * i + 3.0)*(i + j + 1)*(i - j - 1) / (2.0 * i - 1.0) / (i + j + 2.0) / (i - j));
			pg[i][j] = tmp1*pg[i-1][j] - tmp2*pg[i-2][j];
		}
	}
	//计算田谐项的第一列
	int j = 0;
	for (int i = 1; i < n; i++)
	{
		tmp1 = sqrt(1.0*(2 * i + 3)*(2 * i + 1) / (i + j + 2) / (i - j))*sin_lat;
		tmp2 = sqrt(1.0*(2 * i + 3)*(i + j + 1)*(i - j - 1) / (2 * i - 1) / (i + j + 2) / (i - j));
		if (i == 1){ pg[i][j] = tmp1*pg[i - 1][j]; }
		else{ pg[i][j] = tmp1*pg[i - 1][j] - tmp2*pg[i - 2][j];}
	}
 
 }
 //计算归化勒让德多项式带谐项
 void zone_legendre(double lat, double *pz, int n){
	 double sin_lat = sin(lat);
	 pz[0] = 1.0;
	 pz[1] = sqrt(3.0)*sin_lat;
	 double tmp1, tmp2, tmp3;
	 for (int i = 2; i < n + 1; i++)
	 {
		 tmp1 = (2.0 - 1.0 / i)*sin_lat;
		 tmp2 = sqrt((2.0 * i - 1.0) / (2.0* i - 3.0))*(1.0 - 1.0 / i);
		 tmp3 = sqrt((2.0 * i + 1.0) / (2.0 *i - 1.0));
		 pz[i] = tmp3*(tmp1*(pz[i - 1]) - tmp2*pz[i - 2]);
	 }
 }
 //计算归化勒让德函数的值
 void legendre(double lat, double lon, double* legend,int n){//共(1+2L+1)(O+1)/2=(O+1)^2项
	 double**pg = new double*[n];
	 for (int i = 0; i < n; i++){
		 pg[i] = new double[n];
	 }
	 double* pz = new double[n + 1];
	 for (int i = 0; i < n; i++){//初始化
		 pz[i] = 0.0;
		 for (int j = 0; j < n; j++){
			 pg[i][j] =0.0;
		 }
	 }
	 pz[n] = 0.0;
	 //double pg[O][O] = { 0.0 };
	// double pz[O+1] = { 0.0 };
	 geo_legendre(lat, pg,n);
	 zone_legendre(lat, pz,n);
	 cout.setf(ios_base::fixed, ios_base::floatfield);
	 cout.precision(10);
	 for (int i = 0; i < n + 1; i++)
	 {
		 legend[i*i] = pz[i];
		 //cout << "pz[" << i << "]:" << pz[i] << endl;
		 //if (i == O){ break; }
		 for (int j = 0; j < i; j++)
		 {
			legend[i*i+ j*2+1] = pg[i-1][j] * cos((j+1)*lon);
			legend[i*i+ j*2+2] = pg[i-1][j] * sin((j+1)*lon);
			if (legend[i*i + j * 2 + 1] == 0.0 || legend[i*i + j * 2 + 2] == 0.0){ cout << "                                               勒让德值c失败" << endl;}
			//cout << "pg[" << i-1 << "]" << "[" << j << "]:" << pg[i-1][j] << endl;
		 }
	 }
	 /*cout << "纬度：:" <<lat << endl;
	 for (int i = 0; i < O + 1; i++){
		 cout << "legend[" << i << "]:" << legend[i] << endl;
	 }*/
	 for (int i = 0; i < n; i++){
		 delete[] pg[i];
	 }
	 delete[] pg;
	 delete pz;
 }
 //由球谐函数文件插值计算任意点的vtec
 double sh_vtec(psh sh_file, UTC vtime, double lat, double lon,int n){//lat，lon是大地坐标
	 double *legend = new double[(n + 1)*(n + 1)];
	 double mjd;
	 double  dt1, dt2;//dt1是所求时间和前一历元的差值比值，dt2是所求时间和后一历元的差值比值
	 JULIANDAY jul;
	 one_coeff coe1,coe2;
	 BLH ipp;
	 //cout << "观测时间：" << vtime.hour << ":" << vtime.minute << ":" << vtime.second << endl;
	 ipp.latitude = lat; ipp.longitude = lon; ipp.height = hion;
	 double vtec = 0.0;
	 utctojulianday(&vtime, &jul);//把计算穿刺点时的utc转到儒略日
	 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5;//儒略日转到约化儒略日
	 for (int i = 0; i < sh_file->coeff.size(); i++)
	 {
		 if (deltjulianday(&vtime, &sh_file->coeff[i].jul) <= 0.0 || i == (sh_file->coeff.size()-1))//找到插值时间点后面的历元或者最后一个
		 {//暂时没有使用时间插值
			 if (i == 0){ coe1 = sh_file->coeff[i]; coe2 = sh_file->coeff[i]; }//第一个历元
			 else if (i == (sh_file->coeff.size() - 1) && deltjulianday(&vtime, &sh_file->coeff[i].jul) > 0.0){ coe1 = sh_file->coeff[i]; coe2 = sh_file->coeff[i]; }//最后一个历元的数据且其时间不是后一天的零点
			 else{ coe1 = sh_file->coeff[i - 1]; coe2 = sh_file->coeff[i-1]; }//球谐函数系数的有效范围一般为其后1或者2个小时，所以取观测时间前面的历元
			 //dt1 = deltjulianday(&vtime, &coe1.jul) / deltjulianday(&coe2.jul, &coe1.jul); dt2 = 1.0 - dt1;
			 dt1 = 1.0; dt2 = 1.0 - dt1;
			 cout.setf(ios::fixed);
			 cout.precision(10);
			 g2m(&ipp, mjd, coe1.pole_lat*PI / 180.0, coe1.pole_lon*PI / 180.0);//穿刺点由大地坐标转到日固地磁坐标
			 legendre(ipp.latitude, ipp.longitude, legend,n);//由穿刺点计算勒让德多项式的值
			// cout << "code球谐函数系数"<<endl;
			 //cout.precision(4);
			 for (int j = 0; j <coe1.one_coe.size(); j++)
			 {
				 vtec += (coe1.one_coe[j].tec*dt2 + coe2.one_coe[j].tec*dt1) * legend[j];
				// cout <<std::left<<setw(8)<< coe1.one_coe[j].tec;
			 }
			 //cout << endl;
			 return vtec;
			 delete legend;
		 }
	 }
	 cout << "没有找到有效的球谐函数系数" << endl; 
	 delete legend;
	 return 0.0;
 }

 //先预报球谐函数系数，然后计算vtec
 double predict_vtec(UTC vtime, double lat, double lon,double pole_lat,double pole_lon, int n, double coe[][25]) 
 {
	 double vtec = 0.0;
	 double *legend = new double[(n + 1)*(n + 1)];
	 double *sh_coe = new double[(n + 1)*(n + 1)]{ 0.0 };
	 double L[25];
	 double period[12] = { 0.33, 0.5, 1.0, 14.6, 27.0, 121.6, 182.51, 365.25, 1007.18, 1342.9, 2014.35, 4028.71 };//12个主要周期
	 double mjd;
	 BLH ipp;
	 JULIANDAY jul;
	 ipp.latitude = lat; ipp.longitude = lon; ipp.height = hion;
	 utctojulianday(&vtime, &jul);//把计算穿刺点时的utc转到儒略日
	 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5;//儒略日转到约化儒略日
	 g2m(&ipp, mjd,pole_lat*PI / 180.0,pole_lon*PI / 180.0);//穿刺点由大地坐标转到日固地磁坐标
	 legendre(ipp.latitude, ipp.longitude, legend, n);//由穿刺点计算勒让德多项式的值
	 L[0] = 1.0;
	 for (int i = 0; i < 12; i++)
	 {
		 L[i*2+1] = cos(2*PI/period[i]*mjd);
		 L[i*2+2] = sin(2*PI/period[i]*mjd);
	 }
	 //cout << "预报球谐函数系数  "<<endl;
	// cout.precision(4);
	 for (int i = 0; i < (n + 1)*(n + 1); i++) 
	 {
		 for (int j = 0; j < 25; j++) 
		 {
			 sh_coe[i] += L[j] * coe[i][j];
		 }
		 vtec += sh_coe[i] * legend[i];
		// cout << std::left << setw(8) << sh_coe[i];
	 }
	 //cout << endl;
	 delete legend;
	 delete sh_coe;
	 return vtec;
 }
 //矩阵转置
 void Transposition(double**L, int n)
 {
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < i; j++)
			 swap(L[i][j], L[j][i]);
	 }
 }
 //矩阵乘法
 void Multi(double**A, double**B, int n)//AXB->B
 {
	 double **C = new double*[n];
	 for (int i = 0; i < n; i++)
		 C[i] = new double[n];
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < n; j++)
		 {
			 C[i][j] = 0;
			 for (int k = 0; k < n; k++)
				 C[i][j] += A[i][k] * B[k][j];
		 }
	 }
	 for (int i = 0; i < n; i++)
	 {
		 for (int j = 0; j < n; j++)
			 B[i][j] = C[i][j];
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] C[i];
		 C[i] = NULL;
	 }
	 delete C;
	 C = NULL;
 }
 //Cholesky分解方法
 void chol_rf(double**A, double**L, double**d, int n)
 {
	 double s1 = 0.0;
	 double **g = new double*[n]; //开辟行  
	 for (int i = 0; i < n; i++)
		 g[i] = new double[n]; //开辟列  
	 for (int i = 0; i < n; i++){//初始化
		 for (int j = 0; j < n; j++){
			 g[i][j] = 0.0;
		 }
	 }
	 d[0][0] = A[0][0];
	 for (int i = 1; i < n; i++){
		 for (int j = 0; j < i; j++){
			 s1 = 0.0;
			 for (int k = 0; k < j; k++){
				 s1 += g[i][k] * L[j][k];
			 }
			 g[i][j] = A[i][j] - s1;
		 }
		 for (int j = 0; j < i; j++){
			 L[i][j] = g[i][j] / d[j][j];
		 }
		 s1 = 0.0;
		 for (int k = 0; k < i; k++){
			 s1 += g[i][k] * L[i][k];
		 }
		 d[i][i] = A[i][i] - s1;
	 }
	 for (int i = 0; i < n; i++){
		 L[i][i] = 1.0;
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] g[i];
	 }
	 delete[] g;
 }
 /*
 void chol_rf(double**A, double**L, double**d, int n)
 {
	 double s1 = 0.0, s2 = 0.0;;
	 for (int j = 0; j < n; j++)
	 {
		 s1 = 0.0;
		 for (int k = 0; k < j; k++)
		 {
			 s1 += L[j][k] * L[j][k] * d[k][k];
		 }
		 d[j][j] = A[j][j] - s1;
		 for (int i = j + 1; i < n; i++)
		 {
			 s2 = 0.0;
			 for (int k = 0; k < j; k++)
			 {
				 s2 += L[i][k] * L[j][k] * d[k][k];
			 }
			 L[i][j] = (A[i][j] - s2) / d[j][j];
		 }
		 L[j][j] = 1.0;
	 }
 }*/
 //Cholesky分解方法解方程组
 void chol_eq(double**A, double* b, double* x,int n)
 {
	 //cout << "Cholesky分解方法解方程组" << endl;
	 double *y = new double[n];
	 double **L = new double*[n]; //开辟行  
	 double **d = new double*[n]; //开辟行
	 for (int i = 0; i < n; i++)
	{	L[i] = new double[n]; //开辟列 
		d[i] = new double[n]; //开辟列 
	}
	 for (int i = 0; i < n; i++){//初始化
		 y[i] = 0.0;
		 for (int j = 0; j < n; j++){
			 L[i][j] = 0.0;
			 d[i][j] = 0.0;
		 }
	 }
	 chol_rf(A, L, d,n);//Cholesky分解
	 y[0] = b[0];
	 double s1;
	 for (int i = 1; i < n; i++){
		 s1 = 0.0;
		 for (int k = 0; k < i; k++){
			 s1 += L[i][k] * y[k];
		 }
		 y[i] = b[i] - s1;
	 }
	 x[n - 1] = y[n - 1] / d[n - 1][n - 1];
	 for (int i = n - 2; i >= 0; i--){
		 s1 = 0.0;
		 for (int k = i+1; k < n; k++){
			 s1 += L[k][i] * x[k];
		 }
		 x[i] = y[i] / d[i][i] - s1;
	 }
	 for (int i = 0; i < n; i++)
	 {
		 delete[] L[i];
		 delete[] d[i];
	 }
	 delete y;
	 delete[] L;
	 delete[] d;
 }
 /*
 void chol_eq(double**A, double* b, double* x,int n)
 {
	 //cout << "Cholesky分解方法解方程组" << endl;
	 double **L = new double*[n]; //开辟行  
	 double **d = new double*[n]; //开辟行
	 for (int i = 0; i < n; i++)
	 {
		 L[i] = new double[n]; //开辟列 
		 d[i] = new double[n]; //开辟列 
	 }
	 for (int i = 0; i < n; i++){//初始化
		 for (int j = 0; j < n; j++){
			 L[i][j] = 0.0;
			 d[i][j] = 0.0;
		 }
	 }
	 chol_rf(A, L, d,n);//Cholesky分解
	 
 for (int k = 0; k < n; k++)
 {
	 for (int i = 0; i < k; i++)
		 b[k] -= b[i] * L[k][i];
	 b[k] /= L[k][k];
 }
 
 Transposition(L, n);//L^T
 Multi(d, L, n);//D X L^T
 for (int k = n - 1; k >= 0; k--)
 {
	 for (int i = k + 1; i < n; i++)
		 b[k] -= b[i] * L[k][i];
	 b[k] /= L[k][k];
 }
 for (int k = 0; k < n; k++)
 {
	 x[k] = b[k];
 }
 for (int i = 0; i < n; i++)
 {
	 delete[] L[i];
	 delete[] d[i];
 }
 delete[] L;
 delete[] d;
 }
 */
 //寻找chazh星历
 bool find_sp3_ephem(string prn, ptc ut, psp3 sp3file, s_sp3_ephe sp3[])//n：拉格朗日插值法阶数，sp3：返回的差值节点值
 {
	 int num;//差值点前面的节点数
	 if (NN % 2 == 0)num = NN / 2;//差值阶数为偶数
	 else num = (NN + 1) / 2;
	// JULIANDAY jld1, jld2;//jlld1卫星信号发射时间，jlld2卫星信星历时间
	 /*cout << "卫星信号发射时间：" << ut->year<< ":" << ut->month << ":" << ut->day;
	 cout << ":" << ut->hour << ":" << ut->minute << ":" << ut->second << endl; */
	 int kk = 32;
	 for (int i = 0; i < kk; i++)
	 {
		 if (prn.substr(0, 3) == sp3file->all_ephem[i].prn.substr(0, 3))
		 {
			 //cout << "找到对应的卫星号" << endl;
			 // satn = true;
			 int w = sp3file->all_ephem[i].sate_ephem.size();
			 for (int j = 0; j < w; j++)
			 {
				 double delta;
				 delta = deltjulianday(ut, &sp3file->all_ephem[i].sate_ephem[j].utime_n);
				 if (delta > 0.0) continue;//找到差值点前面的差指点时间
				 else{
					 //cout << "找到历元" << endl;
					 if (j == 0) return false;//插指点在所有插值节点之外（之前）
					 else if (j  < num){//j刚好是插值时间点前面的可用插值点数量，从第一个插值点开始取前n+1个节点用来插值
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[k];
							// cout << "前" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else if (w - j < NN + 1 - num){//插值点在最后面，后面节点数不足够
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[w-1 - NN + k];
							 //cout << "后" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else{
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[j - num + k];//j是第num+1个
							// cout << "中" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 return true;
				 }
			 }
		 }
	 }
	 // cout
	 return false;
 }
 
 //精密星历计算卫星坐标
 bool cal_sp3_sate_coor(string prn, ptc ut, psp3 sp3file, pxyz coor){
	 coor->x = 0.0;
	 coor->y = 0.0;
	 coor->z = 0.0;
	 s_sp3_ephe sp3[NN + 1];
	 double s = 1.0;
	 if(find_sp3_ephem(prn, ut, sp3file, sp3))
	 {
		 //cout <<"找到星历"<< endl;
		 for (int i = 0; i < NN + 1; i++)
		 {//拉格朗日法求解卫星坐标
			 s = 1.0;
			 for (int j = 0; j < NN + 1; j++)
			 {
				 if (i == j)continue;
				 //cout << "x-xj:" << deltjulianday(ut, &sp3[j].utime_n) << "xi-xj:" << deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n) << endl;
				 s = s*deltjulianday(ut, &sp3[j].utime_n) / deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n);
			 }
			 coor->x += s*sp3[i].x;
			 coor->y += s*sp3[i].y;
			 coor->z += s*sp3[i].z;
		 }
		 return true;
	 }
	 else{ 
		 return false; 
	 }

 }
 //从dcb文件中找到对应的dcb值
 bool find_dcb(bool b,string prn, string station, pdcb dcbfile, double* dcb_val){//b==0计算卫星dcb，b!=0计算测站dcb
	 //cout << "dcb文件大小：" << dcbfile->dcb_val.size() << endl;
	 for (int i = 0; i < dcbfile->dcb_val.size(); i++)
	 {
		 //cout << "卫星prn:" << dcbfile->dcb_val[i].prn << "  测站:" << dcbfile->dcb_val[i].station << "  obs1:" << dcbfile->dcb_val[i].obs1 << "  obs2:" << dcbfile->dcb_val[i].obs2 << endl;
		 //cout << "prn:" << prn << "  station:" << station <<endl;
		 if (((b==0&&dcbfile->dcb_val[i].prn == prn)||(b!=0&&dcbfile->dcb_val[i].station == station))&&dcbfile->dcb_val[i].obs1 == "C1C"&&dcbfile->dcb_val[i].obs2 == "C2W")
		 {
			 *dcb_val = dcbfile->dcb_val[i].bias*1e-9*c;//返回距离值
			 return true;
		 }
	 }
	 return false;
 }
 //计算c_a0
 double c_a0(double**ctable, UTC ut, double lat, double lon){
	 double legend[36] = { 0.0 };
	 double beta[17] = { 0.0 };//17个通过周期项计算的系数
	 double b[17] = { 0.0 };//17个和穿刺点有关的勒让德系数值
	 double a0 = 0.0;
	 JULIANDAY jul;
	 double mjd;
	 legendre(lat, lon, legend, 5);
	 //mjd 对应当天约化儒略日的奇数整点时刻
	 ut.minute = 0; ut.second = 0.0;
	 if (ut.hour % 2 == 0){ ut.hour += 1; }
	 utctojulianday(&ut, &jul);//把计算穿刺点时的utc转到儒略日
	 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5; 
	 double per[12] = { 1.0, 0.5, 0.33, 14.6, 27.0, 121.6, 182.51, 365.25, 4028.71, 2014.35, 1342.9, 1007.18 };
	 for (int j = 0; j < 17; j++){
		 beta[j] = ctable[0][j];
		 for (int i = 0; i < 12; i++){
			 beta[j] += ctable[2 * i + 1][j] * cos(2 * PI / per[i] * mjd)+ctable[2 * i + 2][j] * sin(2 * PI / per[i] * mjd);
		 }
	 }
	 for (int j = 0; j < 12; j++){
		 b[j] = legend[9 + j];//3度项和4度项前5项
	 }
	 for (int j = 12; j < 17; j++){
		 b[j] = legend[13 + j];//5度项前5项
	 }
	 for (int j = 0; j < 17; j++){
		 a0 += beta[j] * b[j];
	 }
	 return a0;
 }
 //结果输出
 void putresult(pvt pt)
 {
	 ofstream outfile("E:\\GNSS\\结果\\sh_vtec.txt", ios::out);
	 if (!outfile)
	 {
		 cerr << "文件创建失败" << endl;
		 abort();
	 }
	 cout << "开始写文件" << endl;
	 for (unsigned int i = 0; i < pt->allresult.size(); i++)
	 {
		 outfile.setf(ios_base::fixed, ios_base::floatfield);
		 outfile.precision(8);
		 outfile <<setw(4) << pt->allresult[i].rtime.year << setw(3) << pt->allresult[i].rtime.month << setw(3) << pt->allresult[i].rtime.day << 
			 setw(3) << pt->allresult[i].rtime.hour << setw(3) << pt->allresult[i].rtime.minute << setw(15) << pt->allresult[i].rtime.second<<endl;
		 outfile << " 残差                                                                                      个数：" << pt->allresult[i].residul.size()<< endl;
		 for (unsigned int j = 0; j < pt->allresult[i].residul.size(); j++){
			 outfile << std::right << setw(12) << pt->allresult[i].residul[j];
			 if ((j + 1) % 10 == 0)outfile << endl;
		 }
		 outfile << endl;
		 outfile << "COEFFICIENTS" << endl;
		 outfile << "DEGREE   ORDER        VALUE (TECU)        VALUE (TECU)" << endl;
		 for (unsigned int k = 0; k < O + 1; k++){
			 outfile << std::right << setw(4) << k << std::right << setw(9) << 0 << std::right << setw(20) << pt->allresult[i].coe[k*k] << std::right << setw(20) << pt->allresult[i].coe1[k*k] << endl;
			 for (int l = 1; l <= k; l++){
				 outfile << std::right << setw(4) << k << std::right << setw(9) << l << std::right << setw(20) << pt->allresult[i].coe[k*k + l * 2 - 1] << std::right << setw(20) << pt->allresult[i].coe1[k*k + l * 2 - 1] << endl;
				 outfile << std::right << setw(4) << k << std::right << setw(9) << -l << std::right << setw(20) << pt->allresult[i].coe[k*k + l * 2] << std::right << setw(20) << pt->allresult[i].coe1[k*k + l * 2] << endl;
			 }
		 }
	 }
	 outfile.close();
	 cout << "文件写入完毕" << endl;
 }
 //由观测值计算vtec
 double obs_vtec(pobs obsfile, UTC obs_t,int n, string prn, psp3 sp3file, pdcb dcbfile, pblh ipp)
 {
	 UTC day_first;//一天的0时
	 day_first = obs_t;
	 day_first.hour = 0; day_first.minute = 0; day_first.second = 0.0;
	 int ns=0;//一天之中的第n个观测历元,ns是要计算的观测历元之前的起始时间
	 int flag;//判断上个历元是否依然观测到要计算的卫星
	 vector<int>prn_pos;//第i个历元卫星prn存放的位置
	 int p;//位置
	 int qq = 0;//num:穿刺点数量,qq:连续观测一颗卫星的观测历元数,cycle_slip:周跳个数
	 XYZ sta_coor;//测站wgs84坐标
	 BLH sta_coor_blh;//测站和穿刺点的大地坐标
	 pxyz sate_coor = new XYZ;//卫星坐标
	 penupolar pep = new ENUPOLAR;//卫星方位角和高度角
	 double* mf = new double;//投影函数
	 double dt2 = 0.0;//dt1：卫星信号传播时间,dt2：同一颗卫星相邻历元的时间差
	 double w = 1.0;//权
	 double mf0 = 1.0;//上一步的mf
	 double dcb_sat = 0.0, dcb_sta = 0.0;//测站和卫星dcb
	 //周跳探测
	 double n_mw1 = 0.0, n_mw2 = 0.0;//MW组合
	 double n0 = 0.0, n1 = 0.0;//均值
	 double sigma0 = 0.0, sigma1 = 0.0;//标准差
	 int cycle_slip = 0;//周跳次数
	 double vtec1 = 0.0, vtec0 = 0.0;//上一步计算得到的vtec0和当前计算的vtec1
	 if (!find_dcb(0, prn, obsfile->obsheaddata.station, dcbfile, &dcb_sat) || !find_dcb(1, prn, obsfile->obsheaddata.station, dcbfile, &dcb_sta))
	 {
		 cout << "没有找到卫星" << prn << "或者测站dcb" << endl;
		 delete sate_coor;
		 delete pep;
		 delete mf;
		 return 0.0;
	 }
	 sta_coor = obsfile->obsheaddata.approx_coordinate;
	 xyztoblh(&sta_coor, &sta_coor_blh);
	 //计算起始观测时间
	 for (int i = n; i >= 0; i--)
	 {
		 flag = 0;
		 for (int j = 0; j < obsfile->obsdata[i].prn_list.size(); j++)
		 {
			 if (prn == obsfile->obsdata[i].prn_list[j])
			 {
				 flag = 1;
				 prn_pos.push_back(j);
				 break;
			 }
		 }
		 if (flag == 0)//不再有要计算的卫星的观测数据就结束了
		 {
			 ns = i+1;//卫星开始被观测到的历元
			 break;
		 }
		// if (i == 0)ns = 0;//起始时间是一天的开始,初始化ns=0即可
	 }
	;
	 //一颗卫星所有历元按照相位平滑伪距的方法求vtec
	 //int c = 0;
	 for (int j = ns; j <= n; j++)
	 {
		 //cout << "j:" << j << endl;
		// cout << "qq:" << qq<< endl;
		 p = n - j;//第j个历元卫星的观测顺序号对应的观测顺序数组的顺序
		 if (!cal_sp3_sate_coor(prn, &obsfile->obsdata[j].utime_o, sp3file, sate_coor))//求卫星坐标
		 {
			 //cout << "卫星坐标计算失败" << endl;
			 continue;
		 }//卫星坐标计算失败
		 sate_azi_ele(&sta_coor, sate_coor, pep);//求卫星方位角和高度角
		 if (pep->elevation*180.0 / PI < 15.0){ continue; }//卫星高度角要大于15度
		 ipp_pos(&sta_coor_blh, pep, ipp, mf);//求穿刺点纬度经度,单位弧度
		// if (prn == "G04"){
		//	 cout << "p2:" << obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2 << " p1:" << obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 << " dcb_sta:" << dcb_sta << " dcb_sat:" << dcb_sat << " mf:" << *mf << endl;
		 //}
		 if (j != ns)
		 {//第一步不计算dt2
			 dt2 = fabs(deltjulianday(&obsfile->obsdata[j].utime_o, &obsfile->obsdata[j - 1].utime_o));
		 }
		 if (qq == 0 || dt2 > 300.0)//qq=0意味着上一步发现周跳或者第一步，相邻历元时间间隔大于300s，证明卫星落下又升起，需要选为时间起点计算vteco文件个数
		 {
			 vtec1 = 9.52437*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 + dcb_sta + dcb_sat) / (*mf);
			 w = 1.0;
			 qq = 1;
			 //cout << "vtec1:" << vtec1 << endl;
			 n_mw1 = obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1*f1 / c + obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2*f2 / c);
			 n1 = n_mw1;//平均值，第一步直接等于计算值
			 sigma1 = 0.0;//方差，第一步方差为0
			 //continue;
		 }
		 else{
			 n0 = n1;//上一步的平均值
			 sigma0 = sigma1;//上一步的方差
			 n_mw1 = obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1*f1 / c + obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2*f2 / c);
			 if (j<obsfile->obsdata.size()-1){
				 if (j == n)
				 { //到观测历元的时候需要用到后一历元的数据来计算n_mw2
					 n_mw2 = n_mw1;//防止后一历元没有对应的卫星
					 for (int k = 0; k < obsfile->obsdata[n+1].prn_list.size();k++)
					 {
						 if (prn == obsfile->obsdata[n+1].prn_list[k])
						 {
							 n_mw2 = obsfile->obsdata[j + 1].one_obs_data[k].l1 - obsfile->obsdata[j + 1].one_obs_data[k].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j + 1].one_obs_data[k].p1*f1 / c + obsfile->obsdata[j + 1].one_obs_data[k].p2*f2 / c);
							 break;
						 }
					 }
				 }
				 else{
					 n_mw2 = obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].l1 - obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].l2 - (f1 - f2) / (f1 + f2)*(obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].p1*f1 / c + obsfile->obsdata[j + 1].one_obs_data[prn_pos[p - 1]].p2*f2 / c);
				 }
			 }
			 else{ n_mw2 = n_mw1; }//最后一个观测值
			 /*if (qq>2 && fabs(n_mw1 - n0) >= 4.0*sigma0&&fabs(n_mw2 - n_mw1) < 1.0){//发现周跳，最后一个历元n_mw2 == n_mw1，因为没有重新赋值给n_mw2
				 cycle_slip++;
				 qq = 0;
				 cout<<"发现周跳"<<prn<<endl;
				 continue;
			 }*/
			 //下一步迭代值
			 n1 = n0 + 1.0 / qq*(n_mw1 - n0);
			 sigma1 = sqrt(sigma0*sigma0 + 1.0 / qq*((n_mw1 - n0)*(n_mw1 - n0) - sigma0*sigma0));
			 w = 1.0 / qq;
			 vtec1 = 9.52437*((obsfile->obsdata[j].one_obs_data[prn_pos[p]].p2 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].p1 + dcb_sta + dcb_sat)*w + (1.0 - w)*(vtec0 / 9.52437*mf0 +
				 ((obsfile->obsdata[j].one_obs_data[prn_pos[p]].l1*c / f1 - obsfile->obsdata[j].one_obs_data[prn_pos[p]].l2*c / f2) - (obsfile->obsdata[j - 1].one_obs_data[prn_pos[p + 1]].l1*c / f1 - obsfile->obsdata[j - 1].one_obs_data[prn_pos[p + 1]].l2*c / f2)))) / (*mf);
		 }
		 if (vtec1<0.0){
			 qq = 0;
			 continue;
		 }
		 //num++;
		 qq++;
		 mf0 = *mf;
		 vtec0 = vtec1;
	 }
	 delete sate_coor;
	 delete pep;
	 delete mf;
	 return vtec1;
 }
 //计算VTEC
 void vtec(vector<obs>obsfile,double** ctable, double coe[][25],psp3 sp3file, pio ionfile, pdcb dcbfile, psh sh_file)
 {
	 int n = obsfile.size();
	 int m = obsfile[0].obsdata.size();
	 BLH ipp;
	 double mf;
	 double vtec_obs,vtec_sh,vtec_ion,vtec_predict;
	 XYZ sta_coor;//测站wgs84坐标
	 BLH sta_coor_blh;//测站和穿刺点的大地坐标
	 XYZ sate_coor;//卫星坐标
	 ENUPOLAR pep;
	 JULIANDAY jul;
	 double mjd;
	 double pole_lat = 80.33, pole_lon = -72.67;
	 int num = 0;
	 MatrixXd  A1A((O + 1)*(O + 1), (O + 1)*(O + 1)), A2A((O + 1)*(O + 1), (O + 1)*(O + 1)), AA((O + 1)*(O + 1), (O + 1)*(O + 1)), A1B((O + 1)*(O + 1),1), A2B((O + 1)*(O + 1), 1), BB((O + 1)*(O + 1), 1), X((O + 1)*(O + 1), 1);
	 //MatrixXd  atb((O + 1)*(O + 1), 1), ata((O + 1)*(O + 1), (O + 1)*(O + 1)), mx((O + 1)*(O + 1), 1), mx1((O + 1)*(O + 1), 1), inv_ata((O + 1)*(O + 1), (O + 1)*(O + 1));
	 //vector<one_legendre>A;
	 //vector<double>b;
	 //one_legendre a;
	 double sigma1 = 1.0, sigma2 = 1.0;//观测矩阵和预报矩阵单位权中误差
	 double le[(15 + 1)*(15 + 1)];//勒让德函数
	 double a0;//球谐函数系数高阶项的vtec
	 int flag = 0;//判断测站是否在中国境内
	 string s_name;
	 ofstream outfile("E:\\GNSS\\数据\\vtec.txt", ios::out);
	 ofstream coe_update("E:\\GNSS\\数据\\coe.txt", ios::out);
	 outfile.setf(ios_base::fixed, ios_base::floatfield);
	 //outfile.precision(5);
	 coe_update.setf(ios_base::fixed, ios_base::floatfield);
	 //coe_update.precision(12);
	 if (!outfile || !coe_update)
	 {
		 cerr << "文件创建失败" << endl;
		 abort();
	 }
	 //矩阵初始化
	 for (int w = 0;w < (O + 1)*(O + 1);w++)
	 {
		 A1B(w, 0) = 0.0;
		 A2B(w, 0) = 0.0;
		 for (int q = 0;q < (O + 1)*(O + 1);q++)
		 {
			 A1A(w, q) = 0.0;
			 A2A(w, q) = 0.0;
		 }
	 }
	 cout << "观测历元个数：" << m << endl;
	 for (int i = 0; i < m; i++)//一天之中的观测次数，30s一次，共2880次
	 {
		 cout << "历元：" << i << endl;
		 for (int j = 0; j < n; j++)//n：测站数量，所有测站在同一个时间观测的数据组成一个方程组来解球谐函数系数
		 {
			 s_name = obsfile[j].obsheaddata.station;
			 //国内测站
			 if (s_name== "BJFS"||s_name== "BJNM"||s_name=="CHAN"||s_name== "HKSL"||s_name=="HKWS"||s_name=="JFNG"||s_name =="SHAO"||s_name =="URUM"||s_name == "WUH2") 
			 {
				 flag = 1;
			 }
			 if (fabs(deltjulianday(&obsfile[0].obsdata[i].utime_o, &obsfile[j].obsdata[i].utime_o)) > 0.1)
			 {
				 cout << "                                                                      时间有差异" << endl;
				 continue; 
			 }//30s一个观测间隔，一般来说o文件是30s会有一次观测数据
			 sta_coor = obsfile[j].obsheaddata.approx_coordinate;
			 xyztoblh(&sta_coor, &sta_coor_blh);
			 //cout << "测站名称" << obsfile[j].obsheaddata.station << endl;
			 //cout << "时间：" << obsfile[j].obsdata[i].utime_o.hour << ":" << obsfile[j].obsdata[i].utime_o.minute << ":" << obsfile[0].obsdata[i].utime_o.second << endl;
			 for (int k = 0; k < obsfile[j].obsdata[i].prn_list.size(); k++)//一个测站在一个时间观测到的所有卫星
			 {
				 if (!cal_sp3_sate_coor(obsfile[j].obsdata[i].prn_list[k], &obsfile[j].obsdata[i].utime_o, sp3file, &sate_coor))//求卫星坐标
				 {
					 continue;
				 }//卫星坐标计算失败
				 sate_azi_ele(&sta_coor, &sate_coor, &pep);//求卫星方位角和高度角
				 if (pep.elevation*180.0 / PI < 15.0){ continue; }//卫星高度角要大于15度
				 ipp_pos(&sta_coor_blh, &pep, &ipp, &mf);//求穿刺点纬度经度,单位弧度
				 vtec_obs = obs_vtec(&obsfile[j], obsfile[j].obsdata[i].utime_o, i, obsfile[j].obsdata[i].prn_list[k], sp3file, dcbfile, &ipp);
				 vtec_ion = ionex_vtec(ionfile, obsfile[j].obsdata[i].utime_o, ipp.latitude*180.0 / PI, ipp.longitude *180.0 / PI) / 10.0;
				 vtec_sh = sh_vtec(sh_file, obsfile[j].obsdata[i].utime_o, ipp.latitude, ipp.longitude,15);
				 vtec_predict = predict_vtec(obsfile[j].obsdata[i].utime_o, ipp.latitude, ipp.longitude, pole_lat, pole_lon, 15, coe);
				 a0 = c_a0(ctable, obsfile[j].obsdata[i].utime_o, ipp.latitude, ipp.longitude);
				 utctojulianday(&obsfile[j].obsdata[i].utime_o, &jul);//把计算穿刺点时的utc转到儒略日
				 mjd = jul.daynum + (jul.secondnum + jul.secondfrac) / 86400.0 - 2400000.5;//儒略日转到约化儒略日
				 g2m(&ipp, mjd, pole_lat*PI / 180.0, pole_lon*PI / 180.0);//穿刺点由大地坐标转到日固地磁坐标
				 legendre(ipp.latitude, ipp.longitude, le,15);//由穿刺点计算勒让德多项式的值
				 for (int w = 0;w < (O + 1)*(O + 1);w++)
				 {
					 if (flag == 1) { A1B(w, 0) += le[w] * (vtec_sh-a0); }//国内测站用观测值计算的vtec
					 else { A2B(w, 0) += le[w] * (vtec_sh-a0); }//国外测站用预报值计算的vtec
					 for (int q = 0;q < (O + 1)*(O + 1);q++)
					 {
						 if (flag == 1) { A1A(w, q) += le[w] * le[q]; }////国内测站
						 else { A2A(w, q) += le[w] * le[q]; }//国外测站
					 }
				 }
				 outfile << setw(4) << obsfile[0].obsdata[i].utime_o.year << setw(3) << obsfile[0].obsdata[i].utime_o.month << setw(3) << obsfile[0].obsdata[i].utime_o.day <<
					 setw(3) << obsfile[0].obsdata[i].utime_o.hour << setw(4) << obsfile[0].obsdata[i].utime_o.minute << setw(10) << obsfile[0].obsdata[i].utime_o.second;
				 outfile << std::right << setw(10) <<mjd<< std::right << setw(20) << ipp.latitude << std::right << setw(20) << ipp.longitude << std::right << setw(20) << vtec_sh << std::right << setw(20) << vtec_ion <<std::right << setw(20) << vtec_obs << std::right <<
					 setw(20) << vtec_predict <<  std::right << setw(5) <<flag << std::right << setw(10) << obsfile[j].obsdata[i].prn_list[k] << endl;
				 num++;
			 }
		 }
		 //解算过程
		 if (obsfile[0].obsdata[i].utime_o.hour%2 != 0 && obsfile[0].obsdata[i].utime_o.minute == 59 && obsfile[0].obsdata[i].utime_o.second == 30)//两个小时解算一次
		 {
			 cout << "时间：" << obsfile[0].obsdata[i].utime_o.hour << ":" << obsfile[0].obsdata[i].utime_o.minute << endl;
			 for (int w = 0;w < (O + 1)*(O + 1);w++)
			 {
				 BB(w, 0) = A1B(w, 0) / (sigma1*sigma1) + A2B(w, 0) / (sigma2*sigma2);
				 A1B(w, 0) = 0.0;//初始化为0
				 A2B(w, 0) = 0.0;//初始化为0
				 for (int q = 0;q < (O + 1)*(O + 1);q++)
				 {
					 AA(w, q) = A1A(w, q) / (sigma1*sigma1) + A2A(w, q) / (sigma2*sigma2);
					 A1A(w, q) = 0.0;//初始化为0
					 A2A(w, q) = 0.0;//初始化为0
				 }
			 }
			 if (num < (O + 1)*(O + 1)) { num = 0;continue; }//国内测站加上国外总的穿刺点个数要大于球谐函数系数个数
			 X = AA.inverse()*BB;//eigen库解方程
			 //chol_eq(AA, BB, X, (O + 1)*(O + 1));//cholesky分解法解方程
			 coe_update << setw(4) << obsfile[0].obsdata[i].utime_o.year << setw(3) << obsfile[0].obsdata[i].utime_o.month << setw(3) << obsfile[0].obsdata[i].utime_o.day <<
				 setw(3) << obsfile[0].obsdata[i].utime_o.hour << setw(4) << obsfile[0].obsdata[i].utime_o.minute << setw(10) << obsfile[0].obsdata[i].utime_o.second;
			 for (int w = 0;w < (O + 1)*(O + 1);w++)
			 {
				 coe_update <<std::right<< setw(15) << X(w, 0);
			 }
			 coe_update << endl;
		 }
	 }

	 outfile.close();
	 coe_update.close();
	 cout << "计算完毕" << endl;
 }
//找出一个文件夹里面的所有文件
 void getFiles(string path, pdcb dcbfile,vector<string>& files)
 {
	 //文件句柄  
	 long hFile = 0;//这个地方需要特别注意，win10用户必须用long long 类型，win7可以用long类型
	 //文件信息  
	 struct _finddata_t fileinfo;
	 string p,str;
	 double  dcb_val= 0.0;
	 if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	 {
		 do
		 {
			 //如果是目录,迭代之  
			 //如果不是,加入列表  
			 if ((fileinfo.attrib == _A_SUBDIR)&&strcmp(fileinfo.name, ".") && strcmp(fileinfo.name, ".."))
			 {
				getFiles(p.assign(path).append("\\").append(fileinfo.name),dcbfile, files);
			 }
			 else
			 {
				 if (!strcmp(fileinfo.name, ".") || !strcmp(fileinfo.name, ".."))continue;
				 str = fileinfo.name;
				 str = str.substr(0, 4);
				 transform(str.begin(), str.end(), str.begin(), ::toupper);//dcb文件里测站名字都是大写，所以转成大写
				 if (find_dcb(1, "G01", str, dcbfile, &dcb_val)){
					 files.push_back(p.assign(path).append("\\").append(fileinfo.name));
				 }
			 }
		 } while (_findnext(hFile, &fileinfo) == 0);
		 _findclose(hFile);
	 }
 }
