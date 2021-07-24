// BioALg.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "rapidcsv.h"
#include <bitset>
#include <time.h>
#include <random>
#include <cmath>

#include<stdio.h>	//to use the printf function
#include<conio.h>         		//to use the getche function
#include<stdlib.h>         		//to use the rand function

typedef struct Chrom             		// creating the chrom structure
{
	long long num;
	long long fit;
}chrom;                           	// now we have a chrom type that we can use

void* evpop(std::vector<chrom>& popcurrent, long long w, int* v, int poplen);    	//defining the functions that we will use
void* pickchroms(std::vector<chrom>& popcurrent, int poplen);
void* crossover(std::vector<chrom>& popnext, long long w, int* v, int poplen);
void* mutation(std::vector<chrom>& popnext, long long w, int* v, int poplen);
std::vector<double> main2(long long w, int* v, int pop);
void sortChr(std::vector<chrom>& popcurrent, int poplen);

#define _CRT_SECURE_NO_WARNINGS


//int v1[] = {791896991, 566700020, 947619652, 3137662, 615683440, 352054057, 755375718, 813723049, 863575073, 676707716, 930283504, 725982032, 866262633, 569987920, 553897669, 620357772, 614882721, 560519903, 750327497, 1019073458, 337732079, 349078249, 992866332, 897377591};
//int v2[] = {49238907, 30947850, 40087239, 90536176, 24468541, 24663140, 54080654, 69070243, 32347677, 98729654, 14329848, 28243771, 69977268, 53167239, 9965612, 73755224, 12484251, 100753481, 17923564, 94640161, 103786700, 29127456, 4322219, 79836774                       };
//int v3[] = {5101621, 2340206, 2425893, 13978068, 2361075, 9974574, 11314866, 5838561, 11737226, 8940366, 3955156, 3535685, 6662530, 232718, 16301123, 12917229, 7679120, 4646399, 15810499, 9232240, 1507540, 15849850, 10084082, 4372661                                        };
//int v4[] = {1633277, 2783551, 2336439, 1074217, 2975386, 2847124, 2434604, 3233948, 2454907, 852489, 3645634, 823774, 653279, 3254455, 33497, 580163, 157553, 2237722, 2164831, 655294, 2149935, 3157986, 1849632, 1136668                                                       };
//int v5[] = {23827, 497593, 414559, 488033, 837934, 561633, 726146, 224088, 317181, 1036517, 391406, 904974, 671804, 715959, 256468, 67852, 733075, 291434, 575568, 214191, 550982, 647805, 467845, 679641                                                                        };
//int v6[] = {336653, 71148, 254301, 299649, 123069, 326768, 196149, 187585, 88755, 39997, 236721, 356553, 122264, 189366, 179469, 218527, 140199, 207150, 12776, 109414, 112341, 181563, 10699, 70898                                                                             };
//int v7[] = {126890, 4791, 96740, 3164, 59203, 142590, 90544, 85210, 14730, 135818, 45985, 94631, 135634, 3753, 46862, 75932, 80190, 90419, 28007, 104454, 74840, 56601, 76648, 58696                                                                                             };

//int v1[] = { 401732219, 816883405, 853343681, 100807507, 346891520, 315272754, 19747021, 643088904, 1031613989, 1000705629, 987908291, 874960796, 704605826, 695944804, 581917807, 110401604, 1060911056, 301148943, 796734343, 885826844, 625260662, 63153414, 465318086, 805666887 };
//int v2[] = { 57607268, 16822471, 2359640, 20490871, 92245163, 47765374, 34425696, 48439307, 67686419, 5187649, 100530595, 55238513, 31322799, 106474236, 66730735, 66276192, 37649911, 104698908, 17972062, 73797972, 55599893, 88288356, 59389414, 52896714                         };
//int v3[] = { 919744, 1026472, 9964128, 584230, 12344692, 8048107, 774687, 8277535, 3057033, 5452985, 16191563, 9198708, 4259079, 9327451, 14814279, 15973738, 1205116, 1103019, 6681002, 7264291, 12067006, 3075268, 12572393, 12580513                                              };
//int v4[] = { 627104, 2226511, 600200, 1224888, 3372081, 1911316, 2340708, 268962, 3136282, 2725051, 2264816, 1880152, 2892695, 596137, 2336786, 2733730, 813426, 2443241, 864336, 3447328, 1085206, 464868, 2305439, 2703057                                                         };
//int v5[] = { 349297, 999519, 463862, 1014345, 310993, 671503, 436605, 967091, 208121, 395980, 984183, 310005, 554919, 122473, 897987, 323883, 499617, 1028141, 366000, 272610, 964346, 618665, 306131, 263721                                                                        };
//int v6[] = { 156484, 237623, 82623, 90372, 278437, 13678, 219100, 278075, 146616, 321835, 292327, 271866, 55619, 147506, 296275, 299557, 29002, 189761, 188482, 103531, 176117, 298591, 101460, 58837                                                                                };
//int v7[] = { 128834, 78544, 109178, 40529, 43395, 129211, 118050, 98762, 6770, 134868, 41381, 71235, 87240, 30, 14672, 23855, 38252, 15602, 75663, 113806, 34666, 14804, 334, 1084                                                                                                   };

//int v1[] = { 95509355, 28594435, 1049853931, 896821815, 901518516, 684163341, 366968010, 784851318, 64032268, 855617492, 641010614, 384921822, 120071133, 767691824, 453735229, 1064185178, 665070084, 50272341, 620443286, 603537908, 406327814, 809033607, 992568340, 76811112 };
//int v2[] = { 103989081, 102652687, 102428785, 28507510, 31386446, 13446237, 73387060, 98747600, 18582670, 86438057, 60484661, 88628258, 19287857, 103817362, 46879470, 19390770, 8937434, 20284482, 54226437, 61078821, 144888, 46982154, 5471356, 38027516                      };
//int v3[] = { 11303293, 10826650, 128325, 12921521, 15644389, 4477822, 1843930, 4397781, 8649736, 15442454, 2470994, 4117698, 5499943, 10841882, 4835556, 7816189, 4454991, 13364805, 14061760, 8176997, 211177, 16575893, 4593278, 597141                                        };
//int v4[] = { 2047690, 2811735, 2756427, 2471799, 1725381, 2502959, 1854135, 3636417, 3048983, 1086644, 2676755, 1213720, 3055047, 1677414, 3259789, 3080085, 920180, 3234894, 3035151, 942012, 1900680, 158001, 2537040, 1212186                                                 };
//int v5[] = { 566384, 73947, 553208, 563825, 880080, 178150, 939283, 9921, 825490, 152053, 126501, 983313, 1037774, 641811, 721923, 101964, 697236, 558247, 43539, 7281, 516277, 361022, 134646, 830852                                                                           };
//int v6[] = { 352990, 102498, 307868, 230361, 356721, 219163, 187735, 340440, 270949, 323060, 54848, 324383, 2943, 300888, 159669, 317766, 92683, 55197, 164315, 24195, 345303, 174084, 310864, 7533                                                                              };
//int v7[] = { 41319, 90539, 115593, 68520, 100694, 4347, 16295, 51146, 101386, 106277, 54590, 32721, 86833, 59011, 116162, 110993, 104298, 73488, 120407, 105768, 79078, 57978, 120396, 60963                                                                                     };


//int v1[] = { 698150968, 819733317, 776890390, 219271702, 108543181, 172035038, 1029215093, 540084952, 726478917, 1047590380, 996318190, 491311864, 99312855, 108935948, 57219787, 702071622, 54250709, 965966326, 762316212, 722081770, 717554858, 197448242, 812854106, 756528372 };
//int v2[] = { 32604997, 31794795, 40478451, 4407038, 23415780, 40512441, 5903378, 9878446, 93780818, 56286924, 59629162, 151272, 61400493, 64647174, 95167168, 40851301, 40845538, 14750511, 64803430, 93835142, 82956424, 45114109, 27126130, 106254416 };
//int v3[] = { 14879802, 15838382, 14961147, 12094281, 15303817, 11850038, 14970911, 11774765, 4710547, 9126577, 12556314, 5188875, 15545320, 9810215, 13772890, 7142626, 1823216, 1609031, 849338, 10425509, 2743599, 9037414, 11304004, 14397538 };
//int v4[] = { 2057270, 2687261, 3407287, 1969170, 1047871, 2136338, 3556588, 2470359, 2425073, 2049713, 2013692, 3330794, 2581951, 560870, 2693726, 925102, 916776, 183752, 2091185, 3230958, 2506914, 3215335, 295018, 3096572 };
//int v5[] = { 36310, 917268, 776569, 294817, 237732, 83542, 705564, 758709, 866529, 730993, 691243, 476693, 377562, 446696, 415565, 689597, 908453, 197167, 305223, 192018, 699384, 457632, 261464, 307579 };
//int v6[] = { 296178, 308054, 2876, 211377, 360646, 164462, 82857, 31810, 54693, 83978, 269760, 122825, 299939, 210037, 306400, 4306, 75819, 337814, 79777, 172221, 351012, 252214, 26952, 69470 };
//int v7[] = { 133437, 26440, 24307, 60059, 100557, 131280, 238, 117101, 62642, 15727, 33613, 19075, 89470, 23720, 90461, 103510, 70602, 29231, 104078, 78944, 45282, 86723, 32184, 53955 };

int v1[] = { 467354960, 209648272, 658444868, 842358638, 594807848, 295586341, 289228658, 731124395, 764873527, 65792377, 23255844, 666075975, 43477198, 380552658, 731988320, 540736146, 1009463575, 1050805005, 105038835, 651857940, 750321683, 848051054, 999713734, 636320593	};
int v2[] = { 46483381, 96610148, 99457889, 15256039, 51960367, 83221741, 99677022, 78195243, 64566584, 55517358, 82119995, 90788900, 15840823, 54672091, 75113516, 12784879, 8495643, 23283625, 12465434, 23574583, 10119012, 43142758, 22181778, 74304636							};
int v3[] = { 16210853, 9563381, 9827412, 4736256, 8358954, 3435938, 2901405, 3977036, 1660333, 12640581, 10309494, 6459738, 4747892, 2059728, 15471092, 477375, 8270274, 6797824, 11908780, 8098643, 9673238, 11396900, 12187629, 15710494											};
int v4[] = { 2185474, 1288018, 2533424, 2195122, 2730406, 2718761, 1535835, 2606784, 921529, 2432105, 144994, 56035, 473295, 574341, 147511, 2890454, 1473089, 2600098, 3396700, 1954218, 1568893, 2700548, 2597751, 932653															};
int v5[] = { 937103, 820564, 438349, 893944, 459378, 1007231, 581650, 456111, 17209, 73427, 19915, 888855, 993877, 994406, 978024, 870755, 848786, 861183, 953438, 715125, 141730, 521874, 303550, 825493																			};
int v6[] = { 211332, 94213, 8146, 258618, 218022, 250659, 216224, 111149, 311459, 315048, 171626, 18077, 67531, 130171, 15369, 75291, 82854, 326273, 153639, 109002, 78912, 124186, 149064, 153429																					};
int v7[] = { 18550, 142058, 15622, 60366, 2104, 93618, 114745, 65042, 108370, 19420, 12362, 12745, 71974, 25224, 143508, 133681, 126746, 10536, 137489, 120986, 82408, 16573, 140856, 99602																							};



int MAXg = pow(2, 24) - 1;
std::random_device rd; // obtain a random number from hardware
std::default_random_engine e1(rd());
//std::uniform_int_distribution<> distr0(0, MAXg); // define the range
std::uniform_int_distribution<int> uniform_dist0(0, MAXg);
int mean = uniform_dist0(e1);
std::seed_seq seed2{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
std::mt19937 gen(seed2);



std::vector<double> knps(unsigned long long w, int num, int* v)
{
	bool flag = false;
	std::vector<double> ans;

	int max = pow(2, 24) - 1;
	unsigned long long sum = 0;
	int found = 0;
	time_t start, end;
	start = clock();
	double time_taken = 0;

	for (int j = 0; j < max; j++)
	{
		unsigned long long sum = 0;
		std::bitset<24> X(j);
		for (int k = 0; k < 24; k++)
		{
			if (X[k] == 1)
			{
				sum += v[k];
			}
			if (sum > w)
			{
				break;
			}
		}
		if (w == sum)
		{
			if (!flag)
			{
				end = clock();
				time_taken = double(end - start) / double(CLOCKS_PER_SEC);
				flag = true;
			}
			found++;
		}
	}

	ans.push_back(time_taken);
	ans.push_back(found);

	return ans;
}


long long fitness(long long w, int Xn, int* v)
{

	long long sum = 0;
	std::bitset<24> X(Xn);

	for (int k = 0; k < 24; k++)
	{
		if (X[k] == 1)
		{
			sum += v[k];
		}

	}

	long long ans = abs(w - sum);

	return ans;
}

std::vector<int> genRand(int Num)
{
	const int AMOUNT = Num; //amount of random numbers that need to be generated
	int MAX = pow(2, 24) - 1;; //maximum value (of course, this must be at least the same as AMOUNT;

	int MIN = MAX;

	std::vector<int> value(Num);

	std::uniform_int_distribution<>* distr1 = nullptr;
	int d = MAX / Num;

	int k = d;
	//generate random numbers:
	for (int i = 0; i < Num; i++)
	{
		//if (!(i % d))
		//{


		MIN = abs(MIN - d);
		if (MIN < 0)
		{
			MIN = 0;
		}
		std::uniform_int_distribution<> distr(MIN, MAX); // define the range
		distr1 = &distr;
		


		//std::uniform_int_distribution<> distr(0, MAX); // define the range
		//distr1 = &distr;
		//MAX = abs(MAX - d);
		//}

		bool check; //variable to check or number is already used
		int n; //variable to store the number in

		//for (int n = 0; n < 40; ++n)
		//    std::cout <<  << ' '; // generate numbers
		do
		{
			if (distr1 != nullptr)
				n = (*distr1)(gen);
			else
				n = rand() % MAX;



			//check or number is already used:
			check = true;
			for (int j = 0; j < i; j++)
				if (n == value[j]) //if number is already used
				{
					check = false; //set check to false
					break; //no need to check the other elements of value[]
				}
		} while (!check); //loop until new, unique number is found
		value.at(i) = n; //store the generated number in the array
	}
	return value;
}

std::vector<int> genRand2(int Num)
{
	const int AMOUNT = Num; //amount of random numbers that need to be generated
	int MAX = pow(2, 24) - 1;; //maximum value (of course, this must be at least the same as AMOUNT;


	std::vector<int> value(Num);
	//int value[AMOUNT]; //array to store the random numbers in



	//std::normal_distribution<> normal_dist(mean, 16000000);


	//std::mt19937 gen(rd()); // seed the generator




	//std::uniform_int_distribution<> distr1(0, MAX/100000); // define the range
	//std::uniform_int_distribution<> distr2(0, MAX / 10000); // define the range
	//std::uniform_int_distribution<> distr3(0, MAX / 1000); // define the range
	//std::uniform_int_distribution<> distr4(0, MAX / 100); // define the range
	//std::uniform_int_distribution<> distr5(0, MAX / 10); // define the range

	std::uniform_int_distribution<>* distr1 = nullptr;
	std::uniform_int_distribution<> distr(0, pow(2, 24) - 1); // define the range
	distr1 = &distr;
	int d = MAX / Num;

	int k = d;
	//generate random numbers:
	for (int i = 0; i < Num; i++)
	{
		//if (!(i % d))
		//{

		//MAX = abs(MAX - d);
		//}

		bool check; //variable to check or number is already used
		int n; //variable to store the number in

		//for (int n = 0; n < 40; ++n)
		//    std::cout <<  << ' '; // generate numbers
		do
		{
			if (distr1 != nullptr)
				n = (*distr1)(gen);
			else
				n = rand() % MAX;

			//if (i >= 0 && i <= 20)
			//{
			//    n = distr1(gen);
			//}
			//if (i > 20 && i <= 40)
			//{
			//    n = distr2(gen);
			//}
			//if (i > 40 && i <= 60)
			//{
			//    n = distr3(gen);
			//}
			//if (i > 60 && i <= 80)
			//{
			//    n = distr4(gen);
			//}
			//if (i > 80 && i <= 100)
			//{
			//    n = distr5(gen);
			//}

			//check or number is already used:
			check = true;
			for (int j = 0; j < i; j++)
				if (n == value[j]) //if number is already used
				{
					check = false; //set check to false
					break; //no need to check the other elements of value[]
				}
		} while (!check); //loop until new, unique number is found
		value.at(i) = n; //store the generated number in the array
	}
	return value;
}

double coeffMut = 40;
double coeffCross = 40;
double coeffRed = 0.65;

int coeffMutv[] = { 1, 2, 5, 10, 15, 20, 25, 30, 35, 40 };
int coeffCrossv[] = { 1, 2, 5, 10, 15, 20, 25, 30, 35, 40 };
double coeffRedv[] = { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.65, 0.70 };
int popv[] = { 300, 350, 400, 450, 500, 550, 600 };

int main()
{


	//auto v1r = genRand(100);
	//auto v2r = genRand2(100);
	//
	//for (int k = 0; k < 100; k++)
	//{
	//    std::cout << v1r[k] << " ";
	//}
	//std::cout << std::endl;
	//for (int k = 0; k < 100; k++)
	//{
	//    std::cout << v2r[k] << " ";
	//}
	//std::cout << std::endl;

	rapidcsv::Document doc("E:\\prazia\\univer\\БИОАлг\\Лабы\\Problems550.csv");

	std::vector<double> w = doc.GetColumn<double>("Целевой вес");
	std::vector<int> num = doc.GetColumn<int>("Номер вектора");
	std::vector<std::vector<double>> ans;
	std::vector<std::vector<double>> ans2;

	//std::bitset<24> X1(514);

	//main2(w.at(0), v1, 300);
	//std::vector<std::pair<int, int>> ansCoef;
	//std::vector<std::pair<int, int>> ansPop;
	//std::vector<double> ansV;




	//for (int mut = 9; mut >= 0; mut--)
	//{
	//	coeffMut = coeffMutv[mut];
	//	for (int cross = 9; cross >= 0; cross--)
	//	{
	//		coeffCross = coeffCrossv[cross];
	//		for (int red = 10; red >= 0; red--)
	//		{
	//			coeffRed = coeffRedv[red];
	//			for (int pop = 0; pop < 7; pop++)
	//			{
	//				int amm = 0;
	//				double per = 0;
	//				for (int i = 0; i < 50; i++)
	//				{
	//					unsigned long long wN = w.at(i);

	//					auto v = main2(wN, v1, popv[pop]);
	//					if (v.at(0) == 0)
	//					{
	//						amm++;
	//					}
	//					std::cout << "i = " << i + 1 << " Ans = " << v.at(0) << std::endl;
	//				}
	//				ansV.push_back((amm * 100) / 50.0);
	//				std::pair<int, int> Coefp;
	//				std::pair<int, int> Popp;
	//				Coefp.first = mut;
	//				Coefp.second = cross;
	//				Popp.first = red;
	//				Popp.second = pop;
	//				ansCoef.push_back(Coefp);
	//				ansPop.push_back(Popp);
	//			}
	//			std::cout << std::endl << "NEW" << std::endl;
	//		}
	//	}
	//}

	//double maxAns = 0;
	//int ansNum = 0;
	//for (int i = 0; i < ansV.size(); i++)
	//{
	//	if (maxAns < ansV.at(i))
	//	{
	//		ansNum = i;
	//		maxAns = ansV.at(i);
	//	}
	//}

	//std::cout<<"ANSWER: " << std::endl << "Mut: "<< ansCoef.at(ansNum).first << std::endl << "Cross: " << ansCoef.at(ansNum).second << std::endl <<
	//	"Reduction: " << ansPop.at(ansNum).first << std::endl << "Population: " << ansPop.at(ansNum).second << std::endl;

	// GA edition \\

	std::cout << "______________________GA edition______________________" << std::endl;

	for (int i = 0; i < w.size(); i++)
	{
		time_t start, end;
		unsigned long long wN = w.at(i);
		int numN = num.at(i);
		std::vector<double> v;
		std::vector<double> vec2;


		start = clock();
		switch (numN)
		{
		case 1:
			v = main2(wN, v1, 300);
			break;
		case 2:
			v = main2(wN, v2, 300);
			break;
		case 3:
			v = main2(wN, v3, 300);
			break;
		case 4:
			v = main2(wN, v4, 300);
			break;
		case 5:
			v = main2(wN, v5, 300);
			break;
		case 6:
			v = main2(wN, v6, 300);
			break;
		case 7:
			v = main2(wN, v7, 300);
			break;
		default:
			break;
		}

		end = clock();

		std::cout << "i = " << i + 1 << " Time = " << double(end - start) / double(CLOCKS_PER_SEC) << " Ans = " << v.at(0) << std::endl;
		//v.push_back(double(end - start));

		vec2.push_back(double(i + 1));
		vec2.push_back(v.at(0));
		vec2.push_back(double(end - start) / double(CLOCKS_PER_SEC));
		vec2.push_back(v.at(1));
		vec2.push_back(v.at(2));

		ans.push_back(vec2);
	}

	std::ofstream myfile;
	myfile.open("SolveKS5502.csv");
	for (int i = 0; i < ans.size(); i++)
	{
		int num = ans.at(i).at(0);
		double tm1 = ans.at(i).at(1);
		double timeFl = ans.at(i).at(2);
		int ansW = ans.at(i).at(3);
		int ansW2 = ans.at(i).at(4);

		myfile << num << "," << tm1 << "," << timeFl << "," << ansW << "," << ansW2 << "\n";
	}
	myfile.close();

	// BF edition \\


	//std::cout << std::endl << "______________________BF edition______________________" << std::endl;

	//for (int i = 0; i < w.size(); i++)
	//{
	//    time_t start, end;
	//    unsigned long long wN = w.at(i);
	//    int numN = num.at(i);
	//    std::vector<double> v;
	//    std::vector<double> vec2;
	//
	//
	//    start = clock();
	//    switch (numN)
	//    {
	//    case 1:
	//        v = knps(wN, numN, v1);
	//        break;
	//    case 2:
	//        v = knps(wN, numN, v2);
	//        break;
	//    case 3:
	//        v = knps(wN, numN, v3);
	//        break;
	//    case 4:
	//        v = knps(wN, numN, v4);
	//        break;
	//    case 5:
	//        v = knps(wN, numN, v5);
	//        break;        
	//    case 6:
	//        v = knps(wN, numN, v6);
	//        break;        
	//    case 7:
	//        v = knps(wN, numN, v7);
	//        break;
	//    default:
	//        break;
	//    }
	//
	//    end = clock();
	//
	//    std::cout << "i = " << i+1 << " Time = " << double(end - start) / double(CLOCKS_PER_SEC) << std::endl;
	//    //v.push_back(double(end - start));
	//
	//    vec2.push_back(i + 1);
	//    vec2.push_back(v.at(0));
	//    vec2.push_back(double(end - start) / double(CLOCKS_PER_SEC));
	//    vec2.push_back(v.at(1));
	//
	//    ans2.push_back(vec2);
	//}

	//
	//std::ofstream myfile2;
	//myfile2.open("SolveBF550.csv");
	//for (int i = 0; i < ans2.size(); i++)
	//{
	//	int num = ans2.at(i).at(0);
	//	double tm1 = ans2.at(i).at(1);
	//	double timeFl = ans2.at(i).at(2);
	//	int ansW = ans2.at(i).at(3);

	//	myfile2 << num << "," << tm1 << "," << timeFl << "," << ansW << "," << "\n";
	//}
	//myfile2.close();

}






std::vector<double> main2(long long w, int* v, int poplen)								// the main function
{
	int num = 100;								// num is the no. of iterations
	int i, j;

	//printf("\nWelcome to the Genetic Algorithm coded by Loay Al Lusi:http://www.linkedin.com/in/loayallusi in May 2005 \nThe Algorithm is based on the function y = -x ^ 2 + 5 to find the maximum value of the function...\n"); // introduction to the program


	//enter: printf("\nPlease enter the no. of iterations:  ");
	//scanf("%d", &num);             	// enter the no. of iterations in num

	std::vector<chrom> popcurrent(poplen);                        	// we make 4 chromes of popcurrent 
	std::vector<chrom> popnext(poplen);                           	// we make 4 chromes of popnext

	//if (num < 1)                               	// if a -ve number is inserted .. enter num again
	//    goto enter;

	evpop(popcurrent, w, v, poplen);                       	//initialise pop current

	unsigned long long best = 0;
	bool flag = false;
	int am = 0;
	int flSt = -1;

	for (i = 0; i < 2000; i++)                      	// loop num times
	{

		//printf("\ni = %d\n", i);        	// print the iteration number

		for (j = 0; j < poplen; j++)
			popnext[j] = popcurrent[j];            	//copy popcurrent to popnext in order to adjust it

		pickchroms(popnext, popnext.size());                    	//pick best chromes
		crossover(popnext, w, v, popnext.size());                      	//cross over to get children chromes
		mutation(popnext, w, v, popnext.size());                       	//mutate with a low probability

		sortChr(popnext, popnext.size());

		if (popnext[0].fit == 0)
		{
			best = popnext[0].fit;
			flSt = 2;
			break;
		}
		if (!flag)
		{
			flag = !flag;
			best = popnext[0].fit;
		}
		else
		{
			if (best == popnext[0].fit)
			{
				am++;
				if (am >= 2)
				{
					flSt = 1;
					break;

				}

			}
			if (best < popnext[0].fit)
			{
				best = popnext[0].fit;
				am = 0;
			}
		}

		int div = poplen - popnext.size();

		auto vec = genRand(div);

		for (int k = 0; k < vec.size(); k++)
		{
			chrom popN;
			popN.fit = fitness(w, vec[k], v);//y(x(popcurrent[j]));
			popN.num = vec[k];
			popnext.push_back(popN);
		}


		for (j = 0; j < poplen; j++)
			popcurrent[j] = popnext[j];             	//copy the chromes of popnext to popcurrent 

	}

	// loop back until no. of iterations is exceeded
	pickchroms(popnext, popnext.size());
	//printf("\nPress any key to end ! ");

	std::vector<double> ans;
	ans.push_back(best);
	ans.push_back(flSt);
	ans.push_back(i);

	return ans;

	//flushall();                             		// flush the input buffer
	//getche();                                	// wait for a character from the keyboard to end

}                                            	//end of main

void reduction(std::vector<chrom>& popcurrent, int poplen)
{
	int i, j;
	chrom temp;                            	//temp chrome to use in sorting

	for (i = 0; i < poplen - 1; i++)               		//sorting the given set due to fitness
		for (j = 0; j < poplen - 1; j++)
		{
			if (popcurrent[j + 1].num < popcurrent[j].num)
			{
				temp = popcurrent[j + 1];
				popcurrent[j + 1] = popcurrent[j];
				popcurrent[j] = temp;

			}   // end of if
		}                // end of for loop

	/*for()*/
}

void* evpop(std::vector<chrom>& popcurrent, long long w, int* v, int poplen)               	//takes a pointer to a chrom of 4 elements
{
	int i, j, value;
	int random;
	//double amax = pow(2, 24) - 1;

	auto pop = genRand(poplen);

	for (int i = 0; i < pop.size(); i++)
	{
		popcurrent[i].fit = fitness(w, pop[i], v);//y(x(popcurrent[j]));
		popcurrent[i].num = pop[i];
	}




	//for (j = 0; j < 100; j++)                          // loop of j to choose chromes from [0] to [3]
	//{
	//    /*for (i = 0; i < 300; i++)            			// loop of i to choose the gen of the chrom from  [0] to [5]
	//    
	//    {
	//        random = rand();               		// creating random value
	//        random = (random % 2);        			// make the random value 0 or 1 only
	//        popcurrent[j].bit[i] = random;  		// initialising the bit[i] of chrom[j] with random
	//    } */  // end of for(i)
	//
	//    //random = rand();               		// creating random value
	//    //random = (random % 2);        			// make the random value 0 or 1 only
	//    //popcurrent[j].bit[i] = random;  		// initialising the bit[i] of chrom[j] with random
	//
	//    //value = x(popcurrent[j]);               	// get the value of the chrom as integer
	//    	// calcualte the fitness of chrom[j]
	//    //printf("\n popcurrent[%d]=%d%d%d%d%d%d    value=%d    fitness = %d", j,
	//    //    popcurrent[j].bit[5], popcurrent[j].bit[4], popcurrent[j].bit[3], popcurrent[j].bit[2],
	//    //    popcurrent[j].bit[1], popcurrent[j].bit[0], value, popcurrent[j].fit);     // print the chrom[j]
	//
	//}    // end of for(j)                                                              


	return(0);
}                              	//end of evpop function

//int x(chrom popcurrent)        	//x function that evaluate the value of a given chrom
//{
//    unsigned long long z;
//    z = (popcurrent.bit[0] * 1) + (popcurrent.bit[1] * 2) + (popcurrent.bit[2] * 4) + (popcurrent.bit[3] * 8) + (popcurrent.bit[4] * 16);
//    //if (popcurrent.bit[5] == 1)
//    //    z = z * (-1);                  	// z=sum of the ( bits * their weights) with the sign value   
//    return(z);                	//return the value of z
//}                             	// end x function
//
//int y(int x)          		// the y function that we look for it's maximum value takes x value
//{
//    int y;
//    y = -(x * x) + 5;            	// the function is y= - ( x^ 2 ) +5
//    return(y);
//}                             	// end of y function

void sortChr(std::vector<chrom>& popcurrent, int poplen)
{
	int i, j;
	chrom temp;                            	//temp chrome to use in sorting

	for (i = 0; i < poplen - 1; i++)               		//sorting the given set due to fitness
		for (j = 0; j < poplen - 1; j++)
		{
			if (popcurrent[j + 1].fit < popcurrent[j].fit)
			{
				temp = popcurrent[j + 1];
				popcurrent[j + 1] = popcurrent[j];
				popcurrent[j] = temp;

			}   // end of if
		}                // end
}

void* pickchroms(std::vector<chrom>& popcurrent, int poplen)   	// pickchroms takes a pointer to array of chroms
{

	int i, j;
	chrom temp;                            	//temp chrome to use in sorting

	for (i = 0; i < poplen - 1; i++)               		//sorting the given set due to fitness
		for (j = 0; j < poplen - 1; j++)
		{
			if (popcurrent[j + 1].fit < popcurrent[j].fit)
			{
				temp = popcurrent[j + 1];
				popcurrent[j + 1] = popcurrent[j];
				popcurrent[j] = temp;

			}   // end of if
		}                // end of for loop

	//for (i = 0; i < 4; i++)
	//    printf("\nSorting:popnext[%d] fitness=%d", i, popcurrent[i].fit);   	//printing the result

	int perc = coeffRed * poplen;


	std::vector<chrom> popout = std::vector<chrom>(popcurrent.begin(), popcurrent.end() - perc);
	popcurrent = popout;
	//printf("\n");                 //print new line
	//flushall();                                                       //flush the input buffer
	return(0);
}                       //end of pick chromes function

void* crossover(std::vector<chrom>& popnext, long long w, int* v, int poplen) // crossover function takes a pointer to array of chromes
{
	int random1;
	int i;
	random1 = rand();                                  	//random cross over point
	random1 = ((random1 % (int)coeffCross));                    		// cross point should be between (1 - 5)

	std::uniform_int_distribution<int> uniform_dist(0, 24 - 1);
	std::uniform_int_distribution<int> uniform_dist2(0, poplen - 1);


	int crossover_point = uniform_dist(gen); //rand() % 24;

	for (int i = 1; i < random1; i++)
	{
		int random21 = uniform_dist2(gen);                                  	//random cross over point
		int random22 = uniform_dist2(gen);

		int val1 = popnext[random21].num;
		int val2 = popnext[random22].num;

		std::bitset<24> X1(val1);
		std::bitset<24> X2(val2);

		for (int k = 0; k < crossover_point; k++)
		{
			X1[k] = X2[k];
		}
		for (int k = crossover_point; k < 24; k++)
		{
			X2[k] = X1[k];
		}
		//for (int k = 15; k < 24; k++)
		//{
		//    X2[k] = X1[k];
		//}

		popnext[random21].num = X1.to_ulong();
		popnext[random22].num = X2.to_ulong();
	}

	//for (i = 0; i < random; i++)                     	//crossing the bits below the cross point index
	//{
	//    popnext[2].num = popnext[0].num;        	//child 1 cross over
	//    popnext[3].num = popnext[1].num;     	// child 2 cross over
	//} // end of for
	//
	//for (i = random; i < 6; i++)                        	// crossing the bits beyond the cross point index
	//{
	//    popnext[2].num = popnext[1].num;     	// child 1 cross over
	//    popnext[3].num = popnext[0].num;       	// chlid 2 cross over
	//}    // end of for
	//
	for (i = 0; i < poplen; i++)
		popnext[i].fit = fitness(w, popnext[i].num, v);     	// calculating the fitness values for the new set
	//.num
	//for (i = 0; i < 4; i++)
	//    printf("\nCrossOver popnext[%d]=%d    value=%d    fitness = %d", i,
	//        popnext[i].num, fitness(w, popnext[i].num, v), popnext[i].fit);
	// printing the bits, value and fitness for the chromes of the new set

	return(0);
}                  // end crossover function

void* mutation(std::vector<chrom>& popnext, long long w, int* v, int poplen)   // mutation funtion given a pointer to array of chromes
{

	int random;
	int row, col, value;
	//random value is between ( 0 - 49 )
	std::uniform_int_distribution<int> uniform_dist(0, 23);
	std::uniform_int_distribution<int> uniform_dist2(0, poplen - 1);
	std::uniform_int_distribution<int> uniform_dist3(0, 1);

	int amm = 0;

	for (int i = 0; i < poplen; i++)
	{
		random = uniform_dist2(gen);
		if (random > coeffMut)    // Suppusiong Probability of mutation is 50 % 
		{
			col = uniform_dist(gen);                           	// random column (gene) choosing 
			row = uniform_dist2(gen);                           	// random row ( chrome ) choosing

			std::bitset<24> X1(popnext[row].num);


			int mutrand = uniform_dist3(gen);
			if (X1[col] != mutrand)
				amm++;
			X1[col] = mutrand;

			//if (X1[col] == 0)          	// invert the bit to 1 if it was 0              // 
			//    X1[col] = 1;                                                              // umutate bits in num
			//else X1[col] = 0;


			popnext[row].num = X1.to_ulong();
			popnext[row].fit = fitness(w, popnext[row].num, v);   	// calculate the fitness for the mutated chrome
			//value = popnext[row].fit;
			//printf("\nMutation occured in popnext[%d] bit[%d]:=%d   fitness = % d", row, col, popnext[row].num, popnext[row].fit);

			// print the chrome index,bits,value, fintness of the mutated chrome 
		}
	}

	//std::cout << "Population: " << poplen << " Mutated: " << (amm * 100.0) / poplen << " %\n";

	return(0);
}                       //end of mutation

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
