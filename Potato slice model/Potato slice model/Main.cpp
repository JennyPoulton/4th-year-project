#include"Polymer.h"

#include<iostream>
#include<fstream>
#include<ctime>
#include<iomanip>
#include <stdlib.h>

using namespace std;

int main(void)
{

	srand(time(NULL));
	//srand(5);

	ofstream Output1("OutputForce.txt");

	if (!Output1.is_open())
	{
		cout << "Error: output file 1 cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output 1 file opened successfully." << endl;
	}

	ofstream Output2("OutputSpringConstant.txt");

	if (!Output2.is_open())
	{
		cout << "Error: output file 1 cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output 2 file opened successfully." << endl;
	}

	ofstream Output3("OutputPictures.txt");

	if (!Output3.is_open())
	{
		cout << "Error: output file 1 cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output 3 file opened successfully." << endl;
	}

	Polymer StaphAureous;

	double Force = 5;

	for (int n = 0; n < DIMENSION*DIMENSION; n++)
	{
		StaphAureous.Set_Forces_And_Lengths(Force);
		StaphAureous.Calculate_Spring_Constant_Horizontal();
		StaphAureous.Calculate_Spring_Constant_Vertical();

		Output2 << StaphAureous.Return_Spring_Constant_Horizontal() << "\t" << StaphAureous.Return_Spring_Constant_Verticle() << endl;

		StaphAureous.Sort_Lengths_Into_Groups_For_Histogram();

		double Histogram_Glycan[100];
		double Histogram_Peptide[100];

		for (int i = 0; i < 100; i++)
		{
			Histogram_Glycan[i] = 0;
			Histogram_Peptide[i] = 0;
		}

		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				//Histogram_Glycan[(int)((StaphAureous.Return_Length_Glycan(i, j)-StaphAureous.Return_Min_Length_G()) / StaphAureous.Return_Splitter_For_Sorting_G())]++;
				double p = (StaphAureous.Return_Length_Peptide(i, j) - StaphAureous.Return_Min_Length_P());
				double q = StaphAureous.Return_Splitter_For_Sorting_P();

				double t = (StaphAureous.Return_Length_Glycan(i, j) - StaphAureous.Return_Min_Length_G());
				double s = StaphAureous.Return_Splitter_For_Sorting_G();

				for (double r = 0; r < 100; r++)
				{
					if (p >= r*q&&p <= (r + 1)*q)
					{
						Histogram_Peptide[(int)r]++;
					}
				}

				for (double r = 0; r < 100; r++)
				{
					if (t >= r*s&&t <= (r + 1)*s)
					{
						Histogram_Glycan[(int)r]++;
					}
				}

			}
		}


		for (int i = 0; i < 10; i++)
		{
			Output1 << i*StaphAureous.Return_Splitter_For_Sorting_G() << "\t" << Histogram_Glycan[i] << "\t" << i*StaphAureous.Return_Splitter_For_Sorting_P() << "\t" << Histogram_Peptide[i] << endl;
		}

		Output3 << "Peptide" << endl;

		for (int j = 0; j < DIMENSION; j++)
		{
			Output3 << j << "\t";

			for (int k = 0; k < DIMENSION; k++)
			{
				Output3 << StaphAureous.Return_Number_Bonds_Peptide(j, k) << "\t";
			}

			Output3 << endl;
		}

		Output3 << "Glycan" << endl;

		for (int j = 0; j < DIMENSION; j++)
		{
			Output3 << j << "\t";

			for (int k = 0; k < DIMENSION; k++)
			{
				Output3 << StaphAureous.Return_Number_Bonds_Glycan(j, k) << "\t";
			}

			Output3 << endl;
		}

		StaphAureous.Break_Bond();
	}


	//StaphAureous.Set_Forces_And_Lengths(Force);

	//StaphAureous.Sort_Lengths_Into_Groups_For_Histogram();



	system("pause");
	return 0;


}