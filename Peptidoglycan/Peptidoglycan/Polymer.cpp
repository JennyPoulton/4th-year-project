#include"Polymer.h"
#include<iostream>
using namespace std;


Polymer::Polymer()
{
	Max_Length_Glycan = 0;
	Max_Length_Peptide = 0;
	Min_Length_Glycan = 1000000000000000;
	Min_Length_Peptide = 1000000000000000;

	Set_Forces_And_Lengths(10);
	//Delete_Dead_Ends_And_Unjoined();

	Calculate_Spring_Constant_Horizontal();
	Calculate_Spring_Constant_Vertical();

}//in which the forces through each spring set to zero and the spring constants in each direction are found

void Polymer::Delete_Dead_Ends_And_Unjoined()
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			if (Force_Downwards[i][j] == 0 || Force_Upwards[i][j] == 0)
			{
				Murein[i][j].Set_Number_Bonds_Peptide(0);
			}
		}
	}

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			if (Force_Leftwards[i][j] == 0 || Force_Rightwards[i][j] == 0)
			{
				Murein[i][j].Set_Number_Bonds_Glycan(0);
			}
		}
	}

}

void Polymer::Calculate_Spring_Constant_Horizontal()
{
	double Number[DIMENSION];
	for (int i = 0; i < DIMENSION; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			Number[i] = Number[i] + PEPTIDE_SPRING_CONSTANT*Murein[j][i].Return_Number_Bonds_Peptide();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < DIMENSION; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Horizontal = Spring_Constant;
	return;
}

void Polymer::Calculate_Spring_Constant_Vertical()
{
	double Number[DIMENSION];
	for (int i = 0; i < DIMENSION; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			Number[i] = Number[i] + GLYCAN_SPRING_CONSTANT*Murein[i][j].Return_Number_Bonds_Glycan();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < DIMENSION; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Verticle = Spring_Constant;

	return;
}

void Polymer::Set_Forces_And_Lengths(double input_force)
{
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{

			Find_Force_Upwards(input_force, i, j);
			Find_Force_Leftwards(input_force, i, j);

			Murein[i][j].Set_Vertical_Force(Force_Upwards[i][j]);
			Murein[i][j].Set_Length_Peptide(Force_Upwards[i][j] / (double)PEPTIDE_SPRING_CONSTANT);

			Murein[i][j].Set_Horizontal_Force(Force_Leftwards[i][j]);
			Murein[i][j].Set_Length_Glycan(Force_Leftwards[i][j] / (double)PEPTIDE_SPRING_CONSTANT);
		}
	}
}

void Polymer::Find_Force_Upwards(double Input_Force, int p, int q)
{

	int length_of_current_bar_left = 0;
	int length_of_current_bar_right = 0;
	int length_of_above_bar_left = 0;
	int length_of_above_bar_right = 0;
	int length_of_below_bar_left = 0;
	int length_of_below_bar_right = 0;

	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
	{
		Force_Upwards[p][q] = 0;
		return;
	}

	for (int scroll_up = 0; scroll_up <= p; scroll_up++)
	{
		double Numerator = 1;
		double Denominator = 1;
		int join_up_tally = 0;
		int join_down_tally = 0;

		//glycans and q run horizontal, peptides and p run vertical, this set here
		

		while (((q - length_of_current_bar_left) > 0) && ( (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1) || (length_of_current_bar_left <= length_of_below_bar_left)) )
		{
			if (Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() == 1)
			{
				if (Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Glycan() == 1 || Murein[p - scroll_up][q - length_of_current_bar_left + 1].Return_Number_Bonds_Glycan() == 1)
				{
					join_up_tally++;
				}
			}
			
			length_of_current_bar_left++;

		}

		while (((q + length_of_current_bar_right) < DIMENSION-1) && (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1 || length_of_current_bar_right <= length_of_below_bar_right))
		{
			if (Murein[p - scroll_up][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() == 1)
			{
				join_up_tally++;
			}

			length_of_current_bar_right++;
		}

		
		if (p - scroll_up != 0)
		{
			while (((q - length_of_above_bar_left) > 0) && (Murein[p + 1 - scroll_up][q].Return_Number_Bonds_Glycan() == 1 || length_of_above_bar_left <= length_of_current_bar_left))
			{

				if (Murein[p - scroll_up][q - length_of_above_bar_left].Return_Number_Bonds_Peptide() == 1)
				{
					join_down_tally++;
				}

				length_of_above_bar_left++;

			}

			while (((q + length_of_above_bar_right) < DIMENSION-1) && (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1 || length_of_above_bar_right <= length_of_current_bar_right))
			{
				if (Murein[p - scroll_up][q + length_of_above_bar_right].Return_Number_Bonds_Peptide() == 1)
				{
					join_down_tally++;
				}

				length_of_above_bar_right++;
			}
		}
		else
		{
			join_down_tally = DIMENSION; //only true for glycans, for peptides it loops round
		}
		
		join_up_tally = join_up_tally - Murein[p - scroll_up][q].Return_Number_Bonds_Peptide();
		join_down_tally = join_down_tally - Murein[p - scroll_up][q].Return_Number_Bonds_Peptide();

		Numerator = Numerator*join_up_tally;
		Denominator = Denominator*join_down_tally;

		length_of_below_bar_left = length_of_current_bar_left;
		length_of_below_bar_right = length_of_current_bar_right;

		length_of_current_bar_left = 0;
		length_of_current_bar_right = 0;

		length_of_above_bar_left = 0;
		length_of_above_bar_right = 0;

		

		Input_Force = Input_Force*Numerator / Denominator;
	}

	Force_Upwards[p][q] = Input_Force;
	return;

}

void Polymer::Find_Force_Downwards(double Input_Force, int p, int q)
{

	int length_of_current_bar_left = 0;
	int length_of_current_bar_right = 0;
	int length_of_below_bar_left = 0;
	int length_of_below_bar_right = 0;
	int length_of_above_bar_left = 0;
	int length_of_above_bar_right = 0;

	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
	{
		Force_Downwards[p][q] = 0;
		return;
	}


	for (int scroll_down = 0; scroll_down <= DIMENSION-1 - p; scroll_down++)
	{
		double Numerator = 1;
		double Denominator = 1;
		int join_up_tally = 0;
		int join_down_tally = 0;

		//glycans and q run horizontal, peptides and p run vertical, this set here

		while (q - length_of_current_bar_left > 0 && (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 || length_of_current_bar_left <= length_of_above_bar_left))
		{

			if (Murein[p + scroll_down][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() == 1)
			{
				join_up_tally++;
			}

			length_of_current_bar_left++;

		}

		while (q + length_of_current_bar_right < DIMENSION-1 && (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 || length_of_current_bar_right <= length_of_above_bar_right))
		{
			if (Murein[p + scroll_down][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() == 1)
			{
				join_up_tally++;
			}

			length_of_current_bar_right++;
		}

		if (p + scroll_down != DIMENSION-1)
		{
			while (q - length_of_below_bar_left > 0 && (Murein[p + 1 + scroll_down][q].Return_Number_Bonds_Glycan() == 1 || length_of_below_bar_left <= length_of_current_bar_left))
			{

				if (Murein[p + scroll_down][q - length_of_below_bar_left].Return_Number_Bonds_Peptide() == 1)
				{
					join_down_tally++;
				}

				length_of_below_bar_left++;

			}

			while (q + length_of_below_bar_right < DIMENSION-1 && (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 || length_of_below_bar_right <= length_of_current_bar_right))
			{
				if (Murein[p + scroll_down][q + length_of_below_bar_right].Return_Number_Bonds_Peptide() == 1)
				{
					join_down_tally++;
				}

				length_of_below_bar_right++;
			}
		}
		else
		{
			join_down_tally = DIMENSION; //only true for glycans, for peptides it loops round
		}

		join_up_tally = join_up_tally - Murein[p + scroll_down][q].Return_Number_Bonds_Peptide();
		join_down_tally = join_down_tally - Murein[p + scroll_down][q].Return_Number_Bonds_Peptide();

		Numerator = Numerator*join_up_tally;
		Denominator = Denominator*join_down_tally;

		length_of_above_bar_left = length_of_current_bar_left;
		length_of_above_bar_right = length_of_current_bar_right;

		length_of_current_bar_left = 0;
		length_of_current_bar_right = 0;

		length_of_below_bar_left = 0;
		length_of_below_bar_right = 0;

		Input_Force = Input_Force*Numerator / Denominator;
	}

	Force_Downwards[p][q] = Input_Force;
	return;

}



void Polymer::Find_Force_Leftwards(double Input_Force, int q, int p)
{

	if (Murein[q][p].Return_Number_Bonds_Glycan() == 0)
	{
		Force_Leftwards[q][p] = 0;
		return;
	}
	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
	//first we need to find what length levels p-1 and p-2 are

	int leftward_extent_above[DIMENSION]; //n can take values less that p and represents the level
	int rightward_extent_above[DIMENSION]; //n can take values less that p and represents the level

	int total_left = 0; // this total represents the current length of all the bars being taken into account

	for (int n = p - 1; n > 0; n--)
	{
		int m = q;

		leftward_extent_above[n] = 0;

		while ((Murein[n][m].Return_Number_Bonds_Peptide() == 1 || leftward_extent_above[n] < total_left) && m <= DIMENSION - 1 && m >= 0)
		{
			leftward_extent_above[n]++;
			m--;
		}

		total_left = leftward_extent_above[n];



	}

	int total_right = 0; // this total represents the current length of all the bars being taken into account

	for (int n = p - 1; n > 0; n--)
	{

		rightward_extent_above[n] = 0;
		int m = q;

		while ((Murein[n][m].Return_Number_Bonds_Peptide() == 1 || rightward_extent_above[n] < total_right) && m <= DIMENSION - 1 && m >= 0)
		{
			rightward_extent_above[n]++;
			m++;
		}

		total_right = rightward_extent_above[n];


	}

	if (p == 0)
	{
		leftward_extent_above[p] = q;
		rightward_extent_above[p] = DIMENSION - 1 - q;
	}


	//this finds the lengths of all the above levels

	//now the numerator is the number of joins which connects anything within these levels to the one above it

	Numerator[p][q] = 1;

	for (int n = p - 1; n >= 0; n--)
	{
		int m = q;
		int tally_Glycan = 0;

		for (int i = 0; i < leftward_extent_above[n]; i++)
		{
			if (Murein[n][m - i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}

		}

		for (int i = 0; i < rightward_extent_above[n]; i++)
		{
			if (Murein[n][m + i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}
		}

		tally_Glycan = tally_Glycan - Murein[n][q].Return_Number_Bonds_Glycan();

		if (n - 1 >= DIMENSION - 1 || n - 1 <= 0)
		{
			tally_Glycan = DIMENSION;
		}


		Numerator[p][q] = Numerator[p][q] * (double)tally_Glycan;

	}

	Denominator[p][q] = 1;

	for (int n = p - 1; n >= 0; n--)
	{
		int m = q;
		int tally_Glycan = 0;

		for (int i = 0; i < leftward_extent_above[n]; i++)
		{
			if (Murein[n - 1][m - i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}

		}

		for (int i = 0; i < rightward_extent_above[n]; i++)
		{
			if (Murein[n - 1][m + i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}


		}

		tally_Glycan = tally_Glycan - Murein[n][q].Return_Number_Bonds_Glycan();

		if (n + 1 >= DIMENSION - 1 || n <= 0)
		{
			tally_Glycan = DIMENSION;
		}

		Denominator[p][q] = Denominator[p][q] * (double)tally_Glycan;

	}

	int m = q;

	int leftward_extent = 0;

	while (Murein[p][m].Return_Number_Bonds_Peptide() == 1 && m <= DIMENSION - 1 && m >= 0)
	{
		leftward_extent++;
		m--;
	}

	int rightward_extent = 0;
	m = q;

	while (Murein[p][m].Return_Number_Bonds_Peptide() == 1 && m <= DIMENSION - 1 && m >= 0)
	{
		rightward_extent++;
		m++;
	}

	m = q;
	int tally_Glycan = 0;

	for (int i = 0; i < leftward_extent; i++)
	{
		if (Murein[p][m - i].Return_Number_Bonds_Glycan() == 1)
		{
			tally_Glycan++;
		}

	}

	for (int i = 0; i < rightward_extent; i++)
	{
		if (Murein[p][m + i].Return_Number_Bonds_Glycan() == 1)
		{
			tally_Glycan++;
		}
	}

	tally_Glycan = tally_Glycan - Murein[p][q].Return_Number_Bonds_Glycan();

	double Frac = (double)Murein[p][q].Return_Number_Bonds_Glycan() / (double)tally_Glycan;

	Force_Leftwards[q][p] = Frac*Input_Force*Numerator[p][q] / Denominator[p][q];

}

void Polymer::Find_Force_Rightwards(double Input_Force, int q, int p)
{
	if (Murein[q][p].Return_Number_Bonds_Glycan() == 0)
	{
		Force_Leftwards[q][p] = 0;
		return;
	}

	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
	//first we need to find what length levels p-1 and p-2 are

	int leftward_extent_below[DIMENSION]; //n can take values less that p and represents the level
	int rightward_extent_below[DIMENSION]; //n can take values less that p and represents the level

	int total_left = 0; // this total represents the current length of all the bars being taken into account

	for (int n = p + 1; n < DIMENSION; n++)
	{
		int m = q;

		leftward_extent_below[n] = 0;

		while ((Murein[n][m].Return_Number_Bonds_Peptide() == 1 || leftward_extent_below[n] < total_left) && m <= DIMENSION - 1 && m >= 0)
		{
			leftward_extent_below[n]++;
			m--;
		}

		total_left = leftward_extent_below[n];



	}

	int total_right = 0; // this total represents the current length of all the bars being taken into account

	for (int n = p + 1; n < DIMENSION; n++)
	{

		rightward_extent_below[n] = 0;
		int m = q;
		int end = 1;

		while ((Murein[n][m].Return_Number_Bonds_Peptide() || rightward_extent_below[n] < total_right) && m <= DIMENSION - 1 && m >= 0)
		{
			rightward_extent_below[n]++;
			m++;
		}

		total_right = rightward_extent_below[n];


	}

	if (p == DIMENSION - 1)
	{
		leftward_extent_below[p] = q;
		rightward_extent_below[p] = DIMENSION - 1 - q;
	}


	//this finds the lengths of all the above levels

	//now the numerator is the number of joins which connects anything within these levels to the one above it

	Numerator[p][q] = 1;

	for (int n = p + 1; n < DIMENSION; n++)
	{
		int m = q;
		int tally_Glycan = 0;

		for (int i = 0; i < leftward_extent_below[n]; i++)
		{
			if (Murein[n - 1][m - i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}

		}

		for (int i = 0; i < rightward_extent_below[n]; i++)
		{
			if (Return_Number_Bonds_Glycan(n, m + i) == 1)
			{
				tally_Glycan++;
			}

		}

		tally_Glycan = tally_Glycan - Murein[n][q].Return_Number_Bonds_Glycan();

		if (n - 1 >= DIMENSION - 1 || n - 1 <= 0)
		{
			tally_Glycan = DIMENSION;
		}

		Numerator[p][q] = Numerator[p][q] * (double)tally_Glycan;

	}

	Denominator[p][q] = 1;

	for (int n = p + 1; n < DIMENSION; n++)
	{
		int m = q;
		int tally_Glycan = 0;

		for (int i = 0; i < leftward_extent_below[n]; i++)
		{
			if (Murein[n - 1][m - i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}

		}

		for (int i = 0; i < rightward_extent_below[n]; i++)
		{
			if (Murein[n - 1][m + i].Return_Number_Bonds_Glycan() == 1)
			{
				tally_Glycan++;
			}

		}

		tally_Glycan = tally_Glycan - Murein[n][q].Return_Number_Bonds_Glycan();

		if (n + 1 >= DIMENSION - 1 || n <= 0)
		{
			tally_Glycan = DIMENSION;
		}

		Denominator[p][q] = Denominator[p][q] * (double)tally_Glycan;

	}

	Force_Rightwards[q][p] = Input_Force*Numerator[p][q] / Denominator[p][q];

}

void Polymer::Break_Bond()
{
	//peptide=0, glycan=1
	if (Max_Force_Peptide_Or_Glycan == 0)
	{
		Murein[Max_Force_Vertical_Coordinate][Max_Force_Horizontal_Coordinate].Set_Number_Bonds_Peptide(0);
		return;
	}
	else if (Max_Force_Peptide_Or_Glycan == 1)
	{
		Murein[Max_Force_Vertical_Coordinate][Max_Force_Horizontal_Coordinate].Set_Number_Bonds_Glycan(0);
		return;
	}

	cout << "Break Bond Broken" << endl;
}

double Polymer::Return_Spring_Constant_Horizontal()
{
	return Spring_Constant_Horizontal;
}

double Polymer::Return_Spring_Constant_Verticle()
{
	return Spring_Constant_Verticle;
}

void Polymer::Sort_Lengths_Into_Groups_For_Histogram()
{


	int Max_Peptide_Vert;
	int Max_Peptide_Horiz;
	int Max_Glycan_Vert;
	int Max_Glycan_Horiz;

	Max_Length_Peptide = 0;
	Max_Length_Glycan = 0;
	Min_Length_Glycan = 1000000000000000;
	Min_Length_Peptide = 10000000000000000;

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			cout << i << "\t" << j << "\t" << Murein[i][j].Return_Length_Peptide() << "\t" << Murein[i][j].Return_Vertical_Force() << endl;

			if (Murein[i][j].Return_Length_Peptide() > Max_Length_Peptide)
			{
				Max_Length_Peptide = Murein[i][j].Return_Length_Peptide();
				Max_Peptide_Vert = i;
				Max_Peptide_Horiz = j;

			}

			if (Murein[i][j].Return_Length_Peptide() < Min_Length_Peptide)
			{
				Min_Length_Peptide = Murein[i][j].Return_Length_Peptide();
			}


			if (Murein[i][j].Return_Length_Glycan() > Max_Length_Glycan)
			{
				Max_Length_Glycan = Murein[i][j].Return_Length_Glycan();
				Max_Glycan_Vert = i;
				Max_Glycan_Horiz = j;
			}

			if (Murein[i][j].Return_Length_Glycan() < Min_Length_Glycan)
			{
				Min_Length_Glycan = Murein[i][j].Return_Length_Glycan();
			}

		}
	}

	if (Murein[Max_Glycan_Vert][Max_Glycan_Horiz].Return_Horizontal_Force() > Murein[Max_Peptide_Vert][Max_Peptide_Horiz].Return_Vertical_Force())
	{
		Max_Force_Horizontal_Coordinate = Max_Glycan_Horiz;
		Max_Force_Vertical_Coordinate = Max_Glycan_Vert;
		Max_Force_Peptide_Or_Glycan = 1;
	}
	else
	{
		Max_Force_Horizontal_Coordinate = Max_Peptide_Horiz;
		Max_Force_Vertical_Coordinate = Max_Peptide_Vert;
		Max_Force_Peptide_Or_Glycan = 0;
	}


	Splitter_For_Sorting_G = (Max_Length_Glycan - Min_Length_Glycan) / (double)100;
	Splitter_For_Sorting_P = (Max_Length_Peptide - Min_Length_Peptide) / (double)100;

	if (Min_Length_Peptide == Max_Length_Peptide)
	{
		Splitter_For_Sorting_P = 1;
	}

	if (Min_Length_Glycan == Max_Length_Glycan)
	{
		Splitter_For_Sorting_G = 1;
	}

}// ensure if min=max then lengths remain the same

double Polymer::Return_Splitter_For_Sorting_G()
{
	return Splitter_For_Sorting_G;
}

double Polymer::Return_Splitter_For_Sorting_P()
{
	return Splitter_For_Sorting_P;
}

void Polymer::Set_Length_Glycan(double l, int p, int q)
{
	Murein[p][q].Set_Length_Glycan(l);
	return;
}

double Polymer::Return_Length_Glycan(int p, int q)
{
	return Murein[p][q].Return_Length_Glycan();
}

void Polymer::Set_Spring_Constant_Glycan(double k, int p, int q)
{
	Murein[p][q].Set_Spring_Constant_Glycan(k);
}

double Polymer::Return_Spring_Constant_Glycan(int p, int q)
{
	return Murein[p][q].Return_Spring_Constant_Glycan();
}

void Polymer::Set_Number_Bonds_Glycan(int n, int p, int q)
{
	Murein[p][q].Set_Number_Bonds_Glycan(n);
}

int Polymer::Return_Number_Bonds_Glycan(int p, int q)
{
	return Murein[p][q].Return_Number_Bonds_Glycan();
}

void Polymer::Set_Length_Peptide(double l, int p, int q)
{
	Murein[p][q].Set_Length_Peptide(l);
	return;
}

double Polymer::Return_Length_Peptide(int p, int q)
{
	return Murein[p][q].Return_Length_Peptide();
}

void Polymer::Set_Spring_Constant_Peptide(double k, int p, int q)
{
	Murein[p][q].Set_Spring_Constant_Peptide(k);
	return;
}

double Polymer::Return_Spring_Constant_Peptide(int p, int q)
{
	return Murein[p][q].Return_Spring_Constant_Peptide();
}

void Polymer::Set_Number_Bonds_Peptide(int n, int p, int q)
{
	Murein[p][q].Set_Number_Bonds_Peptide(n);
	return;
}

int Polymer::Return_Number_Bonds_Peptide(int p, int q)
{
	return Murein[p][q].Return_Number_Bonds_Peptide();
}

double Polymer::Return_Min_Length_G()
{
	return Min_Length_Glycan;
}

double Polymer::Return_Min_Length_P()
{
	return Min_Length_Peptide;
}