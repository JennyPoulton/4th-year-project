#include"Polymer.h"
#include<iostream>
#include<cmath>
using namespace std;

Polymer::Polymer()
{
	Max_Length_Glycan = 0;
	Max_Length_Peptide = 0;
	Min_Length_Glycan = 1000000000000000;
	Min_Length_Peptide = 1000000000000000;

	Set_Forces_And_Lengths(10);
	Delete_Dead_Ends_And_Unjoined();

	Calculate_Spring_Constant_Horizontal();
	Calculate_Spring_Constant_Vertical();

}//in which the forces through each spring set to zero and the spring constants in each direction are found

void Polymer::Set_Up_Potato()
{
	Circumferance_at_centre= (int)(3.14159*(double)DIMENSION);
	for (int p = 0; p < DIMENSION; p++)
	{
		int Offset = abs(0.5*Circumferance_at_centre - p);
		double Radius_at_offset = sqrt((double)(DIMENSION*DIMENSION) - (double)(Offset*Offset));
		Circumferance_at_offset[p] = (int)(3.14159 * 2 * Radius_at_offset);

		for (int q = 0; q < Circumferance_at_centre; q++)
		{
			if ((double)p > 0.5*(double)(Circumferance_at_centre - Circumferance_at_offset[p]) && p > 0.5*(double)(Circumferance_at_centre + Circumferance_at_offset[p]))
			{
				Set_Number_Bonds_Glycan(-1, p, q);
				Set_Number_Bonds_Peptide(-1, p, q);
			}
		}
	}

}

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
			Number[i] = Number[i] + GLYCAN_SPRING_CONSTANT*abs(Murein[i][j].Return_Number_Bonds_Glycan());
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


		while (((q - length_of_current_bar_left) > 0) && ((Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Glycan() == 1) || (length_of_current_bar_left <= length_of_below_bar_left)))
		{
			if (Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() != 0)
			{
				if (Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Glycan() !=0 || Murein[p - scroll_up][q - length_of_current_bar_left + 1].Return_Number_Bonds_Glycan() == 1)
				{
					join_up_tally++;
				}
			}

			length_of_current_bar_left++;

		}

		while (((q + length_of_current_bar_right) < DIMENSION - 1) && (Murein[p - scroll_up][q + length_of_current_bar_right].Return_Number_Bonds_Glycan() !=0 || length_of_current_bar_right <= length_of_below_bar_right))
		{
			if (Murein[p - scroll_up][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() != 0)
			{
				join_up_tally++;
			}

			length_of_current_bar_right++;
		}


		if (p - scroll_up != 0)
		{
			while (((q - length_of_above_bar_left) > 0) && (Murein[p + 1 - scroll_up][q - length_of_above_bar_left].Return_Number_Bonds_Glycan() !=0 || length_of_above_bar_left <= length_of_current_bar_left))
			{

				if (Murein[p - scroll_up][q - length_of_above_bar_left].Return_Number_Bonds_Peptide() !=0)
				{
					join_down_tally++;
				}

				length_of_above_bar_left++;

			}

			while (((q + length_of_above_bar_right) < DIMENSION - 1) && (Murein[p - scroll_up][q + length_of_above_bar_right].Return_Number_Bonds_Glycan() !=0 || length_of_above_bar_right <= length_of_current_bar_right))
			{
				if (Murein[p - scroll_up][q + length_of_above_bar_right].Return_Number_Bonds_Peptide() != 0)
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

		////join_up_tally = join_up_tally - Murein[p - scroll_up][q].Return_Number_Bonds_Peptide();
		//join_down_tally = join_down_tally - Murein[p - scroll_up][q].Return_Number_Bonds_Peptide();

		Numerator = Numerator*join_up_tally;
		Denominator = Denominator*join_down_tally;

		length_of_below_bar_left = length_of_current_bar_left;
		length_of_below_bar_right = length_of_current_bar_right;

		length_of_current_bar_left = 0;
		length_of_current_bar_right = 0;

		length_of_above_bar_left = 0;
		length_of_above_bar_right = 0;



		Input_Force = Input_Force*Numerator / Denominator;

		if (Numerator == 0)
		{
			Input_Force == 0;
			return;
		}
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


	for (int scroll_down = 0; scroll_down <= DIMENSION - 1 - p; scroll_down++)
	{
		double Numerator = 1;
		double Denominator = 1;
		int join_up_tally = 0;
		int join_down_tally = 0;

		//glycans and q run horizontal, peptides and p run vertical, this set here

		while (q - length_of_current_bar_left > 0 && (Murein[p + scroll_down][q - length_of_current_bar_left].Return_Number_Bonds_Glycan() !=0 || length_of_current_bar_left <= length_of_above_bar_left))
		{

			if (Murein[p + scroll_down][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() !=0)
			{
				join_up_tally++;
			}

			length_of_current_bar_left++;

		}

		while (q + length_of_current_bar_right < DIMENSION - 1 && (Murein[p + scroll_down][q + length_of_current_bar_right].Return_Number_Bonds_Glycan() !=0 || length_of_current_bar_right <= length_of_above_bar_right))
		{
			if (Murein[p + scroll_down][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() !=0)
			{
				join_up_tally++;
			}

			length_of_current_bar_right++;
		}

		if (p + scroll_down != DIMENSION - 1)
		{
			while (q - length_of_below_bar_left > 0 && (Murein[p + 1 + scroll_down][q - length_of_below_bar_left].Return_Number_Bonds_Glycan() !=0 || length_of_below_bar_left <= length_of_current_bar_left))
			{

				if (Murein[p + scroll_down][q - length_of_below_bar_left].Return_Number_Bonds_Peptide() !=0)
				{
					join_down_tally++;
				}

				length_of_below_bar_left++;

			}

			while (q + length_of_below_bar_right < DIMENSION - 1 && (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() !=0 || length_of_below_bar_right <= length_of_current_bar_right))
			{
				if (Murein[p + scroll_down][q + length_of_below_bar_right].Return_Number_Bonds_Peptide() !=0)
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

		//join_up_tally = join_up_tally - Murein[p + scroll_down][q].Return_Number_Bonds_Peptide();
		//join_down_tally = join_down_tally - Murein[p + scroll_down][q].Return_Number_Bonds_Peptide();

		Numerator = Numerator*join_up_tally;
		Denominator = Denominator*join_down_tally;

		length_of_above_bar_left = length_of_current_bar_left;
		length_of_above_bar_right = length_of_current_bar_right;

		length_of_current_bar_left = 0;
		length_of_current_bar_right = 0;

		length_of_below_bar_left = 0;
		length_of_below_bar_right = 0;

		Input_Force = Input_Force*Numerator / Denominator;

		if (Numerator == 0)
		{
			Input_Force == 0;
			return;
		}
	}

	Force_Downwards[p][q] = Input_Force;
	return;

}

void Polymer::Find_Force_Leftwards(double Input_Force, int q, int p)
{

	int length_of_current_bar_up = 0;
	int length_of_current_bar_down = 0;
	int length_of_above_bar_up = 0;
	int length_of_above_bar_down = 0;
	int length_of_below_bar_up = 0;
	int length_of_below_bar_down = 0;

	if (Murein[p][q].Return_Number_Bonds_Glycan() == 0)
	{
		Force_Leftwards[p][q] = 0;
		return;
	}

	for (int scroll_left = 0; scroll_left <= p; scroll_left++)
	{
		double Numerator = 1;
		double Denominator = 1;
		int join_left_tally = 0;
		int join_right_tally = 0;


		while (((q - length_of_current_bar_up) > 0) && ((Murein[p - scroll_left][q - length_of_current_bar_up].Return_Number_Bonds_Peptide() == 1) || (length_of_current_bar_up <= length_of_below_bar_up)))
		{
			if (Murein[p - scroll_left][q - length_of_current_bar_up].Return_Number_Bonds_Glycan() == 1)
			{
				if (Murein[p - scroll_left][q - length_of_current_bar_up].Return_Number_Bonds_Peptide() == 1 || Murein[p - scroll_left][q - length_of_current_bar_up + 1].Return_Number_Bonds_Peptide() == 1)
				{
					join_left_tally++;
				}
			}

			length_of_current_bar_up++;

		}

		while (((q + length_of_current_bar_down) < DIMENSION - 1) && (Murein[p - scroll_left][q + length_of_current_bar_down].Return_Number_Bonds_Peptide() == 1 || length_of_current_bar_down <= length_of_below_bar_down)) //is it here???
		{
			if (Murein[p - scroll_left][q + length_of_current_bar_down].Return_Number_Bonds_Glycan() == 1)
			{
				join_left_tally++;
			}

			length_of_current_bar_down++;
		}


		if (p - scroll_left != 0)
		{
			while (((q - length_of_above_bar_up) > 0) && (Murein[p + 1 - scroll_left][q - length_of_current_bar_up].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_up <= length_of_current_bar_up))
			{

				if (Murein[p - scroll_left][q - length_of_above_bar_up].Return_Number_Bonds_Glycan() == 1)
				{
					join_right_tally++;
				}

				length_of_above_bar_up++;

			}

			while (((q + length_of_above_bar_down) < DIMENSION - 1) && (Murein[p - scroll_left][q + length_of_current_bar_down].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_down <= length_of_current_bar_down))
			{
				if (Murein[p - scroll_left][q + length_of_above_bar_down].Return_Number_Bonds_Glycan() == 1)
				{
					join_right_tally++;
				}

				length_of_above_bar_down++;
			}
		}
		else
		{
			join_right_tally = DIMENSION; //only true for Peptides, for Glycans it loops round
		}

		//join_left_tally = join_left_tally - Murein[p - scroll_left][q].Return_Number_Bonds_Glycan();
		//join_right_tally = join_right_tally - Murein[p - scroll_left][q].Return_Number_Bonds_Glycan();

		Numerator = Numerator*join_left_tally;
		Denominator = Denominator*join_right_tally;

		length_of_below_bar_up = length_of_current_bar_up;
		length_of_below_bar_down = length_of_current_bar_down;

		length_of_current_bar_up = 0;
		length_of_current_bar_down = 0;

		length_of_above_bar_up = 0;
		length_of_above_bar_down = 0;



		Input_Force = Input_Force*Numerator / Denominator;

		if (Numerator == 0)
		{
			Input_Force == 0;
			return;
		}
	}

	Force_Leftwards[p][q] = Input_Force;
	return;

}

void Polymer::Find_Force_Rightwards(double Input_Force, int q, int p)
{

	int length_of_current_bar_up = 0;
	int length_of_current_bar_down = 0;
	int length_of_above_bar_up = 0;
	int length_of_above_bar_down = 0;
	int length_of_below_bar_up = 0;
	int length_of_below_bar_down = 0;

	if (Murein[p][q].Return_Number_Bonds_Glycan() == 0)
	{
		Force_Rightwards[p][q] = 0;
		return;
	}

	for (int scroll_right = 0; scroll_right <= DIMENSION - 1 - p; scroll_right++)
	{
		double Numerator = 1;
		double Denominator = 1;
		int join_left_tally = 0;
		int join_right_tally = 0;


		while (((q - length_of_current_bar_up) > 0) && ((Murein[p + scroll_right][q].Return_Number_Bonds_Peptide() == 1) || (length_of_current_bar_up <= length_of_below_bar_up)))
		{
			if (Murein[p + scroll_right][q - length_of_current_bar_up].Return_Number_Bonds_Glycan() == 1)
			{
				if (Murein[p + scroll_right][q - length_of_current_bar_up].Return_Number_Bonds_Peptide() == 1 || Murein[p + scroll_right][q - length_of_current_bar_up + 1].Return_Number_Bonds_Peptide() == 1)
				{
					join_left_tally++;
				}
			}

			length_of_current_bar_up++;

		}

		while (((q + length_of_current_bar_down) < DIMENSION - 1) && (Murein[p + scroll_right][q + length_of_current_bar_down].Return_Number_Bonds_Peptide() == 1 || length_of_current_bar_down <= length_of_below_bar_down))
		{
			if (Murein[p + scroll_right][q + length_of_current_bar_down].Return_Number_Bonds_Glycan() == 1)
			{
				join_left_tally++;
			}

			length_of_current_bar_down++;
		}


		if (p + scroll_right != 0)
		{
			while (((q - length_of_above_bar_up) > 0) && (Murein[p + 1 + scroll_right][q - length_of_above_bar_up].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_up <= length_of_current_bar_up))
			{

				if (Murein[p + scroll_right][q - length_of_above_bar_up].Return_Number_Bonds_Glycan() == 1)
				{
					join_right_tally++;
				}

				length_of_above_bar_up++;

			}

			while (((q + length_of_above_bar_down) < DIMENSION - 1) && (Murein[p + scroll_right][q].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_down <= length_of_current_bar_down))
			{
				if (Murein[p + scroll_right][q + length_of_above_bar_down].Return_Number_Bonds_Glycan() == 1)
				{
					join_right_tally++;
				}

				length_of_above_bar_down++;
			}
		}
		else
		{
			join_right_tally = DIMENSION; //only true for Peptides, for Glycans it loops round
		}

		//join_left_tally = join_left_tally - Murein[p + scroll_right][q].Return_Number_Bonds_Glycan();
		//join_right_tally = join_right_tally - Murein[p + scroll_right][q].Return_Number_Bonds_Glycan();

		Numerator = Numerator*join_left_tally;
		Denominator = Denominator*join_right_tally;

		length_of_below_bar_up = length_of_current_bar_up;
		length_of_below_bar_down = length_of_current_bar_down;

		length_of_current_bar_up = 0;
		length_of_current_bar_down = 0;

		length_of_above_bar_up = 0;
		length_of_above_bar_down = 0;



		Input_Force = Input_Force*Numerator / Denominator;

		if (Numerator == 0)
		{
			Input_Force == 0;
			return;
		}
	}

	Force_Rightwards[p][q] = Input_Force;
	return;

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
	Min_Length_Glycan = RAND_MAX;
	Min_Length_Peptide = RAND_MAX;

	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			cout << i << "\t" << j << "\t" << Murein[i][j].Return_Length_Peptide() << "\t" << Murein[i][j].Return_Vertical_Force() << "\t" << Murein[i][j].Return_Length_Glycan() << "\t" << Murein[i][j].Return_Horizontal_Force() << endl;

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