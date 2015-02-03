//////we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
//////first we need to find what length levels p-1 and p-2 are
////
////void Polymer::Find_Force_Upwards(int p, int q)
////{
////	
////	int length_of_current_bar_left = 0;
////	int length_of_current_bar_right = 0;
////	int length_of_above_bar_left = 0;
////	int length_of_above_bar_right = 0;
////	int length_of_below_bar_left = 0;
////	int length_of_below_bar_right = 0;
////
////	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
////	{
////		Force_Upwards[p][q] = 0;
////		return;
////	}
////
////
////	for (int scroll_up = 0; scroll_up <= p; scroll_up++)
////	{
////		double Numerator = 1;
////		double Denominator = 1;
////		int join_up_tally = 0;
////		int join_down_tally = 0;
////
////		//glycans and q run horizontal, peptides and p run vertical, this set here
////
////		while (q - length_of_current_bar_left >= 0 || (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_left < length_of_below_bar_left));
////		{
////
////			if (Murein[p - scroll_up][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() == 1)
////			{
////				join_up_tally++;
////			}
////
////			length_of_current_bar_left++;
////
////		}
////
////		while (q + length_of_current_bar_right <= 100 || (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_right < length_of_below_bar_right));
////		{
////			if (Murein[p - scroll_up][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() == 1)
////			{
////				join_up_tally++;
////			}
////
////			length_of_current_bar_right++;
////		}
////
////		if (p - scroll_up != 0)
////		{
////			while (q - length_of_above_bar_left >= 0 || (Murein[p + 1 - scroll_up][q].Return_Number_Bonds_Glycan() == 1 && length_of_above_bar_left < length_of_current_bar_left));
////			{
////
////				if (Murein[p - scroll_up][q - length_of_above_bar_left].Return_Number_Bonds_Peptide() == 1)
////				{
////					join_down_tally++;
////				}
////
////				length_of_above_bar_left++;
////
////			}
////
////			while (q + length_of_above_bar_right <= 100 || (Murein[p - scroll_up][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_right < length_of_below_bar_right));
////			{
////				if (Murein[p - scroll_up][q + length_of_above_bar_right].Return_Number_Bonds_Peptide() == 1)
////				{
////					join_down_tally++;
////				}
////
////				length_of_above_bar_right++;
////			}
////		}
////		else
////		{
////			join_down_tally = 100; //only true for glycans, for peptides it loops round
////		}
////
////		Numerator = Numerator*join_up_tally;
////		Denominator = Denominator*join_down_tally;
////
////		length_of_below_bar_left = length_of_current_bar_left;
////		length_of_below_bar_right = length_of_current_bar_right;
////
////		length_of_current_bar_left = 0;
////		length_of_current_bar_right = 0;
////
////		length_of_below_bar_left = 0;
////		length_of_below_bar_right = 0;
////
////		Input_Force = Input_Force*Numerator / Denominator;
////	}
////
////	Force_Upwards[p][q] = Input_Force;
////	return;
////
////}
////
////void Polymer::Find_Force_Downwards(int p, int q)
////{
////
////	int length_of_current_bar_left = 0;
////	int length_of_current_bar_right = 0;
////	int length_of_below_bar_left = 0;
////	int length_of_below_bar_right = 0;
////	int length_of_above_bar_left = 0;
////	int length_of_above_bar_right = 0;
////
////	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
////	{
////		Force_Downwards[p][q] = 0;
////		return;
////	}
////
////
////	for (int scroll_down = 0; scroll_down <= 100-p; scroll_down++)
////	{
////		double Numerator = 1;
////		double Denominator = 1;
////		int join_up_tally = 0;
////		int join_down_tally = 0;
////
////		//glycans and q run horizontal, peptides and p run vertical, this set here
////
////		while (q - length_of_current_bar_left >= 0 || (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_left < length_of_above_bar_left));
////		{
////
////			if (Murein[p + scroll_down][q - length_of_current_bar_left].Return_Number_Bonds_Peptide() == 1)
////			{
////				join_up_tally++;
////			}
////
////			length_of_current_bar_left++;
////
////		}
////
////		while (q + length_of_current_bar_right <= 100 || (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_right < length_of_above_bar_right));
////		{
////			if (Murein[p + scroll_down][q + length_of_current_bar_right].Return_Number_Bonds_Peptide() == 1)
////			{
////				join_up_tally++;
////			}
////
////			length_of_current_bar_right++;
////		}
////
////		if (p + scroll_down != 100)
////		{
////			while (q - length_of_below_bar_left >= 0 || (Murein[p + 1 + scroll_down][q].Return_Number_Bonds_Glycan() == 1 && length_of_below_bar_left < length_of_current_bar_left));
////			{
////
////				if (Murein[p + scroll_down][q - length_of_below_bar_left].Return_Number_Bonds_Peptide() == 1)
////				{
////					join_down_tally++;
////				}
////
////				length_of_below_bar_left++;
////
////			}
////
////			while (q + length_of_below_bar_right <= 100 || (Murein[p + scroll_down][q].Return_Number_Bonds_Glycan() == 1 && length_of_current_bar_right < length_of_above_bar_right));
////			{
////				if (Murein[p + scroll_down][q + length_of_below_bar_right].Return_Number_Bonds_Peptide() == 1)
////				{
////					join_down_tally++;
////				}
////
////				length_of_below_bar_right++;
////			}
////		}
////		else
////		{
////			join_down_tally = 100; //only true for glycans, for peptides it loops round
////		}
////
////		Numerator = Numerator*join_up_tally;
////		Denominator = Denominator*join_down_tally;
////
////		length_of_above_bar_left = length_of_current_bar_left;
////		length_of_above_bar_right = length_of_current_bar_right;
////
////		length_of_current_bar_left = 0;
////		length_of_current_bar_right = 0;
////
////		length_of_above_bar_left = 0;
////		length_of_above_bar_right = 0;
////
////		Input_Force = Input_Force*Numerator / Denominator;
////	}
////
////	Force_Downwards[p][q] = Input_Force;
////	return;
////
////}
////
////void Polymer::Find_Force_Upwards(double Input_Force, int p, int q)
////{
////
////	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
////	//first we need to find what length levels p-1 and p-2 are
////
////	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
////	{
////		Force_Upwards[p][q] = 0;
////		return;
////	}
////
////	int leftward_extent_above[DIMENSION]; //n can take values less that p and represents the level
////	int rightward_extent_above[DIMENSION]; //n can take values less that p and represents the level
////
////	int total_left = 0; // this total represents the current length of all the bars being taken into account
////
////	for (int n = p - 1; n > 0; n--)
////	{
////		int m = q;
////
////		leftward_extent_above[n] = 0;
////
////		/*while ((Murein[n][m].Return_Number_Bonds_Glycan() == 1 || leftward_extent_above[n] < total_left) && m <= DIMENSION - 1 && m >= 0)
////		{
////			leftward_extent_above[n]++;
////			m--;
////		}*/
////
////		total_left = leftward_extent_above[n];
////
////
////	}
////
////	int total_right = 0; // this total represents the current length of all the bars being taken into account
////
////	for (int n = p - 1; n > 0; n--)
////	{
////
////		rightward_extent_above[n] = 0;
////		int m = q;
////
////		while ((Murein[n][m].Return_Number_Bonds_Glycan() == 1 || rightward_extent_above[n] < total_right) && m <= DIMENSION - 1 && m >= 0)
////		{
////			rightward_extent_above[n]++;
////			m++;
////		}
////
////		total_right = rightward_extent_above[n];
////
////
////	}
////
////	if (p == 0)
////	{
////		leftward_extent_above[p] = q;
////		rightward_extent_above[p] = DIMENSION - 1 - q;
////	}
////
////
////	//this finds the lengths of all the above levels
////
////	//now the numerator is the number of joins which connects anything within these levels to the one above it
////
////	Numerator[p][q] = 1;
////
////	for (int n = p - 1; n >= 0; n--)
////	{
////		int m = q;
////		int tally_peptides = 0;
////
////		for (int i = 0; i < leftward_extent_above[n]; i++)
////		{
////			if (Murein[n][m - i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		for (int i = 0; i < rightward_extent_above[n]; i++)
////		{
////			if (Murein[n][m + i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////		}
////
////		tally_peptides = tally_peptides - Murein[n][q].Return_Number_Bonds_Peptide();
////
////		if (n - 1 >= DIMENSION - 1 || n - 1 <= 0)
////		{
////			tally_peptides = DIMENSION;
////		}
////
////
////		Numerator[p][q] = Numerator[p][q] * (double)tally_peptides;
////
////	}
////
////	Denominator[p][q] = 1;
////
////	for (int n = p - 1; n >= 0; n--)
////	{
////		int m = q;
////		int tally_peptides = 0;
////
////		for (int i = 0; i < leftward_extent_above[n]; i++)
////		{
////			if (Murein[n - 1][m - i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		for (int i = 0; i < rightward_extent_above[n]; i++)
////		{
////			if (Murein[n - 1][m + i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////		}
////
////		tally_peptides = tally_peptides - Murein[n][q].Return_Number_Bonds_Peptide();
////
////		if (n + 1 >= DIMENSION - 1 || n <= 0)
////		{
////			tally_peptides = DIMENSION;
////		}
////
////		Denominator[p][q] = Denominator[p][q] * (double)tally_peptides;
////
////	}
////
////	int m = q;
////
////	int leftward_extent = 0;
////
////	while (Murein[p][m].Return_Number_Bonds_Glycan() == 1 && m <= DIMENSION - 1 && m >= 0)
////	{
////		leftward_extent++;
////		m--;
////	}
////
////	int rightward_extent = 0;
////	m = q;
////
////	while (Murein[p][m].Return_Number_Bonds_Glycan() == 1 && m <= DIMENSION - 1 && m >= 0)
////	{
////		rightward_extent++;
////		m++;
////	}
////
////	m = q;
////	int tally_peptides = 0;
////
////	for (int i = 0; i < leftward_extent; i++)
////	{
////		if (Murein[p][m - i].Return_Number_Bonds_Peptide() == 1)
////		{
////			tally_peptides++;
////		}
////
////	}
////
////	for (int i = 0; i < rightward_extent; i++)
////	{
////		if (Murein[p][m + i].Return_Number_Bonds_Peptide() == 1)
////		{
////			tally_peptides++;
////		}
////	}
////
////	tally_peptides = tally_peptides - Murein[p][q].Return_Number_Bonds_Peptide();
////
////	double Frac = (double)Murein[p][q].Return_Number_Bonds_Peptide() / (double)tally_peptides;
////
////	Force_Upwards[p][q] = Frac*Input_Force*Numerator[p][q] / Denominator[p][q];
////
////}
////
////void Polymer::Find_Force_Downwards(double Input_Force, int p, int q)
////{
////
////	if (Murein[p][q].Return_Number_Bonds_Peptide() == 0)
////	{
////		Force_Downwards[p][q] = 0;
////		return;
////	}
////	//we need to find the number of bonds joining level p-1 and p-2, then p-2 and p-3... 0 and 1
////	//first we need to find what length levels p-1 and p-2 are
////
////	int leftward_extent_below[DIMENSION]; //n can take values less that p and represents the level
////	int rightward_extent_below[DIMENSION]; //n can take values less that p and represents the level
////
////	int total_left = 0; // this total represents the current length of all the bars being taken into account
////
////	for (int n = p + 1; n < DIMENSION; n++)
////	{
////		int m = q;
////
////		leftward_extent_below[n] = 0;
////
////		while ((Murein[n][m].Return_Number_Bonds_Glycan() == 1 || leftward_extent_below[n] < total_left) && m <= DIMENSION - 1 && m >= 0)
////		{
////			leftward_extent_below[n]++;
////			m--;
////		}
////
////		total_left = leftward_extent_below[n];
////
////	}
////
////	int total_right = 0; // this total represents the current length of all the bars being taken into account
////
////	for (int n = p + 1; n < DIMENSION; n++)
////	{
////
////		rightward_extent_below[n] = 0;
////		int m = q;
////		int end = 1;
////
////		while ((Murein[n][m].Return_Number_Bonds_Glycan() || rightward_extent_below[n] < total_right) && m <= DIMENSION - 1 && m >= 0)
////		{
////			rightward_extent_below[n]++;
////			m++;
////		}
////
////		total_right = rightward_extent_below[n];
////
////
////	}
////
////	if (p == DIMENSION - 1)
////	{
////		leftward_extent_below[p] = q;
////		rightward_extent_below[p] = DIMENSION - 1 - q;
////	}
////
////
////	//this finds the lengths of all the above levels
////
////	//now the numerator is the number of joins which connects anything within these levels to the one above it
////
////	Numerator[p][q] = 1;
////
////	for (int n = p + 1; n < DIMENSION; n++)
////	{
////		int m = q;
////		int tally_peptides = 0;
////
////		for (int i = 0; i < leftward_extent_below[n]; i++)
////		{
////			if (Murein[n - 1][m - i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		for (int i = 0; i < rightward_extent_below[n]; i++)
////		{
////			if (Return_Number_Bonds_Peptide(n, m + i) == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		tally_peptides = tally_peptides - Murein[n][q].Return_Number_Bonds_Peptide();
////
////		if (n - 1 >= DIMENSION - 1 || n - 1 <= 0)
////		{
////			tally_peptides = DIMENSION;
////		}
////
////		Numerator[p][q] = Numerator[p][q] * (double)tally_peptides;
////
////	}
////
////	Denominator[p][q] = 1;
////
////	for (int n = p + 1; n < DIMENSION; n++)
////	{
////		int m = q;
////		int tally_peptides = 0;
////
////		for (int i = 0; i < leftward_extent_below[n]; i++)
////		{
////			if (Murein[n - 1][m - i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		for (int i = 0; i < rightward_extent_below[n]; i++)
////		{
////			if (Murein[n - 1][m + i].Return_Number_Bonds_Peptide() == 1)
////			{
////				tally_peptides++;
////			}
////
////		}
////
////		tally_peptides = tally_peptides - Murein[n][q].Return_Number_Bonds_Peptide();
////
////		if (n + 1 >= DIMENSION - 1 || n <= 0)
////		{
////			tally_peptides = DIMENSION;
////		}
////
////		Denominator[p][q] = Denominator[p][q] * (double)tally_peptides;
////
////	}
////
////	Force_Downwards[p][q] = Input_Force*Numerator[p][q] / Denominator[p][q];
////
////}
//
//void Polymer::Find_Force_Rightwards(double Input_Force, int q, int p)
//{
//
//	int length_of_current_bar_up = 0;
//	int length_of_current_bar_down = 0;
//	int length_of_above_bar_up = 0;
//	int length_of_above_bar_down = 0;
//	int length_of_below_bar_up = 0;
//	int length_of_below_bar_down = 0;
//
//	if (Murein[p][q].Return_Number_Bonds_Glycan() == 0)
//	{
//		Force_Rightwards[p][q] = 0;
//		return;
//	}
//
//	for (int scroll_right = 0; scroll_right <= DIMENSION-1-p; scroll_right++)
//	{
//		double Numerator = 1;
//		double Denominator = 1;
//		int join_left_tally = 0;
//		int join_right_tally = 0;
//
//		
//		while (((q - length_of_current_bar_up) > 0) && ((Murein[p +scroll_right][q].Return_Number_Bonds_Peptide() == 1) || (length_of_current_bar_up <= length_of_below_bar_up)))
//		{
//			if (Murein[p +scroll_right][q - length_of_current_bar_up].Return_Number_Bonds_Glycan() == 1)
//			{
//				if (Murein[p +scroll_right][q - length_of_current_bar_up].Return_Number_Bonds_Peptide() == 1 || Murein[p +scroll_right][q - length_of_current_bar_up + 1].Return_Number_Bonds_Peptide() == 1)
//				{
//					join_left_tally++;
//				}
//			}
//
//			length_of_current_bar_up++;
//
//		}
//
//		while (((q + length_of_current_bar_down) < DIMENSION - 1) && (Murein[p +scroll_right][q].Return_Number_Bonds_Peptide() == 1 || length_of_current_bar_down <= length_of_below_bar_down))
//		{
//			if (Murein[p +scroll_right][q + length_of_current_bar_down].Return_Number_Bonds_Glycan() == 1)
//			{
//				join_left_tally++;
//			}
//
//			length_of_current_bar_down++;
//		}
//
//
//		if (p +scroll_right != 0)
//		{
//			while (((q - length_of_above_bar_up) > 0) && (Murein[p + 1 +scroll_right][q].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_up <= length_of_current_bar_up))
//			{
//
//				if (Murein[p +scroll_right][q - length_of_above_bar_up].Return_Number_Bonds_Glycan() == 1)
//				{
//					join_right_tally++;
//				}
//
//				length_of_above_bar_up++;
//
//			}
//
//			while (((q + length_of_above_bar_down) < DIMENSION - 1) && (Murein[p +scroll_right][q].Return_Number_Bonds_Peptide() == 1 || length_of_above_bar_down <= length_of_current_bar_down))
//			{
//				if (Murein[p +scroll_right][q + length_of_above_bar_down].Return_Number_Bonds_Glycan() == 1)
//				{
//					join_right_tally++;
//				}
//
//				length_of_above_bar_down++;
//			}
//		}
//		else
//		{
//			join_right_tally = DIMENSION; //only true for Peptides, for Glycans it loops round
//		}
//
//		join_left_tally = join_left_tally - Murein[p +scroll_right][q].Return_Number_Bonds_Glycan();
//		join_right_tally = join_right_tally - Murein[p +scroll_right][q].Return_Number_Bonds_Glycan();
//
//		Numerator = Numerator*join_left_tally;
//		Denominator = Denominator*join_right_tally;
//
//		length_of_below_bar_up = length_of_current_bar_up;
//		length_of_below_bar_down = length_of_current_bar_down;
//
//		length_of_current_bar_up = 0;
//		length_of_current_bar_down = 0;
//
//		length_of_above_bar_up = 0;
//		length_of_above_bar_down = 0;
//
//
//
//		Input_Force = Input_Force*Numerator / Denominator;
//	}
//
//	Force_Rightwards[p][q] = Input_Force;
//	return;
//
//}
